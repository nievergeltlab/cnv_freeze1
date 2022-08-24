
#Make a place to store temp data
 mkdir temp

###Phenotype and covar file creation (all studies should go in here)

#Phenotype data extraction:
  
 #Affy: Reformat UKB phenotype to the same format as other studies
  cat /mnt/ukbb/adam/ptsd/UKB_ptsd_eur_unrelated_m0_may13_2021.pheno | awk '{print $1,$2,$3,$4,$5,$6,$7}' > phenotype_files/p2_ukbb_eur_ukbb.cov
  cat /mnt/ukbb/adam/ptsd/UKB_ptsd_eur_unrelated_m0_may13_2021.pheno | awk '{print $1,$2,$14}' > phenotype_files/p2_ukbb_ukbb_eur.pheno

 #FOR ALL OTHERS: Phenotype files are already from freeze 3 - an update from f2.5, specifically some studies are now case/control (PTS1, PSY3 Feeny), some are now QT (psy3 NCMH)

   #At some point, convert my IDs to Marieke's, or vice versa..
   
#Append study name to the phenotype file
 for files in $(ls phenotype_files | grep .pheno)
 {
  filename=$(echo $files | sed 's/.pheno//g' | sed 's/p2_//g' | sed 's/_eur//g' )
  awk -v filename=$filename '{print $0,filename}' phenotype_files/$files > temp/$files
 }
 cat temp/*.pheno |  awk '{if(NR == 1) {$4 = "studyname"}; if (NR==1 || ($1 != "FID" && $3 != "NA")) print $1,$2,$3,$4}' > all_pgc2.subjects
 cat phenotype_files/*.cov |  awk '{if (NR==1 || ($1 != "FID" &&  $3 != "NA")) print $0}' > all_pgc2.cov 

###Filter the Jeff datasets, combine with phenotypes, count Ns

R

#Libraries
 library(data.table)
 library(plyr)

#Load phenotypes
 pheno <- fread('all_pgc2.subjects',data.table=F)
  pheno$pheno <- pheno[,3]
  
  #11/8/21 - Made the decision to combine bry2, niut, eacr cohorts, feen, run as case/contorl (all data already use same PCs). 
            #CVC will be run with MRS 1 and 2 datasets (all data already use same PCs).
  pheno$studyname2 <- pheno$studyname
  pheno[which(pheno$studyname %in% c("psy3_bry2", "psy3_eacr", "psy3_feen",  "psy3_niut")),]$studyname2 <- "psy3_comb"
   pheno[which(pheno$studyname %in% c("mrsc_cvc", "mrsc_mrsc")),]$studyname2 <- "mrsc_comb"

#Load PC covariates
 covar <- fread('all_pgc2.cov',data.table=F)

#Load Jeff manifest
 d1a <- read.csv('PTSD_manifest_20210403_am.csv',header=T,stringsAsFactors=F,na.strings=c("-9",NA,"#N/A"))
 #jeff uses my ids
  d1a$FID <- d1a$rpfid
  d1a$IID <- d1a$rpiid
 
#the ids in here are hers, so I set them up to be changed
 marieke <- read.csv('UKBB_manifest_20211116.csv',header=T,stringsAsFactors=F,na.strings=c("-9",NA,"#N/A"))
  marieke$rpfid <- marieke$IID #FID has some crap on it
  marieke$rpiid <- marieke$IID
  #make some space for my ids
  marieke$FID <- NULL
  marieke$IID <- NULL 
 marieke <- subset(marieke, Rare_outlier == 0 )
 marieke <- marieke[,-c("V21","X.batch.1", "identifier.1", "dataset.1", "Case.Control.1" ,"Freq_Calc.1" ,"PLATFORM.1")]
 marieke$NOTES..additional.filter.reason.IQR <- marieke$NOTES..additional.filter.reason.
 
#marieke id update
 marieke_ids <- read.csv('ukbb_cnv_idconversion_dec6_2021.csv',header=T,stringsAsFactors=F,na.strings=c("-9",NA,"#N/A"))
 marieke2 <- merge(marieke,marieke_ids,by=c("rpfid","rpiid"))
 
 #Add in PCL dx, to be compatible with other data
 pcl <- fread('UKBB_pcl_sum_dec6_2021.csv',data.table=F)
 names(pcl)[3] <- "Diagnosis.Phenotype"
 marieke3 <- merge(marieke2,pcl,by=c("FID","IID"))

write.csv(subset(marieke3,select=c('rpiid','Diagnosis.Phenotype')),file="ptsd_status_mariekeids.csv",quote=F,row.names=F)# this is for omar from august 5 2022)
 
#Append UKBB exclusion info to jeffs dataset.
#Load Marieke's CNV IDs. Load my UKBB that passed QC. 
 d1a2 <- rbind.fill(d1a,marieke3) #be sure to add 'dataset' and studyname to the file as well as NOTES..additional.filter.reason.IQR and diagnosis.phenotype 
  
#Which studies are being exclued?
 table(d1a2$dataset,d1a2$NOTES..additional.filter.reason.IQR)
 
#it is fine if it is missing case/control, as long as it has my phenotyping.
 d1b <- subset(d1a2,NOTES..additional.filter.reason.IQR %in% c("pass") )

#Add PCs. European only
 d1m1 <- merge(d1b,covar,by=c("FID","IID"))
 
#Add phenotypes
 d1.1 <- merge(d1m1,pheno,by=c("FID","IID"))


#Subset to only phenotyped. I am removing acouple of people with unknown case/control status - because I will do case/control analysis of all of these, I'd rather have overlapping samples
 d1 <- subset(d1.1,!is.na(pheno) & pheno != "-9"  & Diagnosis.Phenotype != "Unknown" & gender!="U")

#Make dx a simpler variable
 d1$dx <- d1[,"Diagnosis.Phenotype"]
 table(d1$studyname2,d1$dx)

#note sample attrition 
#Sample attrition notes:
#A substantial number of MRS1 subjects failed sample level QC
#GMRF numbers are indeed around 200 - that is what is in freeze 2
#Updated manifest to reflect the fact that the freeze 3 pts1 IDs 

 table(d1a$dataset,d1a$NOTES..additional.filter.reason.IQR)
 table(d1b$dataset)
 table(d1m1$dataset)
 table(d1.1$dataset)
 table(d1$dataset)
 table(d1$studyname)
 table(d1$studyname2)


#This sample list must be merged with the CNV calls:
 #cat cnv_calls/*.txt | awk '{if (NR==1 || $1 != "chrom") print}' > ptsd_10kb_10probes.txt


#Note: CNV call files only note the presence of CNV, ergo not everyone has calls
#Note: Contains more samples than our PTSD data, e.g. sample Mayo-21
 cnvtypes <- fread('ptsd_10kb_10probes.txt',data.table=F)
 
 #Figure out which samples have calls
 d1$found <- 0
 d1[which(d1$identifier %in% cnvtypes$sampleId ),]$found <- 1
 table(d1$studyname,d1$found)
 
#Make a .cnv file
 d1s <- subset(d1,select=c(identifier,FID,IID))
 names(d1s)[1] <- "sampleId"
 cnvtypes_exp2 <- merge(d1s,cnvtypes,by="sampleId") #all.x is unnecesary, because the .fam file is the one that includes subjects with no CNV info. It would make no sense for hte no calls subjects to be in here, because it only lsits when calls are found
 length(unique(cnvtypes_exp2$sampleId)) #this should match n found
 
 cnvtypes_exp2$FID <- cnvtypes_exp2$FID
 cnvtypes_exp2$IID <- cnvtypes_exp2$IID
 cnvtypes_exp2$CHR <- gsub("chr","",cnvtypes_exp2$chrom,)
 cnvtypes_exp2$BP1 <- cnvtypes_exp2$locusStart
 cnvtypes_exp2$BP2 <- cnvtypes_exp2$locusEnd
 cnvtypes_exp2$TYPE <- gsub("DEL","1",cnvtypes_exp2$locusType)
 cnvtypes_exp2$TYPE <- gsub("DUP","3",cnvtypes_exp2$TYPE) #Some people have a dup and a del! they have to be excluded..
 #dup/del exclusion code here
 cnvtypes_exp2$SCORE <- cnvtypes_exp2$Genes.number
 cnvtypes_exp2$SITES <- cnvtypes_exp2$numProbes
 #Just get rid of the deldups
 cnvtypes_exp2 <- subset(cnvtypes_exp2,TYPE != "13")
 
 #CNV files have the following format:   FID	IID	CHR	BP1	BP2	TYPE	SCORE	SITES 
  #I use number of genes as score (but this should beconfidence score)
  #For type, del is coded to 1, dup is 3
 cnvtypes_exp_jeff <- subset(cnvtypes_exp2, select=c(FID,IID,CHR,BP1,BP2,TYPE,SCORE,SITES))
 write.table(cnvtypes_exp_jeff,file="ptsd_10kb_10probes.cnv",quote=F,row.names=F,col.names=T)
 

#Make a .fam file with all subjects
 d1fam <- subset(d1,select=c(FID,IID,gender,dx))
 d1fam$FATHER <- 0
 d1fam$MOTHER <- 0
 d1fam$sex <- gsub("F","2", d1fam$gender)
 d1fam$sex <- gsub("M","1", d1fam$sex)
 d1fam$sex <- gsub("U","0", d1fam$sex)
 d1fam$phenoPTS <- gsub("PTSD",2,d1fam$dx)
 d1fam$phenoPTS <- gsub("CONT",1,d1fam$phenoPTS)
 d1fam_exp <- subset(d1fam,select=c(FID,IID,FATHER,MOTHER,sex,phenoPTS))
 write.table(d1fam_exp,file="ptsd_10kb_10probes.fam",quote=F,row.names=F,col.names=F)

#export study names as well
 write.table(subset(d1,select=c("FID","IID","studyname2")),file="ptsd_10kb_10probes.studies",quote=F,row.names=F,col.names=F)


#Make total phenotype and covariate files

# cat conversion_keys/*.txt | awk '{if (NR == 1 || $1 != "filename") print}' > all_conversion_keys.txt
# cat pcv_covariates/*.txt | awk '{if (NR == 1 || $1 != "File") print}' > all_pgc2.qualitycov

 conversion_keys <- read.table('all_conversion_keys.txt',header=T,stringsAsFactors=F,na.strings=c("-9",NA,"#N/A"))
  names(conversion_keys)[2:3] <- c("FID","IID")
 pc_covars <- read.table('all_pgc2.cov',header=T,stringsAsFactors=F,na.strings=c("-9",NA,"#N/A"))
  
 quality_covars <- read.table('all_pgc2.qualitycov',header=T,stringsAsFactors=F,na.strings=c("-9",NA,"#N/A"))
 quality_covars$subjectID <- quality_covars$File
 
#Merge these data
 cnv_cov1 <- merge(quality_covars,conversion_keys,by=c("subjectID"))
 
 # cnv_cov2 <- merge(cnv_cov1,pc_covars,by=c("FID","IID")) #pcs were already added, no need to add these
 cnv_cov   <- merge(cnv_cov1,d1,by=c("FID","IID"),all.y=TRUE)
 cnv_cov$FID_IID <- paste(cnv_cov$FID,cnv_cov$IID,sep="_")
 
 #remove duplicates, according to which one qced worse (because that sample had to have been the replicate that was removed)
 cnv_cov <- cnv_cov[order(cnv_cov$LRR_SD),]
 cnv_cov <- cnv_cov[!duplicated(cnv_cov$FID_IID),]
 
 #A few UKBB samples are missing due to Sebat's people having called them but me not
  which(!(d1$FID %in% cnv_cov$FID))[1]

 #Anyone who is missing gets median imputed..
  cnv_cov[is.na(cnv_cov$LRR_mean),]$LRR_mean <- median(cnv_cov$LRR_mean,na.rm=T)
  cnv_cov[is.na(cnv_cov$LRR_median),]$LRR_mean <- median(cnv_cov$LRR_median,na.rm=T)
  cnv_cov[is.na(cnv_cov$LRR_SD),]$LRR_SD <- median(cnv_cov$LRR_SD,na.rm=T)
  cnv_cov[is.na(cnv_cov$BAF_mean),]$BAF_mean <- median(cnv_cov$BAF_mean,na.rm=T)
  cnv_cov[is.na(cnv_cov$BAF_median),]$BAF_median <- median(cnv_cov$BAF_median,na.rm=T)
  cnv_cov[is.na(cnv_cov$BAF_SD),]$BAF_SD <- median(cnv_cov$BAF_SD,na.rm=T)
  cnv_cov[is.na(cnv_cov$BAF_drift),]$BAF_drift <- median(cnv_cov$BAF_drift,na.rm=T)
  cnv_cov[is.na(cnv_cov$WF),]$WF <- median(cnv_cov$WF,na.rm=T)
  cnv_cov[is.na(cnv_cov$NumCNV),]$NumCNV<- median(cnv_cov$NumCNV,na.rm=T)
  cnv_cov[is.na(cnv_cov$LRR_SD_standardized),]$LRR_SD_standardized <- median(cnv_cov$LRR_SD_standardized,na.rm=T)
  
 #code sex, ptsd status usefully 
  cnv_cov$sex <- gsub("F","2", cnv_cov$gender)
  cnv_cov$sex <- gsub("M","1", cnv_cov$sex)
  cnv_cov$sex <- gsub("U","0", cnv_cov$sex)
  cnv_cov$phenoPTS <- gsub("PTSD",2,cnv_cov$dx)
  cnv_cov$phenoPTS <- gsub("CONT",1,cnv_cov$phenoPTS)
  
  cnv_cov_exp <- subset(cnv_cov,select=c(FID,IID,C1,C2,C3,C4,C5,
                                         LRR_mean,LRR_median,LRR_SD,
                                         BAF_mean,BAF_median,BAF_SD,
                                         BAF_drift,WF,LRR_SD_standardized,
                                         sex,phenoPTS,pheno,studyname2))
                                         
  write.table(cnv_cov_exp,file="ptsd_10kb_10probes.cov",quote=F,row.names=F)
                                      
 #If you feel like checking, this is a great time to see which CNV feature is correlated to PTSD
 
########################################################
####### convert the file to be run in PLINK (BASH) 
########################################################

#Put all individual studies  into one big file for PLINK
#Multiple FID line repeats because each individual file had a header. This just takes the first instance of the ehader
#cat calls_merged/*.preplink | awk '{if (NR == 1 || $1 != "FID") print}' > calls_merged/GWAS1_2_combined.cnv


#Make a .map file for PLINK
plink-1.07-x86_64/plink --noweb --cnv-list ptsd_10kb_10probes.cnv  --cnv-make-map --out ptsd_10kb_10probes

#Check for overlap (this shouldnt be a problem)
 plink-1.07-x86_64/plink --noweb --cfile ptsd_10kb_10probes  --cnv-check-no-overlap  --out all_overlap

#Get all allele frequencies
 plink-1.07-x86_64/plink --noweb --cfile ptsd_10kb_10probes --cnv-freq-method2  .5   --cnv-write --cnv-write-freq --out cnv_filtered/ptsd_10kb_10probes_freqs.cnv
  
 #Remove only previously implicated CNV
  plink-1.07-x86_64/plink --noweb --cfile ptsd_10kb_10probes         --cnv-exclude previously_implicated_kendall_both.txt  --cnv-region-overlap 0.01  --cnv-write  --out cnv_filtered/ptsd_10kb_10probes_NOimplicated
  plink-1.07-x86_64/plink --noweb --cnv-list cnv_filtered/ptsd_10kb_10probes_NOimplicated.cnv --fam cnv_filtered/ptsd_10kb_10probes_NOimplicated.fam --cnv-make-map  --out cnv_filtered/ptsd_10kb_10probes_NOimplicated

  plink-1.07-x86_64/plink --noweb --cfile ptsd_10kb_10probes         --cnv-exclude previously_implicated_kendall_dups.txt  --cnv-dup  --cnv-region-overlap 0.01  --cnv-write  --out cnv_filtered/ptsd_10kb_10probes_NOimplicateddup
  plink-1.07-x86_64/plink --noweb --cnv-list cnv_filtered/ptsd_10kb_10probes_NOimplicateddup.cnv --fam cnv_filtered/ptsd_10kb_10probes_NOimplicateddup.fam --cnv-make-map  --out cnv_filtered/ptsd_10kb_10probes_NOimplicateddup
 
  plink-1.07-x86_64/plink --noweb --cfile ptsd_10kb_10probes         --cnv-exclude previously_implicated_kendall_dels.txt  --cnv-del   --cnv-region-overlap 0.01  --cnv-write  --out cnv_filtered/ptsd_10kb_10probes_NOimplicateddel
  plink-1.07-x86_64/plink --noweb --cnv-list cnv_filtered/ptsd_10kb_10probes_NOimplicateddel.cnv --fam cnv_filtered/ptsd_10kb_10probes_NOimplicateddel.fam --cnv-make-map  --out cnv_filtered/ptsd_10kb_10probes_NOimplicateddel

 #Include only previously implicated CNV
  plink-1.07-x86_64/plink --noweb --cfile ptsd_10kb_10probes         --cnv-intersect previously_implicated_kendall_both.txt --cnv-write  --out cnv_filtered/ptsd_10kb_10probes_implicated
  plink-1.07-x86_64/plink --noweb --cnv-list cnv_filtered/ptsd_10kb_10probes_implicated.cnv --fam cnv_filtered/ptsd_10kb_10probes_implicated.fam --cnv-make-map  --out cnv_filtered/ptsd_10kb_10probes_implicated

  plink-1.07-x86_64/plink --noweb --cfile ptsd_10kb_10probes         --cnv-intersect previously_implicated_kendall_dups.txt --cnv-write --cnv-dup --out cnv_filtered/ptsd_10kb_10probes_implicateddup
  plink-1.07-x86_64/plink --noweb --cnv-list cnv_filtered/ptsd_10kb_10probes_implicateddup.cnv --fam cnv_filtered/ptsd_10kb_10probes_implicateddup.fam --cnv-make-map  --out cnv_filtered/ptsd_10kb_10probes_implicateddup
  
  plink-1.07-x86_64/plink --noweb --cfile ptsd_10kb_10probes         --cnv-intersect previously_implicated_kendall_dels.txt --cnv-write --cnv-del --out cnv_filtered/ptsd_10kb_10probes_implicateddel
  plink-1.07-x86_64/plink --noweb --cnv-list cnv_filtered/ptsd_10kb_10probes_implicateddel.cnv --fam cnv_filtered/ptsd_10kb_10probes_implicateddel.fam --cnv-make-map  --out cnv_filtered/ptsd_10kb_10probes_implicateddel


 #Include only previously implicated CNV, 50% overlap
  plink-1.07-x86_64/plink --noweb --cfile ptsd_10kb_10probes      --cnv-region-overlap 0.50   --cnv-intersect previously_implicated_kendall_both.txt --cnv-write  --out cnv_filtered/ptsd_10kb_10probes_50implicated
  plink-1.07-x86_64/plink --noweb --cnv-list cnv_filtered/ptsd_10kb_10probes_50implicated.cnv --fam cnv_filtered/ptsd_10kb_10probes_50implicated.fam --cnv-make-map  --out cnv_filtered/ptsd_10kb_10probes_50implicated

  plink-1.07-x86_64/plink --noweb --cfile ptsd_10kb_10probes      --cnv-region-overlap 0.50   --cnv-intersect previously_implicated_kendall_dups.txt --cnv-write --cnv-dup --out cnv_filtered/ptsd_10kb_10probes_50implicateddup
  plink-1.07-x86_64/plink --noweb --cnv-list cnv_filtered/ptsd_10kb_10probes_50implicateddup.cnv --fam cnv_filtered/ptsd_10kb_10probes_50implicateddup.fam --cnv-make-map  --out cnv_filtered/ptsd_10kb_10probes_50implicateddup
  
  plink-1.07-x86_64/plink --noweb --cfile ptsd_10kb_10probes       --cnv-region-overlap 0.50  --cnv-intersect previously_implicated_kendall_dels.txt --cnv-write --cnv-del --out cnv_filtered/ptsd_10kb_10probes_50implicateddel
  plink-1.07-x86_64/plink --noweb --cnv-list cnv_filtered/ptsd_10kb_10probes_50implicateddel.cnv --fam cnv_filtered/ptsd_10kb_10probes_50implicateddel.fam --cnv-make-map  --out cnv_filtered/ptsd_10kb_10probes_50implicateddel




 #Include only CNV overlapping genes
  tail -n+2 glist-hg19_genesymbol.txt > glist-hg19_genesymbol.txt2
  plink-1.07-x86_64/plink --noweb --cfile ptsd_10kb_10probes         --cnv-intersect glist-hg19_genesymbol.txt2 --cnv-write  --out cnv_filtered/ptsd_10kb_10probes_geneoverlap
  plink-1.07-x86_64/plink --noweb --cnv-list cnv_filtered/ptsd_10kb_10probes_geneoverlap.cnv --fam cnv_filtered/ptsd_10kb_10probes_geneoverlap.fam --cnv-make-map  --out cnv_filtered/ptsd_10kb_10probes_geneoverlap
  plink-1.07-x86_64/plink --noweb --cfile cnv_filtered/ptsd_10kb_10probes_geneoverlap --cnv-freq-method2  .5  --cnv-exclude previously_implicated_kendall_both.txt  --cnv-write --cnv-write-freq --out cnv_filtered/ptsd_10kb_10probes_geneoverlap_freqs.cnv
  
  plink-1.07-x86_64/plink --noweb --cfile ptsd_10kb_10probes         --cnv-intersect glist-hg19_genesymbol.txt2 --cnv-write --cnv-dup --out cnv_filtered/ptsd_10kb_10probes_geneoverlapdup
  plink-1.07-x86_64/plink --noweb --cnv-list cnv_filtered/ptsd_10kb_10probes_geneoverlapdup.cnv --fam cnv_filtered/ptsd_10kb_10probes_geneoverlapdup.fam --cnv-make-map  --out cnv_filtered/ptsd_10kb_10probes_geneoverlapdup
   plink-1.07-x86_64/plink --noweb --cfile cnv_filtered/ptsd_10kb_10probes_geneoverlapdup --cnv-exclude previously_implicated_kendall_both.txt --cnv-freq-method2  .5   --cnv-write --cnv-write-freq --out cnv_filtered/ptsd_10kb_10probes_geneoverlapdup_freqs.cnv
 
  plink-1.07-x86_64/plink --noweb --cfile ptsd_10kb_10probes         --cnv-intersect glist-hg19_genesymbol.txt2 --cnv-write --cnv-del --out cnv_filtered/ptsd_10kb_10probes_geneoverlapdel
  plink-1.07-x86_64/plink --noweb --cnv-list cnv_filtered/ptsd_10kb_10probes_geneoverlapdel.cnv --fam cnv_filtered/ptsd_10kb_10probes_geneoverlapdel.fam --cnv-make-map  --out cnv_filtered/ptsd_10kb_10probes_geneoverlapdel
  plink-1.07-x86_64/plink --noweb --cfile cnv_filtered/ptsd_10kb_10probes_geneoverlapdel --cnv-freq-method2  .5 --cnv-exclude previously_implicated_kendall_both.txt  --cnv-write --cnv-write-freq --out cnv_filtered/ptsd_10kb_10probes_geneoverlapdel_freqs.cnv



 #Individual effects of each individual CNV    
 IFS=$'\n'
 for psychiatric_cnv in  $(cat previously_implicated_kendall_numbered.txt )
 do
   flag=$(echo $psychiatric_cnv | awk '{print $5}') #for dup or del
   cnvname=$(echo $psychiatric_cnv | awk '{print $4}') #name here, not used 
   cnvnum=$(echo $psychiatric_cnv | awk '{print $6}') #order number here
   if [ $flag == "DUP" ]
   then
    plinkflag="--cnv-dup"
   fi
   if [ $flag == "DEL" ]
   then
    plinkflag="--cnv-del"
   fi 
   echo $psychiatric_cnv
   
  echo $psychiatric_cnv > temp/$cnvnum
  
  #any overlap
  # plink-1.07-x86_64/plink --noweb --cfile ptsd_10kb_10probes         --cnv-intersect temp/$cnvnum --cnv-write  $plinkflag --out cnv_filtered_neuropsych/ptsd_10kb_10probes_$cnvnum
  # plink-1.07-x86_64/plink --noweb --cnv-list cnv_filtered_neuropsych/ptsd_10kb_10probes_"$cnvnum".cnv --fam cnv_filtered_neuropsych/ptsd_10kb_10probes_"$cnvnum".fam --cnv-make-map  --out cnv_filtered_neuropsych/ptsd_10kb_10probes_"$cnvnum"

  #40% overlap only
   #plink-1.07-x86_64/plink --noweb --cfile ptsd_10kb_10probes    --cnv-region-overlap 0.40     --cnv-intersect temp/$cnvnum --cnv-write  $plinkflag --out cnv_filtered_neuropsych/ptsd_overlap_10kb_10probes_"$cnvnum"
  # plink-1.07-x86_64/plink --noweb --cnv-list cnv_filtered_neuropsych/ptsd_overlap_10kb_10probes_"$cnvnum".cnv --fam cnv_filtered_neuropsych/ptsd_overlap_10kb_10probes_"$cnvnum".fam --cnv-make-map  --out cnv_filtered_neuropsych/ptsd_overlap_10kb_10probes_"$cnvnum"

  #50% overlap only
   plink-1.07-x86_64/plink --noweb --cfile ptsd_10kb_10probes    --cnv-region-overlap 0.50     --cnv-intersect temp/$cnvnum --cnv-write  $plinkflag --out cnv_filtered_neuropsych/ptsd_overlap50_10kb_10probes_"$cnvnum"
   plink-1.07-x86_64/plink --noweb --cnv-list cnv_filtered_neuropsych/ptsd_overlap50_10kb_10probes_"$cnvnum".cnv --fam cnv_filtered_neuropsych/ptsd_overlap50_10kb_10probes_"$cnvnum".fam --cnv-make-map  --out cnv_filtered_neuropsych/ptsd_overlap50_10kb_10probes_"$cnvnum"

  #100% overlap only
   plink-1.07-x86_64/plink --noweb --cfile ptsd_10kb_10probes    --cnv-region-overlap 1     --cnv-intersect temp/$cnvnum --cnv-write  $plinkflag --out cnv_filtered_neuropsych/ptsd_overlap100_10kb_10probes_"$cnvnum"
   plink-1.07-x86_64/plink --noweb --cnv-list cnv_filtered_neuropsych/ptsd_overlap100_10kb_10probes_"$cnvnum".cnv --fam cnv_filtered_neuropsych/ptsd_overlap100_10kb_10probes_"$cnvnum".fam --cnv-make-map  --out cnv_filtered_neuropsych/ptsd_overlap100_10kb_10probes_"$cnvnum"


  #and a flagless version for either
  
 #  plink-1.07-x86_64/plink --noweb --cfile ptsd_10kb_10probes         --cnv-intersect temp/$cnvnum --cnv-write  --out cnv_filtered_neuropsych/ptsd_10kb_10probes_both_$cnvnum
 #  plink-1.07-x86_64/plink --noweb --cnv-list cnv_filtered_neuropsych/ptsd_10kb_10probes_both_"$cnvnum".cnv --fam cnv_filtered_neuropsych/ptsd_10kb_10probes_both_"$cnvnum".fam --cnv-make-map  --out cnv_filtered_neuropsych/ptsd_10kb_10probes_both_"$cnvnum"

 done



#Leftover code for overlaps:
   # R
   # library(data.table)
   # d1 <- fread('all.cnv',data.table=F)
   # overlap <- fread('all_overlap.cnv.overlap',data.table=F)
   # overlap$remove <- 0 
   # overlap[duplicated(overlap[,c("FID","IID","CHR")]),]$remove <- 1

   # d1m <- merge(d1,overlap,by=c("FID","IID","CHR","BP1","BP2"),all=TRUE)

   # d1m2 <- subset(d1m,remove!=1 | is.na(remove),select=c(FID,IID,CHR,BP1,BP2,TYPE,SCORE,SITES))

   # write.table(d1m2, file="all2.cnv",quote=F,row.names=F,sep="\t")

   # cp all.fam all2.fam

#Leftover code for manual maf filtering:

#114383 individuals, MAF 0.01 must be 1,144
 R
 library(data.table)
 cnvf <- fread('cnv_filtered/ptsd_10kb_10probes_freqs.cnv.cnv')
 #The 100%ile is 1019, such that these must already be filtered!!
 
 #ergo, use ptsd_10kb_10probes
 
 
 #Filter to < 1% 
  # awk 'BEGIN{OFS="\t"}{if (NR==1 || $9 <= 1144) print $1,$2,$3,$4,$5,$6,$7,$8}' cnv_filtered/ptsd_10kb_10probes_freqs.cnv > cnv_filtered/ptsd_10kb_10probes_maf01
  #Make a MAP file
  # plink-1.07-x86_64/plink  --noweb --cnv-list  cnv_filtered/ptsd_10kb_10probes_freqs.cnv.cnv --fam cnv_filtered/ptsd_10kb_10probes_freqs.cnv.fam --cnv-make-map  --out cnv_filtered/ptsd_10kb_10probes_maf01
  #Make a CNV file 
 #  plink-1.07-x86_64/plink  --noweb --cnv-list  cnv_filtered/ptsd_10kb_10probes_freqs.cnv.cnv --fam cnv_filtered/ptsd_10kb_10probes_freqs.cnv.fam --map cnv_filtered/ptsd_10kb_10probes_maf01.cnv.map --cnv-write --out cnv_filtered/ptsd_10kb_10probes_maf01

#Get those overlapping psychiatric




#Leftover code to add my extra studies:

# ### add My extra studies
# #Important to make sure that everyone in my data 1) was attempted for CNV calls 2) calling did not fail, its just they just had 0 CNVs

 # table(d1$found,d1$studyname) #peole from each study have been found in this file, so it appears to be fine
 
 # table(cnvtypes$sampleId %in% d1a$identifier)
 # notfound <- which(!cnvtypes$sampleId %in% d1a$identifier)
 # cnvtypes[notfound[1:5],]
 
 
# names(dat) <- c("CHR","BP1","BP2","SUBJ","LENGTH","STRAND","JUNK1","JUNK2","CNVTYPE")
# dat$CHR <- gsub("chr","",dat$CHR)
# dat$SITES <- dat$LENGTH
# dat$TYPE <- sapply(dat$CNVTYPE,cnv_conv)
# dat$SCORE <- 0 #Dummt code for confidence score. Not in use right now.

# #Need to convert the IDs to a PLINK format using a conversion key

# #read in sample IDs conversion sheet (format is: subject genotyping id, subject FID, subject IID)
# sample_ids <- read.table(paste(study,'_conversion_key.txt',sep=''), header=T,na.strings=c(NA,"#N/A"))
# names(sample_ids) <- c("SUBJ", "FID", "IID")

# datA <- merge(dat, sample_ids,by="SUBJ")

# dat_exp <- subset(datA, select=c(FID,IID,CHR,BP1,BP2,TYPE,SCORE,SITES))	

# write.table(dat_exp, paste("calls_merged/", study,"_sample_qced.cnv",sep=""),quote=F,row.names=F)

# q()
# n