
##CNV data handling
 library(Matrix)
 library(data.table)
  
#Read covariates (includes complete sample list)
 cnv_fam <- fread('ptsd_10kb_10probes.cov',data.table=F)
  cnv_fam$FID_IID <- paste(cnv_fam$FID,cnv_fam$IID,sep="_")
 
 #There is a mistake - psy3_comb should be case control. I use the dx phenotype instead
  cnv_fam[which(cnv_fam$studyname2 == "psy3_comb"),]$pheno <- cnv_fam[which(cnv_fam$studyname2 == "psy3_comb"),]$phenoPTS
 #PPDS is rounded to nearest whole number (due to harmonizing, I created some decimals..)
  cnv_fam[which(cnv_fam$studyname2 == "nss1_ppds"),]$pheno <- round(cnv_fam[which(cnv_fam$studyname2 == "nss1_ppds"),]$pheno)
 #GUTS and NHS2 will be analyzed together
  cnv_fam[which(cnv_fam$studyname2 == "psy2_guts"),]$studyname2 <- "psy2_nhsy" 
 
#List all studies (this may be a problem later if a study has no subjects to analyze
 studylist <- unique(cnv_fam$studyname2)
 
#Load effect size conversion sheet to harmonize PTSD
 #rescale each PTSD measure before exporting? logically I can even do the prevalence rescale first! 
  harmonize_scales <-fread('effect_scale_conversion.txt',data.table=F)[,-1]
  names(harmonize_scales)[1] <- "study"
  harmonize_scales$prevalence <- harmonize_scales$N_Cases/harmonize_scales$N
  harmonize_scales$minprevalence <- harmonize_scales$N_Controls/harmonize_scales$N
  
#Read in the CNV call data, with full rank for gene names (see 00_genelist_annotat.txt)
 genesA <- fread("ptsd_10kb_10probes_freqs_annotated_fullrank.cnv",data.table=F)
 names(genesA)[dim(genesA)[2]] <- "geneSymbol"

#'NA' is not actually a gene name. These must be excluded
 genesA <- subset(genesA, !is.na(geneSymbol) || geneSymbol == "NA")
 
#Create a unique identifier for each subject 
 genesA$FID_IID <- paste(genesA$FID,genesA$IID)

#increment by 1000 genes at a time
increment=1000
for (gtype in c("all","dup","del")) # all",
{
 #Subset to dups and dels if needed
  if(!(gtype %in% c("dup","del")))
  {
   genes <- genesA
  }
  if(gtype == "dup")
  {
   genes <- subset(genesA,TYPE==3)
  }
  if(gtype == "del")
  {
   genes <- subset(genesA,TYPE==1)
  }
  
 #List all gene names. Gene matrix columns will be ordered the same as this (i.e. this is the ordered set of column names)
  genelist <- unique(genes$geneSymbol)
 #List of all subject names that are in your CNV analysis.
  #Gene matrix rows will be ordered the same as this (i.e. this is the ordered set of rows)
  subjects <- cnv_fam
 #Create a unique identifier for each subject
  subjects$FID_IID <- paste(subjects$FID,subjects$IID)
 #Assign coordinate entries to all subjects of the gene matrix, to populate a sparse matrix
  genes$coord1 <- match(genes$FID_IID, subjects$FID_IID)
  genes$coord2 <- match(genes$geneSymbol, genelist)
 #Create a sparse matrix to hold the gene data
  genes_matrix <- sparseMatrix(i=genes$coord1,j=genes$coord2,dims=c(nrow(subjects),length(genelist)))
 #1 row per subject, 1 column per gene. Binary matrix for has cnv versus not.
 #
  #use the apply function to get frequencies, we want 0.0001 freq for all genes...
  #11 subjects.
  gmat_freq <- apply(genes_matrix,2,sum)
  maf_filtered <- which(gmat_freq>=11)
  gmat1 <- as.data.frame(genelist[maf_filtered])
  names(gmat1) <- "geneSymbol"
  write.table(gmat1,file=paste("freq_filtered_geneslist_",gtype,".txt",sep=""),row.names=F)
  
  #you could just subset to the above list and be done with it... 
  }
  
  
  #Analysis: split the genes matrix up by 1000s
  for (i in 1:ceiling((dim(genes_matrix)[2]/increment)))
  {
   #have to use a maximum otherwise prolems occur
   limit=minsize <- min(c(i*increment,dim(genes_matrix)[2]))
   
   genes_subset <- as.matrix(genes_matrix[,c(((i-1)*increment+1):limit)]) #
   genelist_subset <- genelist[c(((i-1)*increment+1):limit)] #gene list names
    
   #does global burden need to be a covariate?

   #Create a matrix to store regression coefficients
   #Loop over studies
    for (study in studylist )
    {
     results_matrix <- as.data.frame(matrix(nrow=length(genelist_subset),ncol=11))
     names(results_matrix) <- c("study","Gene","Count","Freq","nsub","beta","se","z","P","beta_scale","se_scale")
     results_matrix$Gene <- genelist_subset
     results_matrix$study <- study
     harmonize_scales_use <- harmonize_scales[which(harmonize_scales$study==study),]
     
     #Chop out the relevant pepolpe
     analyze_these <- which(cnv_fam$studyname2 == study)
     cnv_fam_subset <- cnv_fam[analyze_these,]
     genes_subset_study <- genes_subset[analyze_these,]
     
     nsub <- dim(cnv_fam_subset)[1]
     #loop over ever cnv
     for (geneloop in 1:length(genelist_subset))
     {
      #And over study
      cnv_fam_subset$cnv <- genes_subset_study[,geneloop]
      results <- summary(lm(pheno ~ C1 + C2 + C3 + C4 + C5 + LRR_SD + cnv, data=cnv_fam_subset))$coefficients
      if (study %in% c("psy3_comb","pts1_pts1"))
      {
       results <- summary(glm(pheno-1 ~ C1 + C2 + C3 + C4 + C5 + LRR_SD + cnv, data=cnv_fam_subset,family="binomial"))$coefficients
      }
     
      if(row.names(results)[nrow(results)] != "cnvTRUE")
      {
       results_matrix[geneloop,6:9] <- rep(NA,4)
      } else {   
         results_matrix[geneloop,6:9] <- results[nrow(results),] 
          }
       count <- sum(genes_subset_study[,geneloop] ==TRUE)
       results_matrix[geneloop,"Count"] <- count
       results_matrix[geneloop,"Freq"] <- count / nsub
       results_matrix[geneloop,"nsub"] <- nsub
       
     }
   
     results_matrix$beta_scale <- results_matrix$beta *harmonize_scales_use$conversion_factor
     results_matrix$se_scale <- results_matrix$se *harmonize_scales_use$conversion_factor
     if(results_matrix$study[1] %in%c("psy3_comb","pts1_pts1"))
     {
      results_matrix$beta_scale <- results_matrix$beta_scale   *harmonize_scales_use$prevalence*harmonize_scales_use$minprevalence
      results_matrix$se_scale   <- results_matrix$se_scale     *harmonize_scales_use$prevalence*harmonize_scales_use$minprevalence
     }
    
     write.table(results_matrix,file=paste("gene_basedR/",gtype,"_",study,"_",i,".txt",sep=""),row.names=F)
    }
   }
  
}

#Now in the shell ,combine this stuff together..
for gtype in all dup del
do
 for study in ftca_ftcb mrsc_comb nss1_ppds nss1_nss1 psy2_nhsy pts1_pts1 gali_gali grac_grac psy3_comb psy3_ncmh ukbb_ukbb
 do
  cat gene_basedR/"$gtype"_"$study"_*.txt | sed 's/"//g' | awk '{if (NR==1 || $2!= "Gene") print}' \
  | awk -v study=$study '{if(NR==1) {   $3=$3"_"study ; $4=$4"_"study ; $5=$5"_"study ;$6=$6"_"study ;$7=$7"_"study ;$8=$8"_"study;$9=$9"_"study;$10=$10"_"study;$11=$11"_"study} print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11 }' > gene_basedR_filtered/"$gtype"_"$study".txt
 done
done
      
#It may be easier at this point to just do METAL rather than RMA but lets stick with R...

library(plyr)
library(metafor)


for (gtype in c("all","dup","del"))
{
 #read each file into a data frame with the same name
  for (i in studylist)
  {
   assign(
    i, fread(paste('gene_basedR_filtered/',gtype,"_",i,'.txt',sep=''), data.table=F,header=T, na.strings=c('#N/A','NA','#NULL!'))
    ) 
  }

 #parse the text list of data frame names as a list of data frames
  data_list <- eval( 
     parse( 
      text=paste(
       "list(", paste(studylist, collapse=','), ")" 
       )
      )
     )

 #combine all data frames by gene. Left join, so all genes must be in each file - this should be true by default
  data_premeta <- join_all(data_list,by="Gene", type="left", match="first")

 #results will be stored here:
  meta_results <- as.data.frame(matrix(nrow=dim(data_premeta)[1], ncol=8))
  names(meta_results) <- c("Gene","count","freq","nsub","B","SE","Z","P")

#consider altering the previous analysis to include sample size already, in additon to count and freq?

 #Now for each gene... get the results
 #have to do some indexing tricks
  beta_inc=seq(9,dim(data_premeta)[2],by=9) #index columns, assuming beta starts at col 4, increments by 6
  se_inc=seq(10,dim(data_premeta)[2],by=9)
  count_inc=seq(2,dim(data_premeta)[2],by=9)
  freq_inc=seq(3,dim(data_premeta)[2],by=9)
  nsub_inc=seq(4,dim(data_premeta)[2],by=9)
  
  

   for (incgene in 1:dim(data_premeta))
   {  
    use=data_premeta[incgene,]
    
    betas=as.numeric(use[,beta_inc])
    ses=as.numeric(use[,se_inc])
    slabs=gsub("B_","",names(use[,beta_inc]))
    res<- rma(yi=betas,sei=ses,slab=slabs,method="FE")
    meta_results[incgene,]$Gene <- use$Gene
    meta_results[incgene,5:8] <- c(res$beta,res$se,res$zval,res$pval) 
    meta_results[incgene,2] <- sum(as.numeric(use[,count_inc]))
    meta_results[incgene,3] <- sum(use[,count_inc])/sum(use[,nsub_inc])
   }
   #this needs to include total count, weighted freq avg
   #nsub should only be calculated in those who have valid estimates. maybe also freq.
   write.table(meta_results,paste('gene_basedR_filtered/',gtype,"_meta",'.txt',sep=''),row.names=F)
}
     #if this takes too long, just code in the apply function..
  #Now also save a list of genes that have freq > 0.0001   
     
     
 #elastic net model also an option... not sure what for at this point but it is an option
 