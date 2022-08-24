#Write out frequencies
#plink-1.07-x86_64/plink --noweb --cfile ptsd_10kb_10probes --cnv-freq-method2  .5   --cnv-write --cnv-write-freq --out ptsd_10kb_10probes_freqs
 
#CNV count function only needs to be run on subjects who are in the cnv file, 
#otherwise they will be assigned value 0
 num_cnv <- function(x,dataframe)
 {
  subject_cnvs <- which( dataframe$FID_IID == x )
  if(length(subject_cnvs) >= 1)
  {
   #Find the CNV entries for a given subject
   subject_cnv_data <- dataframe[subject_cnvs,]   
   #Calculate number of CNVs
    cnv_count <- dim(subject_cnv_data)[1]
   #distance spanned by CNV (in KB)
    cnv_span <- sum(subject_cnv_data$span) / 1000000
   #Average span of CNV
    cnv_avg <- cnv_span / cnv_count
   #Return these values
   return(c(cnv_count,cnv_span,cnv_avg))
  } else { return(c(0,0,0))} #If there are no cnvs, return zeros
 }
 
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
 harmonize_scales <-fread('effect_scale_conversion.txt',data.table=F)[,-1]
  names(harmonize_scales)[1] <- "study"
  harmonize_scales$prevalence <- harmonize_scales$N_Cases/harmonize_scales$N
  harmonize_scales$minprevalence <- harmonize_scales$N_Controls/harmonize_scales$N
 
 
 
#Load global CNV parameters
 cnvparms <- fread('cnv_global_parameters_dec17_2021.csv',data.table=F)
 
#Results will be stored in this style of matrix
 resmat0 <- as.data.frame(matrix(nrow=length(studylist),ncol=4))
 names(resmat0) <-c ("beta","se","t","p")
 resmat0$study <- studylist
 resmat <- merge(resmat0,harmonize_scales,by="study")
 
cnv_data <- fread('ptsd_10kb_10probes_freqs.cnv',data.table=F) 
  #FID_IID to link data 
   cnv_data$FID_IID <- paste(cnv_data$FID,cnv_data$IID,sep="_")
  #calculate CNV span now because it has to be done a lot
   cnv_data$span <- cnv_data$BP2 - cnv_data$BP1 
  
  
cnv_data_subset <- subset(cnv_data,TYPE==3)


   #Make a matrix to store cnv count per subject
   
    cnv_fam_subset <- cnv_fam
   #The CNV count detection only gets applied to subjects who appear in the CNV file to save time. 
   #Everyone else is imputed with zero
    cnv_fam_subset$cnv_count <- 0
    cnv_fam_subset$cnv_span <- 0
    cnv_fam_subset$cnv_avg <- 0
   #which subjects have data?
   examine_subjects_list <- unique(cnv_data_subset$FID_IID)
    examine_subjects_rows <- which(cnv_fam$FID_IID %in% examine_subjects_list)
   #Get the parameters for these subjects
    cnv_fam_subset[examine_subjects_rows, c("cnv_count","cnv_span","cnv_avg")] <- t(sapply(cnv_fam_subset[examine_subjects_rows,]$FID_IID,num_cnv,dataframe=cnv_data_subset))
   #write a file to store this data,
    #write.table(subset(cnv_fam_subset,select=c(FID,IID,cnv_count,cnv_span,cnv_avg)),file=paste("global_covar/",cnvparms[i,]$Analysis,"_",cnvparms[i,]$Type,sep=""),row.names=F)
   
   
   #total number of carriers
    by(cnv_fam_subset$cnv_count >0,cnv_fam_subset$studyname2,sum)
    
    #average counts per person 
   mean(cnv_fam_subset$cnv_count )
   sd(cnv_fam_subset$cnv_count )
  
   #total CNV count
   sum(cnv_fam_subset$cnv_count)
   ddply(cnv_fam_subset, ~studyname2, colwise(sum,.cols="cnv_count"))
     

   #Split to just carriers
   cnv_fam_subsetx <- subset(cnv_fam_subset,cnv_count >0)

  #Avg count of CNV among carriers
   as.vector( by(cnv_fam_subsetx$cnv_count ,cnv_fam_subsetx$studyname2,mean))
   as.vector( by(cnv_fam_subsetx$cnv_count,cnv_fam_subsetx$studyname2,sd))
    
   #Avg length of CNV among carriers
    as.vector(by(cnv_fam_subsetx$cnv_span ,cnv_fam_subsetx$studyname2,mean))
     as.vector(by(cnv_fam_subsetx$cnv_span ,cnv_fam_subsetx$studyname2,sd))
      
   #Avg of avg length of CNV among carriers
     as.vector(by(cnv_fam_subsetx$cnv_avg ,cnv_fam_subsetx$studyname2,mean))
    as.vector( by(cnv_fam_subsetx$cnv_avg ,cnv_fam_subsetx$studyname2,sd))
    
    

   #across all studies
   library(plyr)
   

   #N subjects with CNVs
   table(cnv_fam_subset$cnv_count >0)
    
   #among carrirers, mean values
    mean(cnv_fam_subsetx$cnv_span)
    sd(cnv_fam_subsetx$cnv_span)
    mean(cnv_fam_subsetx$cnv_avg)
    sd(cnv_fam_subsetx$cnv_avg)  
    
    
    
    cnv_fam_subset2 <- subset(cnv_fam_subset,studyname2 == study)
    
    
    # cnv_fam_subset_no <- subset(cnv_fam_subset,studyname2 != study & studyname2 != "ukbb")
    
    # #leave out ukbb in this and see what happens)
    # print(study)
    # print(wilcox.test(cnv_fam_subset2$cnv_count,cnv_fam_subset_no$cnv_count))
    # print(ks.test(cnv_fam_subset2$cnv_count,cnv_fam_subset_no$cnv_count))
    # }
  

    library(psych) #apply dsecribe function by study, send results to a single demographic vector
    
    tf <- function(x)
    tf
   ddply(cnv_fam_subset,~studyname2, cnv_count >0)
   
   sort(by(cnv_fam_subsetx$cnv_span,cnv_fam_subsetx$studyname2,mean) )
    by(cnv_fam_subsetx$cnv_count,cnv_fam_subsetx$studyname2,mean) 
 
   
   

    