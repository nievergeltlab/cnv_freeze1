library(data.table) 
library(metafor)
library(ordinal)
 num_cnv <- function(x,dataframe)
 {
  subject_cnvs <- which( dataframe$FID_IID == x )
  if(length(subject_cnvs) >= 1)
  {
   #Find the CNV entries for a given subject
   subject_cnv_data <- dataframe[subject_cnvs,]   
   #Calculate number of CNVs
    cnv_count <- dim(subject_cnv_data)[1]
   #distance spanned by CNV (in MB)
    cnv_span <- sum(subject_cnv_data$span) / 1000000
   #Average span of CNV
    cnv_avg <- cnv_span / cnv_count
   #Return these values
   return(c(cnv_count,cnv_span,cnv_avg))
  } else { return(c(0,0,0))} #If there are no cnvs, return zeros
 }
 
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
   
#In this case, CNV parameters are the file from kendall
 cnvparms <- fread('previously_implicated_kendall_numbered.txt',data.table=F)

#Results will be stored in this style of matrix
 resmat0 <- as.data.frame(matrix(nrow=length(studylist),ncol=7))
 names(resmat0) <-c ("beta","se","t","p","count","freq","Nsub")

 resmat0$study <- studylist
 resmat <- merge(resmat0,harmonize_scales,by="study")
 
 
 #do in both and in the specific CNV
 for (feature in "one")#"both","one"))
 {
 #for main analysis; 26 fails, 47 fails, 50 fails - 1:nrow(cnvparms)[-c(26,47,50)]
 #For one: 
 #main analysis, 0% overlap: c(1:9,11:25,27:46,48:49,51:nrow(cnvparms)) 10 failes,26
 
  #  1:11,
 for (i in 1:nrow(cnvparms)) # ) c(c(1:9,11:25,27:46,48:49,51:nrow(cnvparms))) #
 {
  print(i)
  #store results
   resmat2 <- resmat
  
   #Load CNV data - keep in mind the raw data has insertions and deletions in it no matter what, so must be subste
   #modifying this can make it both or not
   #if (feature == "both")
   #{
   # cnv_data <- fread(paste('cnv_filtered_neuropsych/ptsd_10kb_10probes_both_',i,'.cnv',sep=''),data.table=F) 
   #}
   if(feature =="one")
   {
   # NOTE: SET TO 100% overlap right now! ptsd_overlap_ change to 'ptsd' if 0% overlap needed!
   #ptsd_10kb_10probes #touching
   #ptsd_overlap50_10kb_10probes #50pct
   #ptsd_overlap100_10kb_10probes #100pct
    cnv_data <- fread(paste('cnv_filtered_neuropsych/ptsd_10kb_10probes_',i,'.cnv',sep=''),data.table=F) 
    #cnv_data <- fread(paste('cnv_filtered_neuropsych/ptsd_overlap50_10kb_10probes_',i,'.cnv',sep=''),data.table=F) 
    #cnv_data <- fread(paste('cnv_filtered_neuropsych/ptsd_overlap100_10kb_10probes_',i,'.cnv',sep=''),data.table=F) 
     
   }
   
    cnv_data$FID_IID <- paste(cnv_data$FID,cnv_data$IID,sep="_")
   #calculate CNV span now because it has to be done a lot
    cnv_data$span <- cnv_data$BP2 - cnv_data$BP1 
  
  #Subset to relevant CNVs - in this analysis, all are relevant.
   cnv_data_subset <- subset(cnv_data) 
   
   #FOR THIS BURDEN analysis, we can forgo span and frequency requirements!
   #we can also forgo type subsetting, already done
                        
 
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
  
  #loop analysis over CNV type. 
  for (cnvtype in c("cnv_count"))
  {
   #And over study
   for (study in studylist )
   {
    cnv_fam_subset2 <- subset(cnv_fam_subset,studyname2 == study)
    cnv_fam_subset2$cnv <- cnv_fam_subset2[,cnvtype]

    
   # summary(glm.nb(pheno -6 ~ C1 + C2 + C3 + C4 + C5 + LRR_SD + cnv, data=cnv_fam_subset2))
    #m1 <- lm(pheno ~ C1 + C2 + C3 + C4 + C5 + LRR_SD + cnv, data=cnv_fam_subset2)
    #coeftest(m1,vcov. = vcovHC(m1,type="HC3"))
    #summary(glm(cnv_fam_subset2$pheno>=23 ~ cnv_fam_subset2$cnv,family="binomial"))
    
    if (!(study %in% c("psy3_comb","pts1_pts1")))
    {
     results <- try(summary(clm(as.factor(pheno) ~ C1 + C2 + C3 + C4 + C5 + LRR_SD + cnv, data=cnv_fam_subset2))$coefficients ,silent=TRUE)
     if(class(results) == "try-error")
     {
      results <- matrix(c(NA,NA,NA,NA),nrow=1,ncol=4)
      #row.names(results) <- "cnv"
     }
    }
       
    if (study %in% c("psy3_comb","pts1_pts1"))
    {
     results <- try(summary(glm(pheno-1 ~ C1 + C2 + C3 + C4 + C5 + LRR_SD + cnv, data=cnv_fam_subset2,family="binomial"))$coefficients)
    }
    #save regression results
    if(row.names(results)[nrow(results)] != "cnv")
    { 
     resmat2[which(resmat2$study == study),c("beta","se","t","p")] <- rep(NA,4)
    } else {
    resmat2[which(resmat2$study == study),c("beta","se","t","p")] <- results[nrow(results),]
     }
     resmat2[which(resmat2$study == study),]$count <- sum(cnv_fam_subset2$cnv_count >= 1)
     resmat2[which(resmat2$study == study),]$freq <- sum(cnv_fam_subset2$cnv_count  >= 1) /nrow(cnv_fam_subset2)
     resmat2[which(resmat2$study == study),]$Nsub <- nrow(cnv_fam_subset2)
    
   }
   
   #Rescale everything to odds scale
   # ccstudies <- which(resmat2$study %in% c("psy3_comb","pts1_pts1")) #these will be transformed onto the lienar regression scale again..
   # resmat2$beta_scale <- resmat2$beta *resmat2$conversion_factor
   # resmat2$se_scale <- resmat2$se *resmat2$conversion_factor
   # resmat2[ccstudies,]$beta_scale <- resmat2[ccstudies,]$beta_scale   *resmat2[ccstudies,]$prevalence*resmat2[ccstudies,]$minprevalence
   # resmat2[ccstudies,]$se_scale <- resmat2[ccstudies,]$se_scale *resmat2[ccstudies,]$prevalence*resmat2[ccstudies,]$minprevalence
     
   
   
   #Only do analyses if there are values to test...
   if(! all(is.na(resmat2$p)))
   {
   outfile=paste("neurocnv_0overlap_",i,"_",feature,"_",cnvtype,sep="")
   write.table(resmat2,file=paste("outputs_psychiatric_cnvlevel_jun2022/","bystudy_",outfile,'.txt',sep=''),row.names=F)
   meta_res <- try(rma(yi=resmat2$beta,sei=resmat2$se,slab=resmat2$study,method="FE"))
   capture.output(meta_res,paste("outputs_psychiatric_cnvlevel_jun2022/","meta_",outfile,'.txt',sep=''))
   meta_out <- c(t(meta_res[c("beta","se","zval","pval")]),sum(resmat2$count),sum(resmat2$count)/sum(resmat2$Nsub),sum(resmat2$count>=1)) #its going to be a bit weird since some studies cut
   #count>=1 establishes N studies included
   names(meta_out) <- c("beta","se","zval","pval","count","freq", "nstudy")
   write.table(meta_out,paste("outputs_psychiatric_cnvlevel_jun2022/","meta2_",outfile,'.txt',sep=''),row.names=F)

  pdf(paste("outputs_psychiatric_cnvlevel_jun2022/","meta_",outfile,'.pdf',sep=''),7,7)
   plot(meta_res)
  dev.off()
  pdf(paste("outputs_psychiatric_cnvlevel_jun2022/","forest_",outfile,'.pdf',sep=''),7,7)
   forest(meta_res)
  dev.off()
  }


  }
 }  
 }
   #Export all results to a single file
 for files in $(ls outputs_psychiatric_cnvlevel_jun2022 | grep meta2 | grep 100overlap )
 do
  output=$( tail -n1 outputs_psychiatric_cnvlevel_jun2022/$files)
  echo $files $output >> cnv_burden_psychiatric_cnvlevel_100overlap_or_jun15_2022.txt
 done
 
 
 
  for files in $(ls outputs_psychiatric_cnvlevel_jun2022 | grep meta2 | grep 50overlap )
 do
  output=$( tail -n1 outputs_psychiatric_cnvlevel_jun2022/$files)
  echo $files $output >> cnv_burden_psychiatric_cnvlevel_50overlap_or_jun15_2022.txt
 done
 
  for files in $(ls outputs_psychiatric_cnvlevel_jun2022 | grep meta2 | grep _0overlap )
 do
  output=$( tail -n1 outputs_psychiatric_cnvlevel_jun2022/$files)
  echo $files $output >> cnv_burden_psychiatric_cnvlevel_0overlap_or_jun15_2022.txt
 done
 