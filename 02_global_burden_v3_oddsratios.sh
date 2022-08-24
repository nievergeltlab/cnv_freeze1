#Write out frequencies
#plink-1.07-x86_64/plink --noweb --cfile ptsd_10kb_10probes --cnv-freq-method2  .5   --cnv-write --cnv-write-freq --out ptsd_10kb_10probes_freqs
 library(metafor)
 library(missreg3)
 library(ordinal)
 
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
 
 #odds ratio nb glm likfun
likfun1 <- function(y,x,theta)
{ 
 #the beta parameters will be the first N-1 inputs supplied by optim, dispersion the last
 b=theta[-length(theta)]
 d=theta[length(theta)]
 loglik1 <- sum ( log(gamma(y+d)) -log(gamma(y + 1)) - log(gamma(d)) +
 y*log((1 + exp(x%*%b))^(1/d) - 1) - 
 (1 + y/d)*log(1 + exp(x%*%b)) )
 return(-loglik1)
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
 


#Calculate number of CNVs. This is applied to the CNV data
#CNV data must be subset before hand if you want e.g. specific frequency or length CNVs, or cnv type(dup or del)

 #For each CNV parameter set, do the analysis prescribed.
 #hard coded to 
 for (i in 1:nrow(cnvparms))
 {
  print(i)
  #store results
   resmat2 <- resmat
  
   #Load CNV data 
   
   #main analysis
   if(cnvparms[i,]$Include == "all")
   { cnv_data <- fread('ptsd_10kb_10probes_freqs.cnv',data.table=F) }
   
   #implicated CNVs only
   if(cnvparms[i,]$Type == "Both" & cnvparms[i,]$Include == "Implicated")
   { cnv_data <- fread('cnv_filtered/ptsd_10kb_10probes_implicated_freqs.cnv.cnv',data.table=F) }
   
    if(cnvparms[i,]$Type == "Insertion" & cnvparms[i,]$Include == "Implicated")
   { cnv_data <- fread('cnv_filtered/ptsd_10kb_10probes_implicateddup_freqs.cnv.cnv',data.table=F) }
     
   if(cnvparms[i,]$Type == "Deletion" & cnvparms[i,]$Include == "Implicated")
   { cnv_data <- fread('cnv_filtered/ptsd_10kb_10probes_implicateddel_freqs.cnv.cnv',data.table=F) }
   
    #non implicated CNVs only
   if(cnvparms[i,]$Type == "Both" & cnvparms[i,]$Include == "NoImplicated")
   { cnv_data <- fread('cnv_filtered/ptsd_10kb_10probes_NOimplicated_freqs.cnv.cnv',data.table=F) }
   
    if(cnvparms[i,]$Type == "Insertion" & cnvparms[i,]$Include == "NoImplicated")
   { cnv_data <- fread('cnv_filtered/ptsd_10kb_10probes_NOimplicateddup_freqs.cnv.cnv',data.table=F) }
     
   if(cnvparms[i,]$Type == "Deletion" & cnvparms[i,]$Include == "NoImplicated")
   { cnv_data <- fread('cnv_filtered/ptsd_10kb_10probes_NOimplicateddel_freqs.cnv.cnv',data.table=F) }
     
   
  #implicated CNVs  are removed from the gene only  analysis 
   if(cnvparms[i,]$Type == "Both" & cnvparms[i,]$Include == "genes")
   { cnv_data <- fread('cnv_filtered/ptsd_10kb_10probes_geneoverlap_freqs.cnv.cnv',data.table=F) }
   
    if(cnvparms[i,]$Type == "Insertion" & cnvparms[i,]$Include == "genes")
   { cnv_data <- fread('cnv_filtered/ptsd_10kb_10probes_geneoverlap_freqs.cnv.cnv',data.table=F) }
     
   if(cnvparms[i,]$Type == "Deletion" & cnvparms[i,]$Include == "genes")
   { cnv_data <- fread('cnv_filtered/ptsd_10kb_10probes_geneoverlap_freqs.cnv.cnv',data.table=F) }
   


   
  #FID_IID to link data 
   cnv_data$FID_IID <- paste(cnv_data$FID,cnv_data$IID,sep="_")
  #calculate CNV span now because it has to be done a lot
   cnv_data$span <- cnv_data$BP2 - cnv_data$BP1 
  
  #Subset to relevant CNVs 
   cnv_data_subset <- subset(cnv_data,
                           FREQ >= cnvparms[i,]$freqmin & FREQ <= cnvparms[i,]$freqmax &
                           span >= cnvparms[i,]$start  & span <= cnvparms[i,]$stop) #exampl: singletons
                                 
  #Subset on indels (main dataset was not already subset, this insures that it is)                  
   if(cnvparms[i,]$Type == "Insertion") 
   { cnv_data_subset <- subset(cnv_data_subset,TYPE==3)}
   if(cnvparms[i,]$Type == "Deletion") 
   { cnv_data_subset <- subset(cnv_data_subset,TYPE==1)}
   
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
  #  write.table(subset(cnv_fam_subset,select=c(FID,IID,cnv_count,cnv_span,cnv_avg)),file=paste("global_covar/",cnvparms[i,]$Analysis,"_",cnvparms[i,]$Type,sep=""),row.names=F)
   
  #This and files 03,03b should have dimension error checking
  #loop analysis over CNV type
  for (cnvtype in c("cnv_span","cnv_count","cnv_avg"))
  {
   #And over study
   for (study in studylist )
   {
   
    cnv_fam_subset2 <- subset(cnv_fam_subset,studyname2 == study)
    cnv_fam_subset2$cnv <- cnv_fam_subset2[,cnvtype]
    
    cnv_fam_subset2$sampling <- "uncond"
    cutput <- median(cnv_fam_subset2$pheno,na.rm=T)
    
    if (!(study %in% c("psy3_comb","pts1_pts1")))
    {
      #ordinal logit
       results <- try(summary(clm(as.factor(pheno) ~ C1 + C2 + C3 + C4 + C5 + LRR_SD + cnv, data=cnv_fam_subset2))$coefficients , silent = TRUE)
      # try(summary(glm.nb(pheno -5 ~ C1 + C2 + C3 + C4 + C5 + LRR_SD + cnv, data=cnv_fam_subset2))$coefficients , silent = TRUE)
       
       #Coombes
      #summary(lm(pheno ~ C1 + C2 + C3 + C4 + C5 + LRR_SD + cnv, data=cnv_fam_subset2,yCuts=0.5,lower.tail = FALSE, obstype.name=c("sampling"), xs.includes=FALSE))$coefficients
      
       #Sroka method
        # ymat <- cnv_fam_subset2$pheno - min(cnv_fam_subset2$pheno) #-min is to make sure that the distribution starts at 0, otherwise results are nonsense
        # xmat <- model.matrix (~ C1+C2+C3+C4+C5+LRR_SD+cnv,data=cnv_fam_subset2)
        
        # deltaX <- try ({
          # res1 <- optim(c(rep(0,ncol(xmat)),1), likfun1, y=ymat,x=xmat ,hessian=TRUE,method="BFGS")
        # OI<-solve(res1$hessian)
        # se <- sqrt(diag(OI))
        # tv<-res1$par/se
        # pval<-2*pnorm(tv,lower.tail=F)
        # results <-cbind(res1$par,se,tv,pval)

        # colnames(results)<-c("b","se","t","p")
        # rownames(results) <- c(colnames(xmat),"dispersion")
        # results <- results[-nrow(results),] #ditch the dispersion term
        # },silent=TRUE
      if(class(results) == "try-error")
      {
       results <- matrix(c(NA,NA,NA,NA),nrow=1,ncol=4)
       row.names(results) <- "bcnv"
      }
    }
 
    if (study %in% c("psy3_comb","pts1_pts1"))
    {
     results <- summary(glm(pheno-1 ~ C1 + C2 + C3 + C4 + C5 + LRR_SD + cnv, data=cnv_fam_subset2,family="binomial"))$coefficients
    }
 

    #maybe need an additional one for UKBB to capture batch and stuff covariates
    
    #save regression results, but only if the analysis has been performed
    if(row.names(results)[nrow(results)] != "cnv")
    { 
     resmat2[which(resmat2$study == study),2:5] <- rep(NA,4)
    } else {
     resmat2[which(resmat2$study == study),2:5] <- results[nrow(results),]
     }
   }
   

   outfile=paste(cnvparms[i,]$Analysis,"_",cnvparms[i,]$Type,"_",cnvtype,sep="")
   write.table(resmat2,file=paste("outputs_jun2022/","bystudy_",outfile,'.txt',sep=''),row.names=F)
   try ( { meta_res <- rma(yi=resmat2$beta,sei=resmat2$se,slab=resmat2$study,method="FE") 
   capture.output(meta_res,paste("outputs_jun2022/","meta_",outfile,'.txt',sep=''))
   write.table(t(meta_res[c("beta","se","zval","pval")]),paste("outputs_jun2022/","meta2_",outfile,'.txt',sep=''),row.names=F)
  },silent=TRUE)
  
  pdf(paste("outputs_jun2022/","meta_",outfile,'.pdf',sep=''),7,7)
   plot(meta_res)
  dev.off()
  
  pdf(paste("outputs_jun2022/","forest_",outfile,'.pdf',sep=''),7,7)
   forest(meta_res)
  dev.off()
  }
 }  
 
 #Export all results to a single file
 rm cnv_burden_genes_jun13_2022.txt
 for files in $(ls outputs_jun2022 | grep meta2 )
 do
  output=$( tail -n1 outputs_jun2022/$files)
  echo $files $output >> cnv_burden_genes_jun13_2022.txt
 done
 

