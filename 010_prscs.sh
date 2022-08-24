hna_ptsd_pcs_v4_aug3_2021.fuma.gz

#eur_ptsd_pcs_v4_aug3_2021.fuma.gz
 #European sumstats format
   SNP          A1   A2   BETA      P
 
 zcat TotalPCL_MVP_eur.gz | awk 'BEGIN{OFS=" "}{if (NR==1) {$1="SNP";$4="A1";$5="A2";$8="BETA";$9="P"}; print $1,toupper($4),toupper($5),$8,$10}'| grep -v -P "A T|T A|C G|G C"  > TotalPCL_MVP_eur.gz.prscx

 #This should have all linkable SNPs in it (found in most studies)
 zcat ../eur_ptsd_pcs_v4_aug3_2021.fuma.gz | awk '{if (NR>1) print $1,$3,0,$2,toupper($4),toupper($5)}' | grep -v -P "A T|T A|C G|G C" > eur_ptsd_pcs_v4_aug3_2021.fuma.bim
  
  
python PRScs-master/PRScs.py --ref_dir=/mnt/ukbb/adam/ptsd/prs_csx/ldblk_1kg_eur --bim_prefix=eur_ptsd_pcs_v4_aug3_2021.fuma \
--sst_file=TotalPCL_MVP_eur.gz.prscx --n_gwas=186689  \
--out_dir=mvp_prs_pred --phi=1e-2  

echo "CHR SNP BP A1 A2 BETA" > prs_header.txt

cat prs_header.txt mvp_prs_pred_pst_eff_a1_b0.5_phi1e-02_chr*.txt > mvp_prs_pred_pst_eff_a1_b0.5_phi1e-02_allchr.txt
awk '{print $2,$4,$6}' mvp_prs_pred_pst_eff_a1_b0.5_phi1e-02_allchr.txt > mvp_prs_pred_pst_eff_a1_b0.5_phi1e-02_allchr.txt.score

awk '{print $2}'  mvp_prs_pred_pst_eff_a1_b0.5_phi1e-02_allchr.txt >  mvp_prs_pred_pst_eff_a1_b0.5_phi1e-02_allchr.txt.snplist

 #extract the list of SNPs relevant to the PRS from MVP, and MVP PTSD subjects only
 for chr in {1..22}
 do
  ./plink2 --bed /mnt/ukbb/adam/ptsd/bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".bed \
 --bim /mnt/ukbb/adam/ptsd/bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".bim_rename \
 --fam /mnt/ukbb/adam/ptsd/bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".fam \
 --extract  mvp_prs_pred_pst_eff_a1_b0.5_phi1e-02_allchr.txt.snplist  \
 --keep /mnt/ukbb/adam/ptsd/cnv_aug_2021/phenotype_files/p2_ukbb_ukbb_eur.pheno \
 --make-bed --out ukb_mvp2/mvpsnps_"$chr"  
 done
 
  ls ukb_mvp2/* | grep bed | sed 's/.bed//g' | awk '{print $1".bed",$1".bim",$1".fam"}' > ukb_mvp.mergelist
 
 ../plink --merge-list ukb_mvp.mergelist --make-bed --out prscs_jun21_2022/pts_ukbb_mix_am-qc.hg19.ch.fl.bg

#Get PRS and un average it.
study=ukb
for chr in {1..22}
do
 #./plink2 --bfile ukb_mvp2/mvpsnps_"$chr"  --score  mvp_prs_pred_pst_eff_a1_b0.5_phi1e-02_allchr.txt.score --out prscspreds/ukbb_$chr
 sed 's/#//g' prscspreds/ukbb_"$chr".sscore | awk -v CHR=$chr '{if(NR==1) TSCORE="TSCORE"CHR; if(NR>1) TSCORE=$3*$5; print $1,$2, TSCORE}' >  prscspreds/ukbb2_"$chr".sscore
done

#take the product of the count and number of snps as the prs

R
library(plyr)
library(data.table)
########### COMBINE ALL MRS1 DATA #############

#read each file into a data frame with the same name
for (i in 1:22)
{
	assign(
		paste('ukbb_',i,sep=''), fread(paste('prscspreds/ukbb2_',i,'.sscore',sep=''), header=T, data.table=F)
		) 
}

#parse the text list of data frame names as a list of data frames
data_list <- eval( 
			parse( 
				text=paste(
					"list(", paste(grep("ukbb",ls(),value=TRUE), collapse=','), ")" 
					)
				)
			)


#combine all data frames by id_visit (won't work for subjects missing this variable!!!!!)
datA <- join_all(data_list,by=c("FID","IID"), type="left", match="first")
#sum over all score columns

datA$SCORE <- apply(datA[,-c(1,2)],1,sum)

write.table(datA,'prscspreds/ukbb_pred.profile',quote=F,row.names=F)


for study in ftca nss1 grac psy3 gali mrsc psy2 pts1 ukbb
do
 ../plink --bfile prscs_jun21_2022/pts_"$study"_mix_am-qc.hg19.ch.fl.bg --extract  mvp_prs_pred_pst_eff_a1_b0.5_phi1e-02_allchr.txt.snplist --score  mvp_prs_pred_pst_eff_a1_b0.5_phi1e-02_allchr.txt.score sum  --out prscspreds/"$study"_pred
done


R
library(data.table)
 library(metafor)
 library(missreg3)
 library(ordinal)
 
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
 resmat0 <- as.data.frame(matrix(nrow=length(studylist),ncol=8))
 names(resmat0) <-c ("beta_cnv","se_cnv","t_cnv","p_cnv","beta_prs","se_prs","t_prs","p_prs")
 resmat0$study <- studylist
 resmat <- merge(resmat0,harmonize_scales,by="study")
 
 resmatx <- resmat

 cnv_data <- fread('cnv_filtered/ptsd_10kb_10probes_implicateddel_freqs.cnv.cnv',data.table=F)
 
  
  #FID_IID to link data 
   cnv_data$FID_IID <- paste(cnv_data$FID,cnv_data$IID,sep="_")
  #calculate CNV span now because it has to be done a lot
   cnv_data$span <- cnv_data$BP2 - cnv_data$BP1 
  
   cnv_data_subset <- subset(cnv_data,TYPE==1)                              
  
   #Make a matrix to store cnv count per subject
    cnv_fam_subset <- cnv_fam 
   #The CNV count detection only gets applied to subjects who appear in the CNV file to save time. 
   #Everyone else is imputed with zero
    cnv_fam_subset$cnv_count <- 0
    cnv_fam_subset$cnv_span <- 0
    cnv_fam_subset$cnv_avg <- 0
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
 
   #which subjects have data?
   examine_subjects_list <- unique(cnv_data_subset$FID_IID)
    examine_subjects_rows <- which(cnv_fam$FID_IID %in% examine_subjects_list)
   #Get the parameters for these subjects
    cnv_fam_subset[examine_subjects_rows, c("cnv_count","cnv_span","cnv_avg")] <- t(sapply(cnv_fam_subset[examine_subjects_rows,]$FID_IID,num_cnv,dataframe=cnv_data_subset))
    
   #write a file to store this data,
  #  write.table(subset(cnv_fam_subset,select=c(FID,IID,cnv_count,cnv_span,cnv_avg)),file=paste("global_covar/",cnvparms[i,]$Analysis,"_",cnvparms[i,]$Type,sep=""),row.names=F)
   cnvtype="cnv_span"
   resmat2 <- resmat
  resmatx <- resmat
  resmat3 <- resmat
 for (study in studylist )
   {
      
        cnv_fam_subset2A <- subset(cnv_fam_subset,studyname2 == study)
     studyload <- substr(study,1,4)
    prs <- fread(paste('prscspreds/',studyload,'_pred.profile',sep=''),data.table=F)
    prs$SCALE <- scale(prs$SCORE)


    # d0 <- merge(prs,pcs,by=c("FID","IID"))
    # d1 <- merge(d0,phe,by=c("FID","IID"))

    # summary(lm(Current_PTSD_Continuous ~ C1+ C2+ C3+ C4+ C5+SCALE,data=d1))
    # summary(lm(Current_PTSD_Continuous ~ C1+ C2+ C3+ C4+ C5,data=d1))

  
    cnv_fam_subset2 <- merge(cnv_fam_subset2A,prs,by=c("FID","IID"))
    
    cnv_fam_subset2$cnv <- cnv_fam_subset2[,cnvtype]
    
  
    if (!(study %in% c("psy3_comb","pts1_pts1")))
    {
      #ordinal logit
       results <- try(summary(clm(as.factor(pheno) ~ C1 + C2 + C3 + C4 + C5 + LRR_SD + cnv + SCALE, data=cnv_fam_subset2))$coefficients , silent = TRUE)
       results2 <- try(summary(lm(pheno ~ C1 + C2 + C3 + C4 + C5 + LRR_SD + cnv + SCALE, data=cnv_fam_subset2))$coefficients , silent = TRUE)
       
       print(study)
       print(summary(lm(pheno ~ C1 + C2 + C3 + C4 + C5 + LRR_SD + cnv + SCALE, data=cnv_fam_subset2))$r.squared - summary(lm(pheno ~ C1 + C2 + C3 + C4 + C5 + LRR_SD + SCALE, data=cnv_fam_subset2))$r.squared)
       #CNV explain 0.00033
       print(summary(lm(pheno ~ C1 + C2 + C3 + C4 + C5 + LRR_SD + cnv + SCALE, data=cnv_fam_subset2))$r.squared - summary(lm(pheno ~ C1 + C2 + C3 + C4 + C5 + LRR_SD + cnv, data=cnv_fam_subset2))$r.squared)
       #PRS explains 0.00684196
       
    }
 
    if (study %in% c("psy3_comb","pts1_pts1"))
    {
     results <- summary(glm(pheno-1 ~ C1 + C2 + C3 + C4 + C5 + LRR_SD  + cnv + SCALE, data=cnv_fam_subset2,family="binomial"))$coefficients
     results2 <- results
    }
 
     rCNV <- summary(lm(pheno ~ C1 + C2 + C3 + C4 + C5 + LRR_SD + cnv + SCALE, data=cnv_fam_subset2))$r.squared - summary(lm(pheno ~ C1 + C2 + C3 + C4 + C5 + LRR_SD + SCALE, data=cnv_fam_subset2))$r.squared
     rPRS <- summary(lm(pheno ~ C1 + C2 + C3 + C4 + C5 + LRR_SD + cnv + SCALE, data=cnv_fam_subset2))$r.squared - summary(lm(pheno ~ C1 + C2 + C3 + C4 + C5 + LRR_SD + cnv, data=cnv_fam_subset2))$r.squared
      
      
     resmatx[which(resmat2$study == study),2] <- nrow(cnv_fam_subset2)
     resmatx[which(resmat2$study == study),3] <- rCNV
     resmatx[which(resmat2$study == study),4] <- CI.Rsq(rCNV,nrow(cnv_fam_subset2),9,level=0.95)[2]
     resmatx[which(resmat2$study == study),5] <- rPRS
     resmatx[which(resmat2$study == study),6] <- CI.Rsq(rPRS,nrow(cnv_fam_subset2),9,level=0.95)[2]
     
     resmatx[which(resmatx[,3] == 0),3]  <- NA
     resmatx[which(resmatx[,4] == 0),4]  <- NA
     resmatx[which(resmatx[,5] == 0),3]  <- NA
     resmatx[which(resmatx[,6] == 0),4]  <- NA     
     
    #maybe need an additional one for UKBB to capture batch and stuff covariates
    
    #save regression results, but only if the analysis has been performed
 
     resmat2[which(resmat2$study == study),2:5] <- results[nrow(results)-1,]
     resmat2[which(resmat2$study == study),6:9] <- results[nrow(results),]
     
      resmat3[which(resmat2$study == study),2:5] <- results2[nrow(results2)-1,]
     resmat3[which(resmat2$study == study),6:9] <- results2[nrow(results2),]    
     
   }
   
   ccstudies <- which(resmat3$study %in% c("psy3_comb","pts1_pts1")) #these will be transformed onto the lienar regression scale again..
   resmat3$beta_cnv_scale <- resmat3$beta_cnv *resmat3$conversion_factor
   resmat3$se_cnv_scale <- resmat3$se_cnv *resmat3$conversion_factor
   resmat3[ccstudies,]$beta_cnv_scale <- resmat3[ccstudies,]$beta_cnv_scale   *resmat3[ccstudies,]$prevalence*resmat3[ccstudies,]$minprevalence
   resmat3[ccstudies,]$se_cnv_scale <- resmat3[ccstudies,]$se_cnv_scale *resmat3[ccstudies,]$prevalence*resmat3[ccstudies,]$minprevalence
   
    resmat3$beta_prs_scale <- resmat3$beta_prs *resmat3$conversion_factor
   resmat3$se_prs_scale <- resmat3$se_prs *resmat3$conversion_factor
   resmat3[ccstudies,]$beta_prs_scale <- resmat3[ccstudies,]$beta_prs_scale   *resmat3[ccstudies,]$prevalence*resmat3[ccstudies,]$minprevalence
   resmat3[ccstudies,]$se_prs_scale <- resmat3[ccstudies,]$se_prs_scale *resmat3[ccstudies,]$prevalence*resmat3[ccstudies,]$minprevalence
   
     
   
   
   outfile=paste("cnvdeletions_plus_prs",cnvtype,sep="")
   write.table(resmat2,file=paste("outputs_jun2022/","bystudy_",outfile,'.txt',sep=''),row.names=F)
   meta_cnv <- rma(yi=resmat2$beta_cnv,sei=resmat2$se_cnv,slab=resmat2$study,method="FE")
   meta_prs <- rma(yi=resmat2$beta_prs,sei=resmat2$se_prs,slab=resmat2$study,method="FE")
   
   meta_cnv_beta <- rma(yi=resmat3$beta_cnv_scale,sei=resmat3$se_cnv_scale,slab=resmat3$study,method="FE")
   meta_prs_beta <- rma(yi=resmat3$beta_prs_scale,sei=resmat3$se_prs_scale,slab=resmat3$study,method="FE")
    
   
   meta_cnv_rsq <- rma(yi=resmatx[,3],sei=resmatx[,4],slab=resmat2$study,method="FE")
   meta_prs_rsq <- rma(yi=resmatx[,5],sei=resmatx[,6],slab=resmat2$study,method="FE")
    

   capture.output(meta_cnv,paste("outputs_jun2022/","meta_",outfile,'.txt',sep=''))
   write.table(t(meta_cnv[c("beta","se","zval","pval")]),paste("outputs_jun2022/","meta2_",outfile,'_cnv.txt',sep=''),row.names=F)
   capture.output(meta_prs,paste("outputs_jun2022/","meta_",outfile,'.txt',sep=''))
   write.table(t(meta_prs[c("beta","se","zval","pval")]),paste("outputs_jun2022/","meta2_",outfile,'_prs.txt',sep=''),row.names=F)
    
   #capture.output(meta_cnv,paste("outputs_jun2022/","meta_",outfile,'.txt',sep=''))
   write.table(t(meta_cnv_rsq[c("beta","se","zval","pval")]),paste("outputs_jun2022/","meta2_",outfile,'_cnvrsq.txt',sep=''),row.names=F)
   #capture.output(meta_prs,paste("outputs_jun2022/","meta_",outfile,'.txt',sep=''))
   write.table(t(meta_prs_rsq[c("beta","se","zval","pval")]),paste("outputs_jun2022/","meta2_",outfile,'_prsrsq.txt',sep=''),row.names=F)

  
  
  pdf(paste("outputs_jun2022/","meta_",outfile,'_cnv.pdf',sep=''),7,7)
   plot(meta_cnv)
  dev.off()
  
  pdf(paste("outputs_jun2022/","forest_",outfile,'_cnv.pdf',sep=''),7,7)
   forest(meta_cnv)
  dev.off()
  
   pdf(paste("outputs_jun2022/","meta_",outfile,'_prs.pdf',sep=''),7,7)
   plot(meta_prs)
  dev.off()
  
  pdf(paste("outputs_jun2022/","forest_",outfile,'_prs.pdf',sep=''),7,7)
   forest(meta_prs)
  dev.off()
   
  
  
  }
 }  
 



conda activa