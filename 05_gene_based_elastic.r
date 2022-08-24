
##CNV data handling
 library(Matrix)
 library(data.table)
 library('glmnet')
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
  
  cnv_fam$pheno_rescale <- NA
  cnv_fam$pheno_rescale_resid <- NA
 #Have to rescale here and residualize to make predictions meaningful. not quite optimal.
 for (study in studylist)
 {
  cnv_fam[which(cnv_fam$studyname2 %in% study),]$pheno_rescale = cnv_fam[which(cnv_fam$studyname2 %in% study),]$pheno * harmonize_scales[which(harmonize_scales$study %in% study),]$conversion_factor
  mf3 <- formula( paste ("pheno_rescale ~ C1 +C2 + C3 +C4 +C5 + LRR_SD",sep="")) # +PC1 + PC2 + PC3    
  modeltype=gaussian()
  cnv_fam[which(cnv_fam$studyname2 %in% study),]$pheno_rescale_resid <- resid(glm(mf3, family=modeltype,data= cnv_fam[which(cnv_fam$studyname2 %in% study),],na.action="na.exclude"))       
  }
  
#Read in the CNV call data, with full rank for gene names (see 00_genelist_annotat.txt)
 genesA <- fread("ptsd_10kb_10probes_freqs_annotated_fullrank.cnv",data.table=F)
 names(genesA)[dim(genesA)[2]] <- "geneSymbol"

#'NA' is not actually a gene name. These must be excluded
 genesA <- subset(genesA, !is.na(geneSymbol) && geneSymbol != "NA"  )
 
#Create a unique identifier for each subject 
 genesA$FID_IID <- paste(genesA$FID,genesA$IID)

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
 
 #just do ukbb
  analyze_these <- which(cnv_fam$studyname2 == "ukbb_ukbb")
  cnv_fam_subset <- cnv_fam[analyze_these,]
  cv_these <- which(cnv_fam$studyname2 != "ukbb_ukbb")
  genes_matrix2 <- genes_matrix[analyze_these,]
  
  gmat_freq <- apply(genes_matrix2,2,sum)
  maf_filtered <- which(gmat_freq>=140) #0.001
  genes_matrix3 <- genes_matrix2[,maf_filtered]

 #do cv with remaining data as holdout, phenotype rescaled from 0-1, or just looped and meta'd...
 #phenotype has to be residualized for PCs and LRR SD?
  fit <- cv.glmnet(genes_matrix3, cnv_fam_subset$pheno_rescale_resid)
  
  r2 <- fit$glmnet.fit$dev.ratio[which(fit$glmnet.fit$lambda == fit$lambda.min)]
  
  
  fit$nzero  
  plot(fit, xvar = "lambda", label = TRUE)
  plot(fit, xvar = "dev", label = TRUE)

#which genes are nonzero coefficient. -1 at the end to remove the intercept term
good_values <- (as.numeric(coef(fit,s=fit$lambda.min)[-1]) != 0)
coef(fit,s=fit$lambda.min)[-1][good_values]

write.table(cbind(genelist[maf_filtered][good_values],coef(fit,s=fit$lambda.min)[-1][good_values]),file=paste("elasticnet_",gtype,"_genes.txt"),row.names=F)


coef(fit,s=fit$lambda.min)


  cpred <- predict(fit,genes_matrix[cv_these,maf_filtered],s="lambda.min", type="response")
  cor.test(cpred,cnv_fam[cv_these,]$pheno_rescale_resid)
 
  plot(cpred,cnv_fam[cv_these,]$pheno_rescale_resid)
  dev.off()

  fit2 <- glmnet(genes_matrix3, cnv_fam_subset$pheno)
  fit3 <- bigGlm(genes_matrix3, cnv_fam_subset$pheno)
  fit4 <- lm(cnv_fam_subset$pheno ~ as.matrix(genes_matrix3) )
  
 