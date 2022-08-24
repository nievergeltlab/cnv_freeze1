#CNV GSA package
 library(cnvGSA)
 library(data.table)
 library(org.Hs.eg.db)
 library(AnnotationDbi)
 library(Matrix)
 library(metafor)
#Load banks lists, annotate them to geneSymbol
 load('genesets/gsMain_PGC_2021.RData')

  #His thing needs to be unwrapped into a line style matrix, then annotated
  geneSetMat1 <- as.data.frame(matrix(nrow=length(unlist(gsMain)),ncol=3))
  names(geneSetMat1) <- c("geneset","entrezID","geneSymbol")

  indexer=1
  for (i in names(gsMain))
  {
   #how many items are in this set?
   gsunlist <- unlist(gsMain[i])
   gslength <- length(gsunlist)

   geneSetMat1[c(indexer:c(indexer+gslength-1)),1] <- i
   geneSetMat1[c(indexer:c(indexer+gslength-1)),2] <- as.vector(gsunlist)
   indexer=indexer+gslength
  }

  hs <- org.Hs.eg.db

  geneSetMat1$geneSymbol <- AnnotationDbi::select(hs, 
   keys = geneSetMat1$entrezID,
   columns = c("ENTREZID", "SYMBOL"),
   keytype = "ENTREZID")$SYMBOL
  
  #write.table(geneSetMat,file="genesets/gene_sets_withsymbols.txt",row.names=F)

#write out significan set members: this is here for convenience
# siglist <- scan(what=character())
# Neurof_GoNeuronProj
# Neurof_UnionStringent
# PhHs_NervSys_All
# scRNA_DEGs_ExN
# scRNA_Expressed_ExM
# scRNA_Expressed_ExM_U


 # gesp <- subset(geneSetMat, geneset %in% siglist)
   # write.table(gesp,file="genesets/significant_ gene_sets_withsymbols.txt",row.names=F)

 #modify this to have a minimum PREVALNCE FITLER!
 
#Read phenotype and covariates (includes complete sample list)
 cnv_fam <- fread('ptsd_10kb_10probes.cov',data.table=F)
  cnv_fam$FID_IID <- paste(cnv_fam$FID,cnv_fam$IID,sep="_")
 
 #There is a mistake - psy3_comb should be case control. I use the dx phenotype instead
  cnv_fam[which(cnv_fam$studyname2 == "psy3_comb"),]$pheno <- cnv_fam[which(cnv_fam$studyname2 == "psy3_comb"),]$phenoPTS

 #PPDS is rounded to nearest whole number (due to harmonizing, I created some decimals..)
  cnv_fam[which(cnv_fam$studyname2 == "nss1_ppds"),]$pheno <- round(cnv_fam[which(cnv_fam$studyname2 == "nss1_ppds"),]$pheno)

 #GUTS and NHS2 will be analyzed together
  cnv_fam[which(cnv_fam$studyname2 == "psy2_guts"),]$studyname2 <- "psy2_nhsy"

#List all studies analyzed
 studylist <- unique(cnv_fam$studyname2)
 
#Load effect size conversion sheet to harmonize PTSD
 harmonize_scales <-fread('effect_scale_conversion.txt',data.table=F)[,-1]
  names(harmonize_scales)[1] <- "study"
  harmonize_scales$prevalence <- harmonize_scales$N_Cases/harmonize_scales$N
  harmonize_scales$minprevalence <- harmonize_scales$N_Controls/harmonize_scales$N
  
#Read in the CNV call data, with full rank for gene names (see 00_genelist_annotat.txt)
 genesA <- fread("ptsd_10kb_10probes_freqs_annotated_fullrank.cnv",data.table=F)
 names(genesA)[dim(genesA)[2]] <- "geneSymbol"

#'NA' is not actually a gene name. These must be excluded
 genesA <- subset(genesA, !is.na(geneSymbol))
 
#To run without implicated, you need to have a global covariate for CNV, but only analyze the data itself with the implicated removed..
#Comment as needed.
 ndd_to_remove <- fread('ndd_regions.txt',data.table=F)
 #genesA <-subset(genesA
 #Fuck it, just do it on a loop
 #Cut out anything within 1 million bases of these regions.
 # genesAT <- genesA
 # for (J in 1:dim(ndd_to_remove)[1])
 # {
  # genesAT <- subset(genesAT,!(CHR == ndd_to_remove[J,]$chr & BP1 >= ndd_to_remove[J,]$bp1 -1000000 &  BP2 <=  ndd_to_remove[J,]$bp2 + 1000000))
 # }
 
 genesA <- genesAT
#Create a unique identifier for each subject 
 genesA$FID_IID <- paste(genesA$FID,genesA$IID)


for (gtype in  c("all","dup","del")) # all",
{
 #Subset to dups and dels if needed. + #Load in the global burden data (avg CNV size + CNV count)
  if(!(gtype %in% c("dup","del")))
  {
   genes <- genesA
   cnvGlobal <- fread('global_covar/Genome-wideBurden_Both',data.table=F)
   cnvGlobalPsych <- fread('global_covar/Genome-wideBurdenImplicated_Both',data.table=F)
   
  }
  if(gtype == "dup")
  {
   genes <- subset(genesA,TYPE==3)
   cnvGlobal <- fread('global_covar/Genome-wideBurden_Insertion',data.table=F)
   cnvGlobalPsych <- fread('global_covar/Genome-wideBurdenImplicated_Insertion',data.table=F)
   

  }
  if(gtype == "del")
  {
   genes <- subset(genesA,TYPE==1)
   cnvGlobal <- fread('global_covar/Genome-wideBurden_Deletion',data.table=F)
   cnvGlobalPsych <- fread('global_covar/Genome-wideBurdenImplicated_Deletion',data.table=F)
  

  }
  geneSetMat <- subset(geneSetMat1) # if filtering needs to eb performed.
   #Filter out all genes that overlap NDD regions... must identify these...
 #List all gene names. Gene matrix columns will be ordered the same as this (i.e. this is the ordered set of column names)
  genelist <- unique(genes$geneSymbol)
 #List of all subject names that are in your CNV analysis.
  #Gene matrix rows will be ordered the same as this (i.e. this is the ordered set of rows)
  subjects1 <-  merge(cnv_fam,cnvGlobal,by=c("FID","IID"),suffixes=c("","_extrastuff")) 
  subjects <- merge(subjects1,cnvGlobalPsych,by=c("FID","IID"),suffixes=c("","_psych")) 

 #Create a unique identifier for each subject
  subjects$FID_IID <- paste(subjects$FID,subjects$IID)
 #Assign coordinate entries to all subjects of the gene matrix, to populate a sparse matrix
  genes$coord1 <- match(genes$FID_IID, subjects$FID_IID)
  genes$coord2 <- match(genes$geneSymbol, genelist)
 #Create a sparse matrix to hold the gene data
  genes_matrix <- sparseMatrix(i=genes$coord1,j=genes$coord2,dims=c(nrow(subjects),length(genelist)))
 #1 row per subject, 1 column per gene. Binary matrix for has cnv versus not.
  
 #Create a matrix to store all meta-analyses 
  meta_results <- as.data.frame(matrix(nrow=length(unique(geneSetMat$geneset)), ncol=5))
  names(meta_results) <- c("geneSet","B","SE","Z","P")
  meta_results$geneSet <- unique(geneSetMat$geneset)
    
  #Analysis: Loop over every gene-set, then every study
  for (i in unique(geneSetMat$geneset)) # 1:ceiling((dim(genes_matrix)[2]/increment)))
  {
   #list of geneSymbols to analyze
    geneSetMat_analyze <- unique(subset(geneSetMat,geneset==i)$geneSymbol) #unique is on here just in case the same symbol appears twice in the set.
    geneSetMat_analyze_colnames <- match(geneSetMat_analyze, genelist)
    
    table(is.na(geneSetMat_analyze_colnames))
   #Comput how many are missing from the set 
   
   #THen dump them
    geneSetMat_analyze_colnames <- na.omit(geneSetMat_analyze_colnames)
    
   #Identify the matrix level position of every gene to be analyzed    
    genes_subset <- as.matrix(genes_matrix[,geneSetMat_analyze_colnames]) 
    genelist_subset <- genelist[geneSetMat_analyze_colnames] #gene list names
    
   #Now the matrix has to be collapsed, into the count of genic CNV in the data - this is the sum over each row (subject)
    genes_subset_count <- apply(genes_subset,1,sum) #sum 
  
   #Create a matrix to store regression coefficients. Should probably store some count info in here..
   
    results_matrix <- as.data.frame(matrix(nrow=length(studylist),ncol=5))
    names(results_matrix) <- c("study","beta","se","z","P")
    results_matrix$study <- studylist
    resmat2 <- merge(results_matrix,harmonize_scales,by="study")
   #Loop over studies
    for (study in studylist )
    {
     
     #Chop out the relevant people. notice that the typical cnv_fam is now replaced by subjects since we merged inglobal covariates earlier, 
     #which changed the fundamental subject order of the cnv matrix
     analyze_these <- which(subjects$studyname2 == study)
     cnv_fam_subset <- subjects[analyze_these,]
     genes_subset_count_study <- genes_subset_count[analyze_these]
     
     #CUT OUT NDD CARRIERS
      # no_ndd <- which(cnv_fam_subset$cnv_count_psych == 0) 
      # cnv_fam_subset <- cnv_fam_subset[no_ndd,]
      # genes_subset_count_study <- genes_subset_count_study[no_ndd]


     nsub <- dim(cnv_fam_subset)[1]
    
    #NO NDD VARIABLE DUE TO THESE BEING CUT OUT
     #results <- summary(lm(pheno ~ C1 + C2 + C3 + C4 + C5 + LRR_SD + cnv_count + cnv_avg + cnv_count_psych + cnv_avg_psych  + genes_subset_count_study, data=cnv_fam_subset))$coefficients
     results <- summary(lm(pheno ~ C1 + C2 + C3 + C4 + C5 + LRR_SD + cnv_count + cnv_avg  + genes_subset_count_study, data=cnv_fam_subset))$coefficients
     
     if (study %in% c("psy3_comb","pts1_pts1"))
     {
     # results <- summary(glm(pheno-1 ~ C1 + C2 + C3 + C4 + C5 + LRR_SD +  cnv_count + cnv_avg  + cnv_count_psych + cnv_avg_psych  + genes_subset_count_study, data=cnv_fam_subset,family="binomial"))$coefficients
      results <- summary(glm(pheno-1 ~ C1 + C2 + C3 + C4 + C5 + LRR_SD +  cnv_count + cnv_avg   + genes_subset_count_study, data=cnv_fam_subset,family="binomial"))$coefficients
      
     }
     
     if(row.names(results)[nrow(results)] != "genes_subset_count_study") 
     {
      resmat2[which(resmat2$study == study),2:5] <- rep(NA,4)
     } else {   
       resmat2[which(resmat2$study == study),2:5] <- results[nrow(results),]
       }
      
   
     }
     
     ccstudies <- which(resmat2$study %in% c("psy3_comb","pts1_pts1")) #these will be transformed onto the lienar regression scale again..
     resmat2$beta_scale <- resmat2$beta *resmat2$conversion_factor
     resmat2$se_scale <- resmat2$se *resmat2$conversion_factor
     resmat2[ccstudies,]$beta_scale <- resmat2[ccstudies,]$beta_scale   *resmat2[ccstudies,]$prevalence*resmat2[ccstudies,]$minprevalence
     resmat2[ccstudies,]$se_scale <- resmat2[ccstudies,]$se_scale *resmat2[ccstudies,]$prevalence*resmat2[ccstudies,]$minprevalence
     
     
     write.table(resmat2,file=paste("gene_setR/",gtype,"_",i,"_allstudies_junrev.txt",sep=""),row.names=F)
     res <- rma(yi=resmat2$beta_scale,sei=resmat2$se_scale,slab=resmat2$study,method="FE")
     meta_results[which(meta_results$geneSet ==i),2:5] <- c(res$beta,res$se,res$zval,res$pval) 
    }

    write.table(meta_results,file=paste("gene_setR/",gtype,"_allmetas_junrev.txt",sep=""),row.names=F)
   
}








