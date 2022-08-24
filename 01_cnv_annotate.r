
grep -E "spID|NM_" hg19_genelist_names | LC_ALL=C sort -k1b,1 | sed 's/#kgID/name/g' >  hg19_genelist_NM.txt

LC_ALL=C sort -k1b,1 hg19_genelist | sed 's/#name/name/g' > hg19_genelist_sorted.txt

R

library(data.table)
library(plyr)

#load kgXref and hg19 genelist (only NM_ protein coding genes)
 d1 <- fread('hg19_genelist_NM.txt',data.table=F)
 d2 <- fread('hg19_genelist_sorted.txt',data.table=F)

 d3 <- merge(d1,d2,by="name")

#clean up formatting
 d3$CHROM <- gsub("chr","",d3$chrom)
 d3[which(d3$CHROM =="X"),]$CHROM=23
 d3 <- subset(d3, CHROM %in% c(1:23))

#Get full span of genes. From PLINK manual: These gene lists were downloaded from UCSC table browser for 
#all RefSeq genes on July 24th 2008. Overlapping isoforms of the same gene were combined to form a single 
#full length version of the gene. Isoforms that didn't overlap were left as duplicates of that gene. 

 d4 <- subset(d3,select=c(geneSymbol,refseq,CHROM))
 d4 <- d4[!duplicated(d4$geneSymbol),]

 startstop <- function(x)
 {
  minstart <- min(x$txStart)
  maxend <- max(x$txEnd)
  return(c(minstart,maxend))
 }

#For each gene, sort the start and stop, get the full span
 test <- ddply(.data=d3,~geneSymbol,startstop)
 names(test)[c(2,3)] <- c("BP1","BP2")
 test[,c(2,3)] <- as.numeric(test[,c(2,3)] )

 dm <- merge(d4,test,by="geneSymbol")

# #Get entrez ids, necessary for gene-sets
# #from https://www.ncbi.nlm.nih.gov/gene/?term=homo+sapiens
# cat gene_result.txt  | cut -f3,6 > homosapien_gene_id.txt

# entrezids <- fread('homosapien_gene_id.txt',data.table=F)


library(org.Hs.eg.db)
library(AnnotationDbi)

hs <- org.Hs.eg.db
my.symbols <- dm$geneSymbol
geneslist <- select(hs, 
       keys = my.symbols,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")
 
 #if nothing is returned, this gives na... e.g. with WASH1.. don't use entrez ids for this... they're fucked.
  select(hs, 
       keys = "WASH1",
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")
 
 
 
names(geneslist) <- c("geneSymbol","entrezID")

geneslist[geneslist$geneSymbol %in% geneslist[duplicated(geneslist$geneSymbol),]$geneSymbol,] #A few have extra IDs.. just delete them

geneslist <- geneslist[!duplicated(geneslist$geneSymbol),]

dm2 <- merge(dm,geneslist,by="geneSymbol")

geneslist[w
dmexpA <- subset(dm2,select=c(CHROM,BP1,BP2,geneSymbol))
dmexpB <- subset(dm2,select=c(CHROM,BP1,BP2,entrezID))

write.table(dmexpA,file='glist-hg19_genesymbol.txt',quote=F,row.names=F)
write.table(dmexpB,file='glist-hg19_entrezid.txt',quote=F,row.names=F)
write.table(dm2,   file='glist-hg19_adam.txt',quote=T,row.names=F)

 #Load CNV data
 cnv_data <- fread('ptsd_10kb_10probes_freqs.cnv',data.table=F)
  cnv_data$FID_IID <- paste(cnv_data$FID,cnv_data$IID,sep="_")
 #calculate CNV span now because it has to be done a lot
  cnv_data$span <- cnv_data$BP2 - cnv_data$BP1 
 
 
#Assign all CNVs to genes (if corresponding)
#can be in an expanded format or the other way around

#implies that x column 2 is chr, column 3 is bp1, column4 is bp2
olaptest <- function(x,genematrix)
{
 #Only use genes that are on the same chr
 genematrix$CHROM <- as.numeric( genematrix$CHROM)
 genemat <- subset(genematrix, CHROM==as.numeric(x["CHR"]))
 #CNV data positions
 a <- as.numeric(x["BP1"])
 b  <- as.numeric(x["BP2"])
 w <- as.numeric(genemat$BP1)
 z <- as.numeric(genemat$BP2)
 

 statement1 <- which(b < w )
 statement2 <- which(z < a)
 

 notaligned <- sort(unique(c(statement1,statement2)))
  paste(genemat[-notaligned,]$geneSymbol,collapse=",")
 }

 cnv_data$gene_annots <- apply(cnv_data,1,olaptest,genematrix=dm2)
 #write.table(cnv_data,file="ptsd_10kb_10probes_freqs_annotated.cnv",row.names=F)
 
 
#For doing cnvGSA with bank's data, use bank's symbols
 olaptest2 <- function(x,genematrix)
{
 #Only use genes that are on the same chr
 genematrix$CHROM <- as.numeric( genematrix$CHROM)
 genemat <- subset(genematrix, CHROM==as.numeric(x["CHR"]))
 #CNV data positions
 a <- as.numeric(x["BP1"])
 b  <- as.numeric(x["BP2"])
 w <- as.numeric(genemat$BP1)
 z <- as.numeric(genemat$BP2)
 
 statement1 <- which(b < w )
 statement2 <- which(z < a)
 
 notaligned <- sort(unique(c(statement1,statement2)))
  paste(genemat[-notaligned,]$entrezID,collapse=",")
 }
 
 cnv_data$gene_annots_entrezid <- apply(cnv_data,1,olaptest2,genematrix=dm2)
 
 #Just include both ID forms, no harm here.
 write.table(cnv_data,file="ptsd_10kb_10probes_freqs_annotated.cnv",row.names=F)
 
#Now in reverse, assign all genes to CNVs, to make an expanded matrix for gene-based only results
#this one is slightly different, returns an expanded CNV matrix, where every gene has a row.

olaptest3 <- function(x,genematrix)
{
 #Only use genes that are on the same chr
 genematrix$CHROM <- as.numeric( genematrix$CHROM)
 genemat <- subset(genematrix, CHROM==as.numeric(x["CHR"]))
 #CNV data positions
 a <- as.numeric(x["BP1"])
 b  <- as.numeric(x["BP2"])
 w <- as.numeric(genemat$BP1)
 z <- as.numeric(genemat$BP2)
 

 statement1 <- which(b < w )
 statement2 <- which(z < a )
 

 notaligned <- sort(unique(c(statement1,statement2)))
 if(dim(genemat[-notaligned,])[1] > 0)
 {

  return(
   cbind(x,genemat[-notaligned,]$geneSymbol
   ))
  } else ( return (as.data.frame(x)))
 }
 
#NOTE: NA WILL BE COUNTED AS A GENE SYMBOL! REMOVE THESE BEFORE ANY ANALYSIS OF GENE SETS OR GENES THAT WORK WITH GENE SYMBOLS!
cnv_data_annot <- adply(.data=cnv_data, .margins=1,.fun=olaptest3,genematrix=dm2)
 write.table(cnv_data_annot,file="ptsd_10kb_10probes_freqs_annotated_fullrank.cnv",row.names=F)
 