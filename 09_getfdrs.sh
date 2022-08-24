

d1 <- fread('overallsets.txt',data.table=F)
write.table(p.adjust(d1$V1,"BH"),file="fdr_overallsets.txt",row.names=F)

d1 <- fread('dupsets.txt',data.table=F)
write.table(p.adjust(d1$V1,"BH"),file="fdr_dupsets.txt",row.names=F)

d1 <- fread('delsets.txt',data.table=F)
write.table(p.adjust(d1$V1,"BH"),file="fdr_delsets.txt",row.names=F)

library(data.table)
d1 <- fread('pvlist.txt',data.table=F)
write.table(p.adjust(d1$V1,"BH"),file="fdr_pvlist.txt",row.names=F)



library(data.table) 
d1 <- fread('gene_sets_fdr.csv',data.table=F)
d1$fdrboth <- p.adjust(d1$pboth,"BH")
d1$fdrdel <- p.adjust(d1$pdel,"BH")
d1$fdrdup <- p.adjust(d1$pdup,"BH")

write.csv(d1,file="gene_sets_fdr.fdrs.csv",row.names=F)



d1 <- fread('overallsets.txt',data.table=F)
write.table(p.adjust(d1$V1,"BH"),file="fdr_overallsets.txt",row.names=F)

d1 <- fread('dupsets.txt',data.table=F)
write.table(p.adjust(d1$V1,"BH"),file="fdr_dupsets.txt",row.names=F)

d1 <- fread('delsets.txt',data.table=F)
write.table(p.adjust(d1$V1,"BH"),file="fdr_delsets.txt",row.names=F)
