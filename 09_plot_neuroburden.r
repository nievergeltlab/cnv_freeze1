library(data.table)


library(Hmisc)
library(fmsb)
library(data.table)
library(plotrix)
library(lmtest)

library(ggplot2)


d1 <- fread('cnv_psychiatric_individual_v2.csv',data.table=F)

d1$colorz <- NA

d1[which(d1$Type=="Deletion"),]$colorz <- "blue"
d1[which(d1$Type=="Duplication"),]$colorz <- "red"
d1$position <- 1:dim(d1)[1]

pdf('cnv_neuro_jun21_2022.pdf',7,7)
zp1 <- ggplot(d1, aes(x=position,colour = Type)) 
zp1 <- zp1 + geom_hline(yintercept = 1, colour = gray(1/2), lty = 2) 
zp1 <- zp1 + geom_linerange(aes(x = Name, ymin = CI_lower,
                                ymax = CI_upperT),
                            lwd = 1, position = position_dodge(width = 1/4))
zp1 <- zp1 + geom_pointrange(aes(x =Name, y = OR, ymin = CI_lower,
                                 ymax = CI_upperT),
                             lwd = 1/2, position = position_dodge(width = 1/4),
                             shape = 21, fill = "WHITE")
                             
zp1 <- zp1 + coord_flip() + theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
zp1 <- zp1 + ggtitle("Comparing several models") + ylim(0,12)

print(zp1)
dev.off()


