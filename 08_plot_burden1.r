library(data.table)


library(Hmisc)
library(fmsb)
library(data.table)
library(plotrix)
library(lmtest)


d1 <- fread('cnv_psychiatric_burden_v2.csv',data.table=F) #[1:27,]


d1$colorz <- NA
d1[which(d1$Type =="ZEither"),]$colorz <- "black"
d1[which(d1$Type =="ZDuplication"),]$colorz <- "blue"
d1[which(d1$Type =="ZDeletion"),]$colorz <- "red"
d1$Name4 <- d1[,1]



library(ggplot2)
pdf('cnv_burden_jun21_2022.pdf',7,4)
zp1 <- ggplot(d1, aes(color = Type)) + scale_color_manual(values = c( "ZEither"= 'black', "ZDuplication" = "#00BFC4" ,  "ZDeletion" = "#F8766D" ))
zp1 <- zp1 + geom_hline(yintercept = 1, colour = gray(1/2), lty = 2)
zp1 <- zp1 + geom_linerange(aes(x = Name4, ymin = CI_lower,
                                ymax = CI_upper),
                            lwd = 1, position = position_dodge(width = 1/4))
zp1 <- zp1 + geom_pointrange(aes(x = Name4, y = OR, ymin = CI_lower,
                                 ymax = CI_upper),
                             lwd = 1/2, position = position_dodge(width = 1/4),
                             shape = 21, fill = "WHITE")
                             
zp1 <- zp1 + coord_flip() + theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
zp1 <- zp1 + ggtitle("Global burden test")  + ylim(0,3)

print(zp1)
dev.off()



          
