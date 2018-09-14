#for bmi
#

load("./SNP+p_bmi_c/posteriorSummary.RData")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(grid)
library(gridExtra)
library(ggplot2)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(Other)
genexprobe<-as.matrix(Other$UCSC_RefGene_Name)
rownames(genexprobe)<-rownames(Other)

#result<-process_posterior(files,bmi_names,genexprobe,1:5)

tmp<-data.frame(M1_G=100*result$varSNP1,M2_G=100*result$varSNP2,M3_G=100*result$varSNP3,M1_phi=100*result$varP1,M2_phi=100*result$varP2,M3_phi=100*result$varP3,sample=1:length(result$varSNP1))
tmp<-reshape(tmp,varying=1:6,sep = "_",direction = "long")
tmp$sample<-NULL
tmp$id<-NULL
tmp<-melt(tmp)
colnames(tmp)<-c("Group","variable","value")
tmp$mixVar<-ifelse(tmp$Group=="G" & tmp$variable =="M1","0.0001",0)
tmp$mixVar<-ifelse(tmp$Group=="G" & tmp$variable =="M2","0.001",tmp$mixVar)
tmp$mixVar<-ifelse(tmp$Group=="G" & tmp$variable =="M3","0.01",tmp$mixVar)
tmp$mixVar<-ifelse(tmp$Group=="phi" & tmp$variable =="M1","0.01",tmp$mixVar)
tmp$mixVar<-ifelse(tmp$Group=="phi" & tmp$variable =="M2","0.1",tmp$mixVar)
tmp$mixVar<-ifelse(tmp$Group=="phi" & tmp$variable =="M3","1",tmp$mixVar)
library(ggplot2)
library(reshape2)
lg<-ggplot()
lg <- lg 
lg <- lg + geom_violin(data=tmp,aes(x=mixVar,y=value,color=Group,size=1.2),scale="width")
lg <- lg + xlab("Variance of mixture")+ylab("PVE (%)")+ facet_grid(Group ~ .,labeller = label_parsed,scales="free")

lg <- lg + theme_light()
lg <- lg + theme(legend.position="bottom",axis.text.x = element_text(size=20,face="plain"),
                                                axis.text.y = element_text(size=20,face="plain"),  
                                                axis.title.x = element_text(size=25,face="plain",hjust = 0.01),
                                                axis.title.y = element_text(size=25,face="plain"),
                                                legend.text =element_text(size=20,face="plain"),
                                                strip.text.y = element_text(size = 20),
                                                  plot.title = element_text(size=30)
                 )
lg <- lg + scale_color_discrete(labels=c("Genetic","Methylation"),name="")
lg<- lg + scale_size_continuous(labels=c(""),name="",guide = "none")
lg<- lg + guides(colour = guide_legend(override.aes = list(size=1.2)))
lg <- lg + ggtitle("BMI")
lg <- arrangeGrob(lg, top = textGrob("a", x = unit(0, "npc")
                                     , y   = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=30, fontfamily="Helvetica")))

#plot of PVE95
tmp3<-data.frame(TOTAL=(100*result$VEepi),IP95=100*result$VE95)
tmp3<-melt(tmp3)
colnames(tmp3)<-c("Markers","PVE")
lq<-ggplot(tmp3, aes(x="", y=PVE,color=Markers,size=1.2))
lq<- lq + geom_violin(position=position_dodge(0),scale="width")
lq<- lq + labs(x="",y="PVE (%)")

lq<- lq + ylim(0,100)
lq<- lq + theme_light()
lq<- lq +scale_size_continuous(labels=c(""),name="",guide = "none")
lq<- lq +theme(legend.position="bottom",axis.text.x = element_text(size=20,face="plain"),
               axis.text.y = element_text(size=20,face="plain"),  
               axis.title.x = element_text(size=25,face="plain",hjust = 0.01),
               axis.title.y = element_text(size=25,face="plain"),
               legend.text =element_text(size=20,face="plain"),
               strip.text.y = element_text(size = 20),
               plot.title = element_text(size=30),
               legend.title=element_blank()
               )
lq<- lq + scale_color_manual(labels = c("All methylation markers", "Methylation markers with 95% IP"), values = c("blue", "red"))
lq<- lq + guides(colour = guide_legend(override.aes = list(size=1.2)))
lq <- arrangeGrob(lq, top = textGrob("b", x = unit(0, "npc")
                                     , y   = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=30, fontfamily="Helvetica")))

library(gridExtra)
gp<-grid.arrange(lg,lq,nrow=1)
ggsave("VE_bmig.pdf",grid.draw(gp),scale = 2)




#for smoking
#

library(ggplot2)
library(reshape2)
library(gridExtra)
library(readr)
#We load the posterior summary object.
load("./SNP+p_smk/posteriorSummary.RData")
#We load the annotation for  the Ilumina 450k array
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(Other)
genexprobe<-as.matrix(Other$UCSC_RefGene_Name)
rownames(genexprobe)<-rownames(Other)
#tmp<-getMappedEntrezIDs(rownames(genexprobe),array.type = "450k")
#result<-process_posterior(files,bmi_names,genexprobe,1:5)


# In tmp we get a data frame that takes the 
tmp<-data.frame(M1_G=100*result$varSNP1,M2_G=100*result$varSNP2,M3_G=100*result$varSNP3,M1_phi=100*result$varP1,M2_phi=100*result$varP2,M3_phi=100*result$varP3,sample=1:length(result$varSNP1))
tmp<-reshape(tmp,varying=1:6,sep = "_",direction = "long")
tmp$sample<-NULL
tmp$id<-NULL
tmp<-melt(tmp)
colnames(tmp)<-c("Group","variable","value")
tmp$mixVar<-ifelse(tmp$Group=="G" & tmp$variable =="M1","0.0001",0)
tmp$mixVar<-ifelse(tmp$Group=="G" & tmp$variable =="M2","0.001",tmp$mixVar)
tmp$mixVar<-ifelse(tmp$Group=="G" & tmp$variable =="M3","0.01",tmp$mixVar)
tmp$mixVar<-ifelse(tmp$Group=="phi" & tmp$variable =="M1","0.01",tmp$mixVar)
tmp$mixVar<-ifelse(tmp$Group=="phi" & tmp$variable =="M2","0.1",tmp$mixVar)
tmp$mixVar<-ifelse(tmp$Group=="phi" & tmp$variable =="M3","1",tmp$mixVar)

#Plot of mixtures and variance explained
tg <- ggplot()
tg <- tg + geom_violin(data=tmp,aes(x=mixVar,y=value,color=Group,size=1.2),scale="width") 
tg <- tg + xlab("Variance of mixture")
tg <- tg + ylab("PVE %")
tg <- tg + facet_grid(Group ~ .,labeller = label_parsed,scales="free")+scale_color_discrete(labels=c("Genetic","Methylation"),name="")
tg <- tg + theme_light()
tg <- tg + theme(legend.position="none",axis.text.x = element_text(size=20,face="plain"),
                 axis.text.y = element_text(size=20,face="plain"),  
                 axis.title.x = element_text(size=25,face="plain",hjust = 0.01),
                 axis.title.y = element_text(size=25,face="plain"),
                 legend.text =element_text(size=20,face="plain"),
                 strip.text.y = element_text(size = 20),
                 plot.title = element_text(size=30)
)
tg<- tg + scale_size_continuous(labels=c(""),name="",guide = "none")
tg <- tg + ggtitle("Smoking")
tg <- arrangeGrob(tg, top = textGrob("c", x = unit(0, "npc")
                                     , y   = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=30, fontfamily="Helvetica")))

#Plot of PVE total vs 95% IP
tmp3<-data.frame(TOTAL=(100*result$VEepi),IP95=100*result$VE95)
tmp3<-melt(tmp3)
colnames(tmp3)<-c("Markers","PVE")
tq<- ggplot(tmp3, aes(x="", y=PVE,color=Markers,size=1.2))
tq<- tq + geom_violin(position=position_dodge(0),scale="width")
tq<- tq + labs(x="",y="PVE (%)")
tq<- tq + scale_color_manual(labels = c("All epigenetic markers", "95% IP"), values = c("blue", "red"))
tq<- tq +ylim(0,100)
tq<- tq+theme_light()
tq<- tq +scale_size_continuous(labels=c(""),name="",guide = "none")
tq<- tq +theme(legend.position="none",axis.text.x = element_text(size=20,face="plain"),
               axis.text.y = element_text(size=20,face="plain"),  
               axis.title.x = element_text(size=25,face="plain",hjust = 0.01),
               axis.title.y = element_text(size=25,face="plain"),
               legend.text =element_text(size=20,face="plain"),
               strip.text.y = element_text(size = 20),
               plot.title = element_text(size=30),
               legend.title=element_blank()
)

tq <- arrangeGrob(tq, top = textGrob("d", x = unit(0, "npc")
                                     , y   = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=30, fontfamily="Helvetica")))

## we add the genomic coverage plot
library(gridExtra)
gp<-grid.arrange(tg,tq,nrow=1)
ggsave("VE_smoking.pdf",grid.draw(gp),scale = 2)

load("./SNP+p_bmi_c/posteriorSummary.RData")
result_bmi<-result
load("./SNP+p_smk/posteriorSummary.RData")
result_smk<-result
Trait<-c(rep("BMI",length(result_bmi$genesPerSample)),rep("smoking",length(result_smk$genesPerSample)))
GP=c(100*result_bmi$genesPerSample,100*result_smk$genesPerSample)
q<-ggplot()
q<-q + geom_density(data=data.frame(GP=GP,Trait<-Trait),aes(x=GP,fill=Trait),alpha=0.25)
q<-q + xlab("Prop. of probes mapped to genes (%)") + ggtitle("Prop. of probes in model mapped to genes")
q<-q + theme_light()
q<-q +theme(legend.position="right",axis.text.x = element_text(size=20,face="plain"),
             axis.text.y = element_text(size=20,face="plain"),  
             axis.title.x = element_text(size=25,face="plain",hjust = 0.01),
             axis.title.y = element_text(size=25,face="plain"),
             legend.text =element_text(size=20,face="plain"),
             strip.text.y = element_text(size = 20),
             plot.title = element_text(size=30),
             legend.title=element_blank()
             )
q<-q + scale_fill_manual( values = c("red","blue"))
q <- arrangeGrob(q, top = textGrob("e", x = unit(0, "npc")
                                   , y   = unit(1, "npc"), just=c("left","top"),
                                   gp=gpar(col="black", fontsize=30, fontfamily="Helvetica")))


gp<-grid.arrange(lg,lq,tg,tq,q,nrow=1,layout_matrix=rbind(c(1,2),c(3,4),c(5,5)))
#ggsave("VE.pdf",plot=grid.draw(gp),scale=2)
ggsave("VE.pdf",plot=grid.draw(gp),scale=3.2)