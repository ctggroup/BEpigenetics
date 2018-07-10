#for bmi
#

load("../SNP+p_bmi_c/posteriorSummary.RData")

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
library(ggplot2)
library(reshape2)
lg<-ggplot()
lg <- lg 
lg <- lg + geom_violin(data=tmp,aes(x=variable,y=value,color=Group),scale="width")
lg <- lg + xlab("Mixture")+ylab("PVE %")+ facet_grid(Group ~ .,labeller = label_parsed)
lg <- lg + scale_color_discrete(labels=c("Genetic","Methylation"),name="")
lg <- lg + theme(legend.position="bottom")
lg <- lg + theme_light()
lg <- lg + ggtitle("BMI")
lg <- arrangeGrob(lg, top = textGrob("a", x = unit(0, "npc")
                                     , y   = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=18, fontfamily="Helvetica")))


tmp3<-data.frame(TOTAL=(100*result$VEepi),IP95=100*result$VE95)
tmp3<-melt(tmp3)
colnames(tmp3)<-c("Markers","PVE")
lq<-ggplot(tmp3, aes(x="", y=PVE,color=Markers))
lq<- lq + geom_violin(position=position_dodge(0),scale="width")
lq<- lq + labs(x="")
lq<- lq + scale_color_manual(labels = c("All epigenetic markers", "95% IP"), values = c("blue", "red"))
lq<- lq + ylim(0,100)
lq<- lq + theme_light()
lq <- arrangeGrob(lq, top = textGrob("b", x = unit(0, "npc")
                                     , y   = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=18, fontfamily="Helvetica")))

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
load("../SNP+p_smk/posteriorSummary.RData")
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


#Plot of mixtures and variance explained
tg <- ggplot()
tg <- tg + geom_violin(data=tmp,aes(x=variable,y=value,color=Group),scale="width") 
tg <- tg + xlab("Mixture")
tg <- tg + ylab("PVE %")
tg <- tg + facet_grid(Group ~ .,labeller = label_parsed)+scale_color_discrete(labels=c("Genetic","Methylation"),name="")
tg <- tg + theme(legend.position="bottom")+theme_light()
tg <- tg + ggtitle("Smoking")
tg <- arrangeGrob(tg, top = textGrob("c", x = unit(0, "npc")
                                     , y   = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=18, fontfamily="Helvetica")))

#Plot of PVE total vs 95% IP
tmp3<-data.frame(TOTAL=(100*result$VEepi),IP95=100*result$VE95)
tmp3<-melt(tmp3)
colnames(tmp3)<-c("Markers","PVE")
tq<- ggplot(tmp3, aes(x="", y=PVE,color=Markers))
tq<- tq + geom_violin(position=position_dodge(0),scale="width")
tq<- tq + labs(x="")
tq<- tq + scale_color_manual(labels = c("All epigenetic markers", "95% IP"), values = c("blue", "red"))
tq<- tq +ylim(0,100)
tq<- tq+theme_light()
tq <- arrangeGrob(tq, top = textGrob("d", x = unit(0, "npc")
                                     , y   = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=18, fontfamily="Helvetica")))

## we add the genomic coverage plot
library(gridExtra)
gp<-grid.arrange(tg,tq,nrow=1)
ggsave("VE_smoking.pdf",grid.draw(gp),scale = 2)

load("~/Documents/Methylation/SNP+p_bmi_c/posteriorSummary.RData")
result_bmi<-result
load("~/Documents/Methylation/SNP+p_smk/posteriorSummary.RData")
result_smk<-result
Trait<-c(rep("BMI",length(result_bmi$genesPerSample)),rep("smoking",length(result_smk$genesPerSample)))
GP=c(100*result_bmi$genesPerSample,100*result_smk$genesPerSample)
q<-ggplot()
q<-q + geom_density(data=data.frame(GP=GP,Trait<-Trait),aes(x=GP,fill=Trait),alpha=0.25)
q<-q + xlab("Genome coverage % (methylation probes)")+ggtitle("Per sample genomic coverage")
q<-q + theme_light()
q<-q + scale_fill_manual( values = c("red","blue"))
q <- arrangeGrob(q, top = textGrob("d", x = unit(0, "npc")
                                   , y   = unit(1, "npc"), just=c("left","top"),
                                   gp=gpar(col="black", fontsize=18, fontfamily="Helvetica")))


gp<-grid.arrange(lg,lq,tg,tq,q,nrow=1,layout_matrix=rbind(c(1,2),c(3,4),c(5,5)))
#ggsave("VE.pdf",plot=grid.draw(gp),scale=2)
ggsave("VE.pdf",plot=grid.draw(gp),scale=3.2)