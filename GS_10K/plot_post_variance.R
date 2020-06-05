library(ggplot2)
library(reshape2)
#TODO still change the component labels given that bmi and smoking used different
plot_post_variance <- function(result, genexprobe, plot_title, components = c("0.0001", "0.001", "0.01", "0.01", "0.1", "1.0"),lettersL){
  tmp<-data.frame(M1_G=100*result$varSNP1,M2_G=100*result$varSNP2,M3_G=100*result$varSNP3,M1_phi=100*result$varP1,M2_phi=100*result$varP2,M3_phi=100*result$varP3,sample=1:length(result$varSNP1))
  tmp<-reshape(tmp,varying=1:6,sep = "_",direction = "long")
  tmp$sample<-NULL
  tmp$id<-NULL
  tmp<-melt(tmp)
  colnames(tmp)<-c("Group","variable","value")
  tmp$mixVar<-ifelse(tmp$Group=="G" & tmp$variable =="M1",components[1],0)
  tmp$mixVar<-ifelse(tmp$Group=="G" & tmp$variable =="M2",components[2],tmp$mixVar)
  tmp$mixVar<-ifelse(tmp$Group=="G" & tmp$variable =="M3",components[3],tmp$mixVar)
  tmp$mixVar<-ifelse(tmp$Group=="phi" & tmp$variable =="M1",components[4],tmp$mixVar)
  tmp$mixVar<-ifelse(tmp$Group=="phi" & tmp$variable =="M2",components[5],tmp$mixVar)
  tmp$mixVar<-ifelse(tmp$Group=="phi" & tmp$variable =="M3",components[6],tmp$mixVar)
  tmp$mixVar <- as.character(tmp$mixVar)
 
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
  lg <- lg + ggtitle(plot_title)
  lg <- arrangeGrob(lg, top = textGrob(lettersL[[1]], x = unit(0, "npc")
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
  lq <- arrangeGrob(lq, top = textGrob(lettersL[[2]], x = unit(0, "npc")
                                       , y   = unit(1, "npc"), just=c("left","top"),
                                       gp=gpar(col="black", fontsize=30, fontfamily="Helvetica")))
  
  list(lg=lg, lq=lq)  
}
