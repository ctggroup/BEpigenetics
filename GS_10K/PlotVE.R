library(grid)
library(gridExtra)
library(ggplot2)
source("post_process_gs_10k.R")
source('plot_post_variance.R')

load("posteriorSummaryBMI.RData")
plot_bmi <- plot_post_variance(bmi_result, genexprobe, "BMI (CpG + SNP)",c("0.00001","0.0001","0.001","0.0001","0.001","0.01"),lettersL=c("a","b"))

load("posteriorSummarySM.RData")
plot_sm <- plot_post_variance(sm_result, genexprobe, "Smoking(CpG + SNP)",lettersL=c("c","d"))



Trait<-c(rep("BMI",length(bmi_result$genesPerSample)),rep("smoking",length(sm_result$genesPerSample)))
GP=100*c(bmi_result$genesPerSample,sm_result$genesPerSample)
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



gp<-grid.arrange(plot_bmi$lg,plot_bmi$lq,plot_sm$lg,plot_sm$lq,q,nrow=1,layout_matrix=rbind(c(1,2),c(3,4),c(5,5)))
ggsave("VE.pdf",plot=grid.draw(gp),scale=3.2)


# now for the PC and cell counts
# 
load("posteriorSummaryBMI_pc.RData")
plot_bmi_pc <- plot_post_variance(bmi_result_pc, genexprobe, "BMI  ( CG + SNP + PC + CC)")

load("posteriorSummarySM_pc.RData")
plot_sm_pc <- plot_post_variance(sm_result_pc, genexprobe, "SMK ( CG + SNP + PC + CC)")

#we need to recover the sigmaC from the thinned posteriors



# we read the thinned chains for smoking with pc and cell counts
C1_sm_pc <- fread(file = "./sm_thinned/C1_sm_pc.csv_thinned.csv")
C2_sm_pc <- fread(file = "./sm_thinned/C2_sm_pc.csv_thinned.csv")
C3_sm_pc <- fread(file = "./sm_thinned/C3_sm_pc.csv_thinned.csv")
C4_sm_pc <- fread(file = "./sm_thinned/C4_sm_pc.csv_thinned.csv")

C_sm_pc <- rbindlist( list(C1_sm_pc, C2_sm_pc, C3_sm_pc, C4_sm_pc))

C_sm_pc <- rbindlist( list(C1_sm_pc, C2_sm_pc, C3_sm_pc, C4_sm_pc))
if(any(grepl("gamma",colnames(C_sm_pc))))
  C_sm_pc$sigmaC <- 46*rowSums(C_sm_pc[,grep("gamma",colnames(C_sm_pc)),with=F]^2)

summary(C_sm_pc$sigmaC)

