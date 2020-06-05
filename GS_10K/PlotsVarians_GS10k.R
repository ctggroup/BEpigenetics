#for bmi
#
#code ocean directory
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(grid)
library(gridExtra)
library(ggplot2)
source("post_process_gs_10k.R")
source('plot_post_variance.R')

load("./bmi_thinned/bmi_colnames_GS10k.RData")


data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(Other)
genexprobe<-as.matrix(Other$UCSC_RefGene_Name)
rownames(genexprobe)<-rownames(Other)


save_mean_model() <- function(chains, filename){
  load("./bmi_thinned/bmi_colnames_GS10k.RData")
  
  meanBetas <- colMeans(C_bmi[,grep("beta",colnames(C_bmi)),with=F])
  names(meanBetas) <- bmi_colnames
  meanMu <- mean(unlist(C_bmi[,"mu",with=F]))
  save(list=c("meanBetas","meanMu"),file = filename)
  
}


# we read the thinned chains for thinned chains of bmi with pc and cell counts
C1_bmi_pc <- fread(file = "./bmi_thinned/C1_bmi_pc.csv_thinned.csv")
C2_bmi_pc <- fread(file = "./bmi_thinned/C2_bmi_pc.csv_thinned.csv")
C3_bmi_pc <- fread(file = "./bmi_thinned/C3_bmi_pc.csv_thinned.csv")
C4_bmi_pc <- fread(file = "./bmi_thinned/C4_bmi_pc.csv_thinned.csv")

C_bmi_pc <- rbindlist( list(C1_bmi_pc, C2_bmi_pc, C3_bmi_pc, C4_bmi_pc))
if(any(grepl("gamma",colnames(C_bmi_pc))))
  C_bmi_pc$sigmaC <- 46*rowSums(C_bmi_pc[,grep("gamma",colnames(C_bmi_pc)),with=F]^2)

save_mean_model(C_bmi, "meanModelBMI_pc.RData")
bmi_result_pc<-thinned_process_posterior(C_bmi,bmi_colnames,genexprobe)
save(list="bmi_result_pc",file="posteriorSummaryBMI_pc.RData")

# we read the thinned chains for bmi
C1_bmi <- fread(file = "./bmi_thinned/C1bmi.csv_thinned.csv")
C2_bmi <- fread(file = "./bmi_thinned/C2bmi.csv_thinned.csv")
C3_bmi <- fread(file = "./bmi_thinned/C3bmi.csv_thinned.csv")
C4_bmi <- fread(file = "./bmi_thinned/C4bmi.csv_thinned.csv")

C_bmi <- rbindlist( list(C1_bmi,C2_bmi,C3_bmi,C4_bmi))
save_mean_model(C_bmi, "meanModelBMI.RData")
bmi_result<-thinned_process_posterior(C_bmi,bmi_colnames,genexprobe)
save(list="bmi_result",file="posteriorSummaryBMI.RData")

# we read the thinned chains for smoking
C1_sm <- fread(file = "./sm_thinned/C1_sm.csv_thinned.csv")
C2_sm <- fread(file = "./sm_thinned/C2_sm.csv_thinned.csv")
C3_sm <- fread(file = "./sm_thinned/C3_sm.csv_thinned.csv")
C4_sm <- fread(file = "./sm_thinned/C4_sm.csv_thinned.csv")

C_sm <- rbindlist( list(C1_sm,C2_sm,C3_sm,C4_sm))
save_mean_model(C_sm, "meanModelSM.RData")
sm_result<-thinned_process_posterior(C_sm,bmi_colnames,genexprobe)
save(list="sm_result", file="posteriorSummarySM.RData")

# we read the thinned chains for smoking with pc and cell counts
C1_sm_pc <- fread(file = "./sm_thinned/C1_sm_pc.csv_thinned.csv")
C2_sm_pc <- fread(file = "./sm_thinned/C2_sm_pc.csv_thinned.csv")
C3_sm_pc <- fread(file = "./sm_thinned/C3_sm_pc.csv_thinned.csv")
C4_sm_pc <- fread(file = "./sm_thinned/C4_sm_pc.csv_thinned.csv")

C_sm_pc <- rbindlist( list(C1_sm_pc, C2_sm_pc, C3_sm_pc, C4_sm_pc))
if(any(grepl("gamma",colnames(C_sm_pc))))
  C_sm_pc$sigmaC <- 46*rowSums(C_sm_pc[,grep("gamma",colnames(C_sm_pc)),with=F]^2)

save_mean_model(C_sm_pc, "meanModelSM_pc.RData")
sm_result_pc<-thinned_process_posterior(C_sm_pc,bmi_colnames,genexprobe)
save(list="sm_result_pc", file="posteriorSummarySM_pc.RData")

