source("TableSummary.R")
#library(xtable)
load("posteriorSummaryBMI.RData")
table.summary.bmi<-table.summary(bmi_result)
#print(xtable(table.summary.bmi),file="res/table_summary_bmi.tex")
library(kableExtra)
tmp <-kable(table.summary.bmi, "latex", booktabs = T,digits=2) %>%
  kable_styling(latex_options = "striped")
cat(tmp,file = "res/table_summary_bmi.tex")
load("posteriorSummarySM.RData")
table.summary.smk<-table.summary(sm_result)
kable(table.summary.smk, "latex", booktabs = T,digits=2) %>%
  kable_styling(latex_options = "striped") %>% cat(file="res/table_summary_smk.tex")

#print(xtable(table.summary.smk),file="/results/table_summary_smk.tex")