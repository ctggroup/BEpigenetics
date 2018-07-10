source("TableSummary.R")
load("../SNP+p_smk/posteriorSummary.RData")
result$PIP
table.summary.smk<-table.summary(result)
print(xtable(table.summary.smk),file="../SNP+p_smk/table_summary_smk.tex")
load("../SNP+p_bmi_c/posteriorSummary.RData")
table.summary.bmi<-table.summary(result)
print(xtable(table.summary.bmi),file="../SNP+p_bmi_c/table_summary_bmi.tex")