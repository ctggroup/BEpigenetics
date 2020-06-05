load("posteriorSummaryBMI.RData")
load("posteriorSummaryBMI_pc.RData")
load("posteriorSummarySM.RData")
load("posteriorSummarySM_pc.RData")

bmi95 <- bmi_result$PIP[bmi_result$PIP>0.95]
bmi_pc95 <- bmi_result_pc$PIP[bmi_result_pc$PIP>0.95]
bmi_result_pc$PIP[bmi_result_pc$PIP>0.95]
bmi_result$PIP[bmi_result$PIP>0.95]

setdiff(names(bmi95), names(bmi_pc95))
setdiff(names(bmi_pc95), names(bmi95))

sm95 <- sm_result$PIP[sm_result$PIP>0.95]
sm_pc95 <- sm_result_pc$PIP[sm_result_pc$PIP>0.95]
setdiff(names(sm95), names(sm_pc95))
setdiff(names(sm_pc95), names(sm95))

save(list = c("bmi95","sm95"), file= "pip95.RData")
EWAS_Catalog_03_07_2019 <- read_delim("EWAS_Catalog_03-07-2019.txt","\t",escape_double = FALSE, trim_ws = TRUE)

library(tidyverse)
bmi<- EWAS_Catalog_03_07_2019 %>% filter(CpG %in% names(bmi95)) 
prev_bmi <- bmi %>% filter(str_detect(Trait,"Body mass index")) %>% select(CpG) %>% unique() 

not_bmi <- bmi %>%  filter(CpG %in% setdiff(CpG,prev_bmi$CpG)) %>% select(Trait,CpG) 


smk<- EWAS_Catalog_03_07_2019 %>% filter(CpG %in% names(sm95)) 
prev_smk <- smk %>% filter(str_detect(str_to_upper(Trait),"SMOKING")) %>% select(CpG) %>% unique() 


not_smk <- smk %>%  filter(CpG %in% setdiff(CpG,prev_smk$CpG)) %>% select(Trait,CpG) 
unique(not_smk$CpG)
table(not_smk$Trait)
not_smk_gene <-  smk %>%  filter(CpG %in% setdiff(CpG,prev_smk$CpG)) %>% select(Gene) %>% unique()

bmi_not_in <- names(sm95) %in% EWAS_Catalog_03_07_2019$CpG

kable(topGSA(gometh(names(bmi95))), "latex", booktabs = T,digits=2) %>%
  kable_styling(latex_options = "striped") %>% cat(file="GO95_bmi.tex")

kable(topGSA(gometh(names(sm95))), "latex", booktabs = T,digits=2) %>%
  kable_styling(latex_options = "striped") %>% cat(file="GO95_smoking.tex")
