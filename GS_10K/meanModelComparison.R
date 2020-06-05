load("meanModelBMI.RData")
bmi_effects <- meanBetas

load("meanModelBMI_pc.RData")
bmi_pc_effects <- meanBetas

load("meanModelSM.RData")
sm_effects <- meanBetas


load("meanModelSM_pc.RData")
sm_pc_effects <- meanBetas

fit.bmi <- lm(bmi_effects ~ bmi_pc_effects)
fit.sm <- lm(sm_effects ~ sm_pc_effects)

library(jtools)
export_summs(fit.bmi,fit.sm, model.names= c("BMI effects","smoking effects"), to.file = "latex", file.name = "lm_pc.tex")

