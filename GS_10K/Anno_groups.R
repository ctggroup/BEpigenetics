#for bmi
#
#code ocean directory
library(tidyverse)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(grid)
library(gridExtra)
library(ggplot2)
source("post_process_gs_10k.R")
source('plot_post_variance.R')

load("./bmi_thinned/bmi_colnames_GS10k.RData")

sigmas_names  <- c("age and sex",
                   "Adipose Subcutaneous",
                   "Adipose Visceral (Omentum)",
                   "Adrenal Gland",
                   "Artery Aorta",
                   "Artery Tibial",
                   "Brain Amygdala",
                   "Brain Ant. cing. cortex",
                   "Brain Caudate",
                   "Brain Cerebellar Hemisphere",
                   "Brain Cerebellum",
                   "Brain Cortex",
                   "Brain Frontal Cortex",
                   "Brain Hippocampus",
                   "Brain Hypothalamus",
                   "Brain Nucleus accumbens",
                   "Brain Putamen",
                   "Brain Spinal cord",
                   "Brain Substantia nigra",
                   "Breast Mammary Tissue",
                   "Cells EBV--transf. lymphocytes",
                   "Cells Transformed fibroblasts",
                   "Colon Sigmoid",
                   "Colon Transverse",
                   "Esophagus Mucosa",
                   "Esophagus Muscularis",
                   "Heart Atrial Appendage",
                   "Heart Left Ventricle",
                   "Kidney Cortex",
                   "Liver",
                   "Lung",
                   "Minor Salivary Gland",
                   "Muscle Skeletal",
                   "Nerve Tibial",
                   "Ovary" ,
                   "Pancreas", 
                   "Pituitary", 
                   "Prostate",
                   "Skin Not Sun Exposed",
                   "Skin Sun Exposed",
                   "Small Intestine Terminal Ileum",
                   "Spleen",
                   "Stomach",
                   "Testis", 
                   "Thyroid",
                   "Uterus",
                   "Whole Blood",
                   "others",
                   "no information",
                   "not mapped",
                   "SNP")
group_ids <- data.frame(groups = sigmas_names, group_id = 0:50)
groupxprobe <- fread("./tissueSpecN1/tissueSpecN1_groups_only.csv")
groupxprobet <- table(groupxprobe) 
names(groupxprobet) <- sigmas_names
groupxprobet <- as.matrix(groupxprobet)
groupxprobet <- t(groupxprobet)
groupxprobet <- as_tibble(groupxprobet)

colnames(groupxprobe) <- "group_id"
probexname <- merge(groupxprobe, group_ids)

# we read the thinned chains for thinned chains of bmi with pc and cell counts
C1_bmi <- fread(file = "./C_tsn1/C1_tsn1.csv_thinned.csv")
C2_bmi  <- fread(file = "./C_tsn1/C2_tsn1.csv_thinned.csv")
C3_bmi <- fread(file = "./C_tsn1/C3_tsn1.csv_thinned.csv")
C4_bmi <- fread(file = "./C_tsn1/C4_tsn1.csv_thinned.csv")



C_bmi <- rbindlist( list(C1_bmi, C2_bmi, C3_bmi, C4_bmi),fill=T)

sigmas_bmi <- C_bmi[,grep("sigma",colnames(C_bmi)) ,with=F]
sigmas_bmi <- as_tibble(sigmas_bmi)
sigmas_bmi[sigmas_bmi == Inf] <- 0
sigmas_bmi$`sigmaG[1]` <- 0
sigmas_bmi <- sigmas_bmi %>% mutate(sigma_tot = rowSums(.)) 
sigmas_bmi_norm <- sigmas_bmi %>% mutate_all(funs(./sigma_tot)) %>% select(-sigma_tot) 




colnames(sigmas_bmi_norm)[-1] <-sigmas_names
sigmas_bmi_norm <- sigmas_bmi_norm  %>% select(-sigmaE)


groups_plot_bmi <- sigmas_bmi_norm %>% gather( variable, value) %>% ggplot() + geom_violin(aes(variable,value))+theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
#lets try with a ridge plot instead
sigmas_bmi_norm %>% gather( variable, value) %>% filter(variable != "SNP") %>% ggplot() + geom_density_ridges_gradient(aes(value,variable, fill= stat(x))) + theme_light()


ggsave(filename="groups_bmi.pdf")

sigmas_bmi_norm %>% mutate_all(funs(./groupxprobet$.))  %>% gather( variable, value) %>% ggplot() + geom_violin(aes(variable,value))+theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("bminormalised to probes per group")
ggsave(filename= "groups_bmi_norm.pdf")

sqrdnorm <- function(x){
  sum(x^2)
}

plot_groups_mixture <- function( chains){
  betas <- chains[,grep("beta",colnames(chains)) ,with=F] 
  comp <- chains[,grep("comp",colnames(C_bmi)) ,with=F]
  mask_1 <- ifelse(comp==1,1,0)
  mask_2 <- ifelse(comp==2,1,0)
  mask_3 <- ifelse(comp==3,1,0)
  betas_1 <- mask_1 * betas
  betas_2 <- mask_2 * betas
  betas_3 <- mask_3 * betas
  betas_1 <- as_tibble(cbind(probes = colnames(bmi_colnames), t(betas_1)))
  betas_2 <- as_tibble(cbind(probes = colnames(bmi_colnames), t(betas_2)))
  betas_3 <- as_tibble(cbind(probes = colnames(bmi_colnames), t(betas_3)))
  betas_1$group <- probexname$groups
  betas_2$group <-  probexname$groups
  betas_3$group <-  probexname$groups
  VE_1 <- betas_1 %>%  group_by(group) %>%  summarise_all(sqrdnorm)
  VE_2 <- betas_2 %>%  group_by(group) %>%  summarise_all(sqrdnorm)
  VE_3 <- betas_3 %>%  group_by(group) %>%  summarise_all(sqrdnorm)
  VE_1_long <- VE_1 %>% gather(key = sample,value =  VE, -group)
  VE_1_long$mixture <- 0.0001
  VE_2_long <- VE_2 %>% gather(key = sample,value =  VE, -group)
  VE_2_long$mixture <- 0.001
  VE_3_long <- VE_3 %>% gather(key = sample,value =  VE, -group)
  VE_3_long$mixture <- 0.01 
  VE_mix <- bind_rows(VE_1_long,VE_2_long,VE_3_long)
  VE_mix <- VE_mix %>% mutate_at(vars(group),funs(factor))
  VE_plot <- VE_mix %>% ggplot() #+ geom_violin(aes(x=group, y= VE)) + facet_grid(mixture ~ .) +theme(axis.text.x = element_text(angle = 90, hjust = 1))
  VE_mix
}

gxmplot_bmi <- plot_groups_mixture(C_bmi)
tmp <- gxmplot_bmi %>% filter(group != "SNP") %>% inner_join(gxmplot_bmi %>% filter(group != "SNP") %>% group_by(group) %>% summarise( meanVE = mean(VE)  )) %>% arrange(meanVE) 
gxmplot_bmi <- tmp %>% filter(group!= "not mapped") %>% filter(group!= "no information") %>%  dplyr::rename(Tissue = group) %>% ggplot() + geom_density_ridges(aes(VE,Tissue)) + facet_grid( . ~ mixture, scales ="free_x") + theme_ridges() + xlab("Variance explained")

# we read the thinned chains for thinned chains of bmi with pc and cell counts
C1_smk <- fread(file = "./C_tsn1/C1smk_tsn1.csv_thinned.csv")
C2_smk  <- fread(file = "./C_tsn1/C2smk_tsn1.csv_thinned.csv")
C3_smk <- fread(file = "./C_tsn1/C3smk_tsn1.csv_thinned.csv")
C4_smk <- fread(file = "./C_tsn1/C4smk_tsn1.csv_thinned.csv")


C_smk <- rbindlist( list(C1_smk, C2_smk, C3_smk, C4_smk),fill=T)
sigmas_smk <- C_smk[,grep("sigma",colnames(C_smk)) ,with=F]
sigmas_smk <- as_tibble(sigmas_smk)
sigmas_smk[sigmas_smk == Inf] <- 0
sigmas_smk$`sigmaG[1]` <- 0

sigmas_smk <- sigmas_smk %>% mutate(sigma_tot = rowSums(.)) 
sigmas_smk_norm <- sigmas_smk %>% mutate_all(funs(./sigma_tot)) %>% select(-sigma_tot) 
colnames(sigmas_smk_norm)[-1] <-sigmas_names
sigmas_smk_norm <- sigmas_smk_norm %>% select(-sigmaE)


groups_plot_smk <- sigmas_smk_norm %>% gather( variable, value) %>% ggplot() + geom_violin(aes(variable,value))+theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave(filename="groups_smk.pdf")

sigmas_smk_norm %>% mutate_all(funs(./groupxprobet$.)) %>% gather( variable, value) %>% ggplot() + geom_violin(aes(variable,value))+theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("smoking normalised by probes per group")
ggsave(filename = "groups_smk_norm.pdf")

gxmplot_smk <- plot_groups_mixture(C_smk)
tmp_smk <- gxmplot_smk %>% filter(group != "SNP") %>% inner_join(gxmplot_smk %>% filter(group != "SNP") %>% group_by(group) %>% summarise( meanVE = mean(VE)  )) %>% arrange(meanVE) 
tmp_smk$group <- factor(tmp_smk$group, levels= unique(tmp_smk$group))
gxmplot_smk <- tmp_smk %>% filter(group!= "not mapped") %>% filter(group!= "no information") %>%  dplyr::rename(Tissue = group) %>% ggplot() + geom_density_ridges(aes(VE,Tissue)) + facet_grid( . ~ mixture, scales ="free_x") + theme_ridges() + xlab("Variance explained") 

save(list = "gxmplot_smk", file = "gxmplot_smk.RData")

source("./plotGenes.R")

#Here lets do plot aesthethics
#
groups_plot_bmi_sub <- groups_plot_bmi + theme_light() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) 
groups_plot_smk_sub <- groups_plot_smk + theme_light() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) 
gxmplot_bmi_sub <- gxmplot_bmi + theme_light() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
gxmplot_smk_sub <- gxmplot_smk + theme_light() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))

#gp<-grid.arrange(label_subfigure(groups_plot_bmi_sub,paste("a      ","V.E. by probes in regions dif. exp. in tissues BMI")), 
#                 label_subfigure(groups_plot_smk_sub,paste("b      ","V.E. by probes in regions dif. exp. in tissues Smoking")), 
#                 label_subfigure(gxmplot_bmi_sub,paste("c      ","Mixtures contribution BMI")), 
#                 label_subfigure(gxmplot_smk_sub,paste("d      ","Mixtures contribution Smoking")), 
#                 ncol=2)
#                 
gp <- grid.arrange(label_subfigure(gxmplot_bmi,paste("a      ","Mixtures contribution BMI")),
                   label_subfigure(gxmplot_smk,paste("b      ","Mixtures contribution Smoking")))

ggsave("VE_tissue.pdf",plot = gp, scale = 2)
