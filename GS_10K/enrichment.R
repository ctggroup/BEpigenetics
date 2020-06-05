source("./flattenAnn.R")
source("./PIPdf.R")
source("./plotGenes.R")
source("./plot_GO.R")
source("./EwasCor.R")
library(gridExtra)
library(grid)
library(data.table)

flatAnn<-flattenAnn("450k")


load("posteriorSummaryBMI.RData")
PIPdf.bmi <- PIPdf(bmi_result, flatAnn,filename = "res/PIPtable_bmi.tex" )

load("posteriorSummarySM.RData")
PIPdf.smk<- PIPdf(sm_result, flatAnn,filename = "res/PIPtable_bmi.tex" )

PIPdf<-rbind("BMI"=PIPdf.bmi,"smoking"=PIPdf.smk)
colnames(PIPdf)<-c("0<5%",'5%<50%','50%<95%','95%>')
#print(xtable(PIPdf),file="res/PIPdf.tex")

##############
##############
##############

##DANIEL continue here for the enrichment
##use the probes from genes.GO to get the squared effect sizes

titles_plots <- function( plot_type){
  titles <- list()
  titles[["Ewas"]] <- paste("Previous probe associations from EWAS catalog", plot_type)
  titles[["Effects"]] <- paste("Effects of methylation probes with 95% IP", plot_type)
  titles[["Go"]] <- "GO enrichment for the methylation probes with 95% IP"
  titles
}

    
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(Other)

genexprobe<-as.matrix(Other$UCSC_RefGene_Name)
rownames(genexprobe)<-rownames(Other)

# Enrichment plots plots
EWAS_Catalog_03_07_2019 <- read_delim("EWAS_Catalog_03-07-2019.txt","\t",escape_double = FALSE, trim_ws = TRUE)

load("posteriorSummaryBMI.RData")
titles_bmi <- titles_plots("BMI")
subfigures <- list(Effects = "a", Ewas = "c")
plots_bmi <- plotGenes(bmi_result, genexprobe = genexprobe, EwasCatalog = EWAS_Catalog_03_07_2019, titles = list(Ewas = "", Effects = "", Go = ""), subfigures)

load("posteriorSummarySM.RData")
titles_smk <- titles_plots("Smoking")
subfigures <- list(Effects = "b", Ewas = "d")
plots_smk <- plotGenes(sm_result, genexprobe = genexprobe, EwasCatalog = EWAS_Catalog_03_07_2019, titles = list(Ewas = "", Effects = "", Go = ""), subfigures)



load("enrichmentGO_bmi.RData")
enrich_bmi.df<- rbindlist(enrich_bmi)
title_GO_bmi <- "Enrichment of terms according to representation in the model BMI"
output_file <- "enrichmentxproportion_bmi.pdf"
GO.bmi <- plot_GO(enrich_bmi.df, title_GO_bmi, output_file,bmi_result$PIP,"bmi")


load("enrichmentGO_smk.RData")
enrich_sm.df<- rbindlist(enrich_smk)
title_GO_smk <- "Enrichment of terms according to representation in the model smoking"
output_file <- "enrichmentxproportion_smk.pdf"
GO.smk <- plot_GO(enrich_sm.df, title_GO_smk, output_file,sm_result$PIP,"smk")

GO.bmi <- GO.bmi + theme(legend.position = "none")
GO.smk <- GO.smk + theme(legend.postion = "none")

GO.smk <- arrangeGrob(GO.smk, top = textGrob("f", x = unit(0, "npc")
                                           , y   = unit(1, "npc"), just=c("left","top"),
                                           gp=gpar(col="black", fontsize=18, fontfamily="Helvetica")))
GO.bmi<-arrangeGrob(GO.bmi, top = textGrob("e", x = unit(0, "npc")
                                           , y   = unit(1, "npc"), just=c("left","top"),
                                           gp=gpar(col="black", fontsize=18, fontfamily="Helvetica")))



# Enrichment plots plots
EWAS_Catalog_03_07_2019 <- read_delim("EWAS_Catalog_03-07-2019.txt","\t",escape_double = FALSE, trim_ws = TRUE)
load("meanModelBMI.RData")

load("posteriorSummaryBMI.RData")
bmi_betas <- meanBetas
bmi_betas   <- tibble(CpG = names(bmi_betas), effect = as.numeric(bmi_betas))
bmi_pip <- tibble(CpG = names(bmi_result$PIP), PIP = as.numeric(bmi_result$PIP))

sum_bmi <- build_summary(bmi_betas,bmi_pip, "body mass index")

library(tidyverse)
load("meanModelSM.RData")
load("posteriorSummarySM.RData")
sm_betas <- meanBetas
sm_betas   <- tibble(CpG = names(sm_betas), effect = as.numeric(sm_betas))
sm_pip <- tibble(CpG = names(sm_result$PIP), PIP = as.numeric(sm_result$PIP))

sum_sm <-  build_summary(sm_betas,sm_pip, "smoking")


Text1 <-"circle:Median;\tthick line:5-95 percentile;\tthin line:2.5-97.5 percentile"


plot_ewas_bmi <- plot_effects_cat(sum_bmi)
plot_ewas_sm <- plot_effects_cat(sum_sm)


library(ggplotify)

#fs<-fs+theme(legend.position="none")
#g2 and g2s are the tables
#
library(cowplot)
library(ggtext)
#p <- plot_grid(trait.effects, trait.pip, align = "v", ncol = 1, rel_heights = c(0.75, 0.25))
tissues <- unique(c(unique(sum_bmi$Tissue),unique(sum_sm$Tissue)))

plots_smk$Effects = plots_smk$Effects + theme(axis.title.y = element_blank())
plots_smk$EWAS_cat = plots_smk$EWAS_cat + theme(axis.title.y = element_blank())



 
 
legend.label= "blue      Same sign\nred       Different sign\ncross    Posterior mean\ncircle    Blood\ntriangle Other tissue"
plot_ewas_bmi$effects = plot_ewas_bmi$effects   + scale_shape(
  guide = guide_legend(
    direction = "horizontal",
    title.position = "top"
  )
)  + theme(legend.position = "none") + scale_color_discrete(guide=F)  +geom_text(aes(label = legend.label, x = 12, y = 0.000002), hjust = 0,size=4) 
plot_ewas_bmi$effect <- ggplot_gtable(ggplot_build(plot_ewas_bmi$effect))
plot_ewas_bmi$effect$layout$clip[gt$layout$name == "panel"] <- "off"


#+ theme(legend.justification = c("left","bottom"), legend.position= c(0.01,0.0), legend.box = "horizontal", legend.spacing.x = unit(0.0,"cm")) 


plot_ewas_sm$effects = plot_ewas_sm$effects   + scale_shape(
  guide = guide_legend(
    direction = "horizontal",
    title.position = "top"
  )
) + theme(legend.position = "none", axis.title.y=element_blank()) + scale_color_discrete(guide=F)  
  
  
  
#+ theme(legend.justification = c("left","bottom"), legend.position= c(0.01,0.0), legend.box = "horizontal", legend.spacing.x = unit(0.0,"cm"), axis.title.y=element_blank()) + scale_color_discrete(guide=F)
#plot_ewas_sm$effects = plot_ewas_sm$effects + theme(legend.position = "none") + theme(axis.title.y = element_blank()) + scale_shape_manual(values = tissues)
  #theme(axis.title.y = element_blank(), legend.text = element_text(color = "white"),legend.title = element_text(color = "white"), legend.key = element_rect(fill = "white") )  + 
  #scale_color_discrete( guide = guide_legend(override.aes = list(color = "white")))
plot_ewas_sm$pip = plot_ewas_sm$pip + theme(axis.title.y = element_blank()) 


plot_ewas_bmi_title <- "Mean effects vs effects in EWAS catalog BMI"
plot_ewas_sm_title <- "Mean effects vs effects in EWAS catalog Smoking"

ewas_effects_bmi <- plot_grid(plot_ewas_bmi$effects, plot_ewas_bmi$pip, align = "v", ncol = 1, rel_heights = c(0.75, 0.25))
ewas_effects_sm <- plot_grid(plot_ewas_sm$effects, plot_ewas_sm$pip, align = "v", ncol = 1, rel_heights = c(0.75, 0.25))

# gp <- plot_grid(plots_bmi$Effects,plots_smk$Effects,plots_bmi$EWAS_cat+coord_flip(),plots_smk$EWAS_cat+coord_flip(),
#                 ewas_effects_bmi,ewas_effects_sm,
#                 align = 'h',
#   labels = c("a", "b", "c","d","e","f"),
#   hjust = 0,
#   nrow = 3)
# 

 gp<-grid.arrange(label_subfigure(add_sub(plots_bmi$Effects,Text1,size=10),paste("a      ",titles_bmi[["Effects"]])), 
                  label_subfigure(add_sub(plots_bmi$Effects,Text1,size=10),paste("b      ",titles_smk[["Effects"]])), 
                  label_subfigure(plots_bmi$EWAS_cat+coord_flip(),paste("c      ",titles_bmi[["Ewas"]])), 
                  label_subfigure(plots_smk$EWAS_cat+coord_flip(),paste("d      ",titles_smk[["Ewas"]])), 
                  label_subfigure(ewas_effects_bmi,paste("e      ",plot_ewas_bmi_title)),
                  label_subfigure(ewas_effects_sm,paste("f      ",plot_ewas_sm_title)),
                  ncol=2)

ggsave(plot = gp, filename = "geneticSummary.pdf",scale = 3.0)


#ggsave(filename = "geneticSummary.pdf",gp)
#dev.off()


