library(xtable) 
library(missMethyl)
library(gridExtra)
library(grid)
source("./mytopGo.R")

label_subfigure <- function(this.plot,subfigure){
  fs <- arrangeGrob(this.plot, top = textGrob(subfigure, x = unit(0, "npc")
                                      , y   = unit(1, "npc"), just=c("left","top"),
                                      gp=gpar(col="black", fontsize=18, fontfamily="Helvetica")))
  fs
}

plotTopGo <- function(go, title){
  tgo<-mytopGo(go, number=10)
  tgo$N<-NULL
  tgo$Ont<-NULL
  tgo$DE<-NULL
  tgo$FDR<-NULL
  title <- textGrob(title,gp=gpar(fontsize=12))
  padding <- unit(0.5,"line")
  print(xtable(tgo),file="GO95_bmi.tex")
  mytheme <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = 1.0)),
    colhead = list(fg_params=list(cex = 1.0)),
    rowhead = list(fg_params=list(cex = 1.0)))
  g2 <- tableGrob(tgo, rows=NULL,theme=mytheme)
  g2 <- gtable::gtable_add_rows(g2, 
                                heights = grobHeight(title)+ padding,
                                pos = 0)
  g2<-gtable::gtable_add_grob(g2,title,t=1,l =1,r=ncol(g2))
 
  g2
}

plotEwasCatalog <- function(result, EWASCatalog, title, subfigure) {
  genesAM<- merge(data.frame(probes=names(result$PIP)[result$PIP>0.95 & grepl("cg",names(result$PIP))]),Other,by.x="probes",by.y="row.names")
 
  genesEwas<-merge(genesAM,EWASCatalog,by.x="probes",by.y="CpG")
  tmp<-subset(genesAM,!(genesAM$probes %in% genesEwas$probes))
  
  theme_new <- theme_set(theme_light())
  theme_new <- theme_update(axis.text.x = element_text(angle = 45, hjust = 1,size=5,face="bold"))
  
  g3<-as_tibble(genesEwas) %>% 
    dplyr::select(Trait) %>% 
    gather(key= Trait, value = Trait) %>% 
    add_count(Trait) %>% 
    filter(n > 10) %>% 
    ggplot()+geom_bar(mapping=aes(x=Trait))+geom_point(mapping=aes(x=Trait,y=n))   +ggtitle( title )+xlab("")
 
  g3
}

plotEffect95 <- function(result, title, subfigure){
  g4<-result$gBeta$data
  probesL<-grepl("cg",g4$Parameter)
  g4$Marker<-ifelse(probesL,"Methylation","Genetic")
  f <- ggplot(g4, aes(x=median, y=reorder(Parameter, median),color=Marker)) + geom_point(size=3) +
    geom_segment(aes(x=Low, xend=High, yend=reorder(Parameter, median),color=Marker), size=1.5) +
    geom_segment(aes(x=low, xend=high, yend=reorder(Parameter, median),color=Marker), size=0.5) +
    xlab("HPD") + ylab("Parameter")+theme(axis.text=element_text(size=5),legend.position = "none")+ggtitle(title)+ylab("Methylation probe")+scale_color_discrete(name = "Marker type")+theme_light()+theme(legend.position = "none")
  
  f
}

plotGenes <- function(posteriorSummary, genexprobe, EwasCatalog, titles, subfigures ) {
  result <- posteriorSummary
  output <-list()
  go<-gometh(names(result$PIP)[result$PIP>0.95 & grepl("cg",names(result$PIP))],names(result$PIP)[grep("cg",names(result$PIP))])
  
  plot_GO_table <- tryCatch(plotTopGo(go, titles[["Go"]]),
                            error=function(e) print(e))
  if(!inherits(plot_GO_table, "error")){
    output[["GO_table"]] <- plot_GO_table
  }
  
  plot_Ewas <- tryCatch(plotEwasCatalog(result, EwasCatalog, titles[["Ewas"]], subfigures[["Ewas"]]),
                        error=function(e) print(e))
  
  if(!inherits(plot_Ewas, "error")){
    output[["EWAS_cat"]] <- plot_Ewas
  }
  
  plot_effects <- tryCatch(plotEffect95(result, titles[["Effects"]], subfigures[["Effects"]]),
                           error= function(e) print(e))
  if(!inherits(plot_effects, "error")){
    output[["Effects"]] <- plot_effects
  }
  
 output
}
  
