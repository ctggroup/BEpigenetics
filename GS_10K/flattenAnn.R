library(org.Hs.eg.db)
library(plyr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(missMethyl)
library(limma)
flattenAnn <- function(array.type)
  # flatten 450k or EPIC array annotation
  # Belinda Phipson
  # 10 February 2016
  # Updated 7 July 2016
{
  #if(array.type=="450K")    
  anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  #else
  # anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  
  # get rid of the non-CpG sites
  strlen<-str_length(rownames(anno))
  ann.keep<-anno[strlen==10,]
  
  # get rid of CpGs that are not annotated
  missing<-ann.keep$UCSC_RefGene_Name==""
  ann.keep<-ann.keep[!missing,]
  
  # get individual gene names for each CpG
  geneslist<-strsplit(ann.keep$UCSC_RefGene_Name,split=";")
  names(geneslist)<-rownames(ann.keep)
  
  grouplist<-strsplit(ann.keep$UCSC_RefGene_Group,split=";")
  names(grouplist)<-rownames(ann.keep)
  
  flat<-data.frame(symbol=unlist(geneslist),group=unlist(grouplist))
  flat$symbol<-as.character(flat$symbol)
  flat$group <- as.character(flat$group)
  
  flat$cpg<- substr(rownames(flat),1,10)
  
  flat$alias <- alias2SymbolTable(flat$symbol)
  
  eg <- toTable(org.Hs.egSYMBOL2EG)
  m <- match(flat$alias,eg$symbol)
  flat$entrezid <- eg$gene_id[m]
  flat <- flat[!is.na(flat$entrezid),]
  
  # keep unique cpg by gene name annotation
  id<-paste(flat$cpg,flat$entrezid,sep=".")
  d <- duplicated(id)
  flat.u <- flat[!d,]
  flat.u
}