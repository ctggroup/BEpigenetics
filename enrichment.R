library(org.Hs.eg.db)
library(plyr)
flattenAnn <- function(array.type)
  # flatten 450k or EPIC array annotation
  # Belinda Phipson
  # 10 February 2016
  # Updated 7 July 2016
{
  if(array.type=="450K")    
    anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  else
    anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  
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

tmp<-flattenAnn("450k")


load("../SNP+p_bmi_c/posteriorSummary.RData")
result.bmi<-result
genes.bmi<-merge(data.frame(cpg=names(result$PIP),PIP=as.numeric(result$PIP)),tmp,by.x="cpg", by.y="cpg")
genes.bmi<-genes.bmi[order(genes.bmi$PIP,decreasing = T),]
genes.bmi<-genes.bmi[genes.bmi$PIP>0,]
#print(xtable(genes.bmi),file="PIPtable_bmi.tex")

PIPdf.bmi<-data.frame("0<5%"=length(which(result$PIP>0 & result$PIP<0.05)),'5%<50%'=length(which(result$PIP>0.05 & result$PIP<0.5)),'50%<95%'=length(which(result$PIP>0.5 & result$PIP<0.95)),'95%>'=length(which(result$PIP>0.95)))


load("../SNP+p_smk/posteriorSummary.RData")
result.smk<-result
genes.smk<-merge(data.frame(cpg=names(result$PIP),PIP=as.numeric(result$PIP)),tmp,by.x="cpg", by.y="cpg")
genes.smk<-genes.smk[order(genes.smk$PIP,decreasing = T),]
genes.smk<-genes.smk[genes.smk$PIP>0,]
#print(xtable(genes.smk),file="PIPtable_smk.tex")

PIPdf.smk<-data.frame("0<5%"=length(which(result$PIP>0 & result$PIP<0.05)),'5%<50%'=length(which(result$PIP>0.05 & result$PIP<0.5)),'50%<95%'=length(which(result$PIP>0.5 & result$PIP<0.95)),'95%>'=length(which(result$PIP>0.95)))

PIPdf<-rbind("BMI"=PIPdf.bmi,"smoking"=PIPdf.smk)
colnames(PIPdf)<-c("0<5%",'5%<50%','50%<95%','95%>')
print(xtable(PIPdf),file="PIPdf.tex")

egGO2ALLEGS<-getFromNamespace("org.Hs.egGO2ALLEGS","org.Hs.eg.db")
GeneID.PathID <- AnnotationDbi::toTable(egGO2ALLEGS)[,c("gene_id","go_id","Ontology")]
genes.GO.smk<-merge(genes.smk,GeneID.PathID,by.x="entrezid",by.y="gene_id")
genes.GO.bmi<-merge(genes.smk,GeneID.PathID,by.x="entrezid",by.y="gene_id")


head(genes.GO[order(genes.GO$go_id),])
library(plyr)
count(genes.GO.bmi,"go_id")
load("~/Documents/Methylation/SNP+p_bmi_c/mean_model_bmi_c.RData")
betas.bmi<-meanBetas
betas.bmi<-data.frame("cpg"=names(betas.bmi),"effect"=as.numeric(betas.bmi))
load("~/Documents/Methylation/SNP+p_smk/mean_model_smk.RData")
betas.smk<-meanBetas
betas.smk<-data.frame("cpg"=names(betas.smk),"effect"=as.numeric(betas.smk))

betas.bmi<-merge(genes.GO.bmi,betas.bmi)
betas.smk<-merge(genes.GO.smk,betas.smk)

VE.betas.bmi<-aggregate(effect ~ go_id,data=betas.bmi,function(x){sum(x^2)})/mean(result.bmi$VEepi)
VE.betas.smk<-aggregate(effect ~ go_id,data=betas.smk,function(x){sum(x^2)})/mean(result.bmi$VEepi)
VE.betas.bmi<-VE.betas.bmi[!is.na(VE.betas.bmi$effect),]
VE.betas.smk<-VE.betas.smk[!is.na(VE.betas.smk$effect),]

VE.bmi<-count(genes.GO.bmi,"go_id")
VE.bmi$VE<-VE.bmi$freq/length(meanBetas)

VE.plot.bmi<-merge(VE.bmi,VE.betas.bmi)
VE.plot.bmi<-VE.plot.bmi[VE.plot.bmi$effect!=0,]
VE.plot.bmi$effect<-VE.plot.bmi$effect/VE.plot.bmi$VE
TERM <- suppressMessages(AnnotationDbi::select(GO.db::GO.db,keys=VE.plot.bmi$go_id,columns="TERM"))
VE.plot.bmi$TERM<-TERM$TERM

p<-ggplot(data=VE.plot.bmi, mapping=aes(x=VE.plot.bmi$VE,y=VE.plot.bmi$effect)) +geom_point() +geom_line(data= VE.plot.bmi,aes(x=VE.plot.bmi$VE,y=VE.plot.bmi$VE))
plot(p)
head(VE.plot.bmi[order(VE.plot.bmi$effect,decreasing=T),])


VE.smk<-count(genes.GO.smk,"go_id")
VE.smk$VE<-VE.smk$freq/length(meanBetas)


VE.plot.smk<-merge(VE.smk,VE.betas.smk)
VE.plot.smk<-VE.plot.smk[VE.plot.smk$effect!=0,]
VE.plot.smk$effect<-VE.plot.smk$effect/VE.plot.smk$VE
TERM <- suppressMessages(AnnotationDbi::select(GO.db::GO.db,keys=VE.plot.smk$go_id,columns="TERM"))
VE.plot.smk$TERM<-TERM$TERM

p<-ggplot(data=VE.plot.smk, mapping=aes(x=VE.plot.smk$VE,y=VE.plot.smk$effect)) +geom_point() +geom_line(data= VE.plot.smk,aes(x=VE.plot.smk$VE,y=VE.plot.smk$VE))
plot(p)
head(VE.plot.smk[order(VE.plot.smk$effect,decreasing=T),],10)

##Sample wise enrichment
##



