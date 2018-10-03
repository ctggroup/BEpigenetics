library(plyr)
library(org.Hs.eg.db)
library(coda)
library(ggmcmc)
library(bigmemory)
library(missMethyl)
library(limma)
#source("https://bioconductor.org/biocLite.R")
#biocLite("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(stringr)
flattenAnn <- function(array.type)
  # flatten 450k or EPIC array annotation
  # Belinda Phipson
  # 10 February 2016
  # Updated 7 July 2016
{
 # if(array.type=="450K")    
    anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  #else
#    anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  
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
extract_parameter <- function(chain,parameter){
  as.mcmc(chain[,grep(parameter,colnames(chain))])
}
changeNamesSigma<-function(x){
  colnames(x)[grep("sigma",colnames(x))]<-c("sigma[epsilon]^2","sigma[G]^2","sigma[phi]^2")
  x
}
enrichment<-function(chains,GeneID.PathID,probe_names,flatAnnon=flattenAnn("450k"),subsample){
  #for some reason, if I do this using lapply I run into heap problems, probably the list of bigmemory objects is not well allocated
  file<-chains[1]
  print(paste("Reading file: ",file))
  C1<-read.big.matrix(file,header=T,type="double")
  file<-chains[2]
  print(paste("Reading file: ",file))
  C2<-read.big.matrix(file,header=T,type="double")
  file<-chains[3]
  print(paste("Reading file: ",file))
  C3<-read.big.matrix(file,header=T,type="double")
  file<-chains[4]
  print(paste("Reading file: ",file))
  C4<-read.big.matrix(file,header=T,type="double")
  
  
  print("subsampling and concatenating chains")
  
  
  chainsComp <- as.mcmc.list(lapply(list(C1[subsample,],C2[subsample,],C3[subsample,],C4[subsample,]),function(x){tmp<-extract_parameter(x,"comp")
  colnames(tmp)<-probe_names
  tmp
  }))
  
  print("extracting beta parameters")
  chainsBeta<- as.mcmc.list(lapply(list(C1[subsample,],C2[subsample,],C3[subsample,],C4[subsample,]),function(x){
    tmp<-extract_parameter(x,"beta")
    colnames(tmp)<-probe_names
    tmp[,grep("cg",probe_names)]
  }))
  print("extracting sigma parameter")
  chainsSigma<-as.mcmc.list(lapply(list(C1[subsample,],C2[subsample,],C3[subsample,],C4[subsample,]),function(x){
    changeNamesSigma(extract_parameter(x,"sigma"))}))
  tmp<-do.call(rbind,chainsSigma)[,"sigma[phi]^2"]/rowSums(do.call(rbind,chainsSigma)) 
  #M<-length(probe_names[grep("cg",probe_names)])
  #print(M)
#  enrich<-apply(X = diag(1/sqrt(tmp)) %*% do.call(rbind,chainsBeta)    ,MARGIN = 1,FUN=function(betas){
#            
    enrich<-apply(X = do.call(rbind,chainsBeta)    ,MARGIN = 1,FUN=function(betas){
          M<-length(betas[betas!=0])
           annonBeta<-merge(data.frame(cpg=names(betas)[betas!=0],effect=as.numeric(betas[betas!=0])),flatAnnon,by.x="cpg", by.y="cpg")
           annonTerm<-merge(annonBeta,GeneID.PathID,by.x="entrezid",by.y="gene_id")
           Mb<-sum(betas^2)
           VE.betas<-aggregate(effect ~ go_id,data=annonTerm,function(x){sum(x^2)})
           VE.betas<-VE.betas[!is.na(VE.betas),]
           VE.betas$P.VE<-VE.betas$effect/Mb
           VE.betas$effect<-NULL
           print("number of non-zero effects")
           print(M)
           annonTerm<-dplyr::distinct(annonTerm,go_id,cpg)
           P.probes<-plyr::count(annonTerm,"go_id")
         
           print("max number of probes per term")
          
           P.probes$P.probes<-P.probes$freq/M
           print(max(P.probes$P.probes))
           
           VE.plot<-merge(P.probes,VE.betas)
           VE.plot<-VE.plot[VE.plot$P.probes!=0 & VE.plot$P.VE!=0,]
           VE.plot$enrichment<-VE.plot$P.VE/VE.plot$P.probes
           VE.plot<-VE.plot[VE.plot$freq>2,]
            annonTerm
         })
   enrich
   #save(list="enrich",file="enrichmentGO.RData")
}
load("bmi_col_names.RData")
egGO2ALLEGS<-getFromNamespace("org.Hs.egGO2ALLEGS","org.Hs.eg.db")
GeneID.PathID <- AnnotationDbi::toTable(egGO2ALLEGS)[,c("gene_id","go_id","Ontology")]
files=c("C1bmi.csv","C2bmi.csv","C3bmi.csv","C4bmi.csv")
enrich<-enrichment(files,GeneID.PathID,probe_names=bmi_names,flatAnnon=flattenAnn("450k"),subsample=1:2) 
tmp<-enrich[[8]]

 tmp.count<-plyr::count(tmp,"go_id")
 head(tmp.count[order(tmp.count$freq,decreasing = T),])


