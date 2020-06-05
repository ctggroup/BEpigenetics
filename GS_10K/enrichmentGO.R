library(plyr)
library(org.Hs.eg.db)
library(coda)
library(ggmcmc)
library(missMethyl)
library(limma)
#source("https://bioconductor.org/biocLite.R")
#biocLite("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(stringr)
source("./flattenAnn.R")

extract_parameter <- function(chain,parameter){
  as.data.frame(chain[,grep(parameter,colnames(chain)),with=F])
}
changeNamesSigma<-function(x){
  colnames(x)[grep("sigma",colnames(x))]<-c("sigma[epsilon]^2","sigma[phi]^2","sigma[G]^2")
  x
}
enrichment<-function(C1,C2,C3,C4,GeneID.PathID,probe_names,flatAnnon=flattenAnn("450k"),subsample){
  #for some reason, if I do this using lapply I run into heap problems, probably the list of bigmemory objects is not well allocated
  
  print("subsampling and concatenating chains")
  
  
  chainsComp <- lapply(list(C1,C2,C3,C4),function(x){tmp<-extract_parameter(x,"comp")
  colnames(tmp)<-probe_names
  tmp
  })
  
  print("extracting beta parameters")
  chainsBeta<- lapply(list(C1,C2,C3,C4),function(x){
    tmp<-extract_parameter(x,"beta")
    colnames(tmp)<-probe_names
    tmp[,grep("cg",probe_names)]
  })
  print("extracting sigma parameter")
  chainsSigma<-lapply(list(C1,C2,C3,C4),function(x){
    changeNamesSigma(extract_parameter(x,"sigma"))})
  tmp<-rbindlist(chainsSigma)[,"sigma[phi]^2"]/rowSums(do.call(rbind,chainsSigma)) 
  #M<-length(probe_names[grep("cg",probe_names)])
  #print(M)
  #  enrich<-apply(X = diag(1/sqrt(tmp)) %*% do.call(rbind,chainsBeta)    ,MARGIN = 1,FUN=function(betas){
  #            
  enrich<-apply(X = rbindlist(chainsBeta)    ,MARGIN = 1,FUN=function(betas){
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
    
    VE.plot<-merge(merge(P.probes,VE.betas),annonTerm)
    VE.plot<-VE.plot[VE.plot$P.probes!=0 & VE.plot$P.VE!=0,]
    VE.plot$enrichment<-VE.plot$P.VE/VE.plot$P.probes
    VE.plot<-VE.plot[VE.plot$freq>2,]
    VE.plot
  })
  enrich
  #
}
