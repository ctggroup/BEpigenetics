#First we load the thinned chains
#
library(data.table)
library(tidyverse)
library(coda)
library(ggmcmc)

#function to extract parametres from chains, we turn a big.matrix object into normal matrices objects.
extract_parameter <- function(chain,parameter){
  as.mcmc(chain[,grep(parameter,colnames(chain))])
}
changeNamesSigma<-function(x){
  if(length(grep("sigma",colnames(x)))==3)
    colnames(x)[grep("sigma",colnames(x))]<-c("sigma[epsilon]^2","sigma[phi]^2","sigma[G]^2")
  else
    colnames(x)[grep("sigma",colnames(x))]<-c("sigma[epsilon]^2","sigma[phi]^2","sigma[G]^2","sigma[C]^2")
  x
}

thinned_process_posterior <- function(chains, probe_names, genexprobe){
  chains <- as_tibble(chains)
  print("extracting sigma parameter")
  chainsSigma<- changeNamesSigma(extract_parameter(chains,"sigma"))

  print("extracting beta parameters")
  chainsBeta<-extract_parameter(chains,"beta")
  colnames(chainsBeta)<-probe_names
    
  print("extracting mu parameter")
  chainsMu <- extract_parameter(chains,"mu")
  
  print("extracting components parameters")
  chainsComp <- extract_parameter(chains,"comp")
  colnames(chainsComp) <- probe_names
  
  print("calculating markers in model per sample")
  postMarkers<-apply(X = chainsComp, MARGIN=2, function(x){as.numeric(x!=0)} )
  colnames(postMarkers)<-probe_names
  
  
  print("Calculating posterior inclusion probability")
  postIncl<-apply(X = postMarkers, MARGIN=2, function(x){sum(x)/length(x)} )
  names(postIncl)<-probe_names
  print("Calculating the variance explained by the probes with 95% IP")
  VE95<- apply(X=chainsBeta[,which(postIncl>0.95)],MARGIN=1,FUN=function(x)var(x)*length(which(postIncl>0.95)))/rowSums(chainsSigma)
  
  print("Creating grob for the probes with 95% IP")
  sBeta<-ggs(chainsBeta[,which(postIncl>0.95),drop=FALSE])
  gBeta<- ggs_caterpillar(sBeta)
  
  
  print("calculating markers in model summaries")
  markersInModel<-apply(X=chainsComp,MARGIN=1,FUN = function(x){sum(x!=0)})
  probesInModel<-apply(X=chainsComp,MARGIN=1,FUN = function(x){sum(x[grep("cg",probe_names)]!=0)})
  snpsInModel<-apply(X=chainsComp,MARGIN=1,FUN = function(x){sum(x[grep("rs",probe_names)]!=0)})

  
  print("calculating snps in mixture per sample")
  SNPper1<-apply(X=chainsComp,MARGIN = 1,FUN=function(x){sum(as.numeric(x[grep("rs",probe_names)]==1))})
  SNPper2<-apply(X=chainsComp,MARGIN = 1,FUN=function(x){sum(as.numeric(x[grep("rs",probe_names)]==2))})
  SNPper3<-apply(X=chainsComp,MARGIN = 1,FUN=function(x){sum(as.numeric(x[grep("rs",probe_names)]==3))})
  print("calculating probes in mixture per sample")
  Pper1<-apply(X=chainsComp,MARGIN = 1,FUN=function(x){sum(as.numeric(x[grep("cg",probe_names)]==1))})
  Pper2<-apply(X=chainsComp,MARGIN = 1,FUN=function(x){sum(as.numeric(x[grep("cg",probe_names)]==2))})
  Pper3<-apply(X=chainsComp,MARGIN = 1,FUN=function(x){sum(as.numeric(x[grep("cg",probe_names)]==3))})  
  
  mask1<-apply(X=chainsComp,MARGIN = 1,FUN=function(x){as.numeric(x==1)})
  mask2<-apply(X=chainsComp,MARGIN = 1,FUN=function(x){as.numeric(x==2)})
  mask3<-apply(X=chainsComp,MARGIN = 1,FUN=function(x){as.numeric(x==3)})
  
  
  #variance explained computed as variance(Betas_in_Mixture)/total(variance)
  
  print("calculating explained variance per mixture SNP")
  
  varExplained1S=apply(X=t(mask1)*chainsBeta,MARGIN=1,FUN=function(x)var(x[grepl("rs",probe_names) & x!=0]))*SNPper1/rowSums(chainsSigma)
  varExplained2S=apply(X=t(mask2)*chainsBeta,MARGIN=1,FUN=function(x)var(x[grepl("rs",probe_names) & x!=0]))*SNPper2/rowSums(chainsSigma)
  varExplained3S=apply(X=t(mask3)*chainsBeta,MARGIN=1,FUN=function(x)var(x[grepl("rs",probe_names) & x!=0]))*SNPper3/rowSums(chainsSigma)
  
  
  
  print("calculating explained variance per mixture probes")
  #variance explained computed as variance(Betas_in_Mixture)/total(variance)
  
  varExplained1=apply(X=t(mask1)*chainsBeta,MARGIN=1,FUN=function(x)var(x[grepl("cg",probe_names) & x!=0]))*Pper1/rowSums(chainsSigma)
  varExplained2=apply(X=t(mask2)*chainsBeta,MARGIN=1,FUN=function(x)var(x[grepl("cg",probe_names) & x!=0]))*Pper2/rowSums(chainsSigma)
  varExplained3=apply(X=t(mask3)*chainsBeta,MARGIN=1,FUN=function(x)var(x[grepl("cg",probe_names) & x!=0]))*Pper3/rowSums(chainsSigma)
  
  
  
  
  varExplained1[is.na(varExplained1)] <- 0
  varExplained2[is.na(varExplained2)] <- 0
  varExplained3[is.na(varExplained3)] <- 0
  varExplained1S[is.na(varExplained1S)] <- 0
  varExplained2S[is.na(varExplained2S)] <- 0
  varExplained3S[is.na(varExplained3S)] <- 0
  
  print("calculating proportion of genes mapped to probes per sample")
  numberOfGenes<-length(unique(genexprobe))
  genesPerSample<-apply(postMarkers,MARGIN=1,FUN=function(x){length(unique(genexprobe[rownames(genexprobe) %in% colnames(postMarkers)[which(abs(x)>0)],1]))/numberOfGenes})
  list(PIP=postIncl, #posterior inclusion probabilities
       genesPerSample=genesPerSample,
       perSampleMarkers=markersInModel, #number of markers in model per sample
       perSampleProbes=probesInModel, #number of probes in model per sample
       perSampleSNP=snpsInModel, #number of SNP in model per sample
       SNPper1=SNPper1, #number of snps in mixture one
       SNPper2=SNPper2, #number of snps in mixture 2
       SNPper3=SNPper3, #number of snps in mixture3
       Pper1=Pper1, #number of probes in mixture 1
       Pper2=Pper2, #number of probes in mixture2
       Pper3=Pper3, #number of probes in mixture 3
       varP1=varExplained1, #variance explained by probes in mixture 1
       varP2=varExplained2, #variance explained by probes in mixture 2
       varP3=varExplained3, #variance explained by probes in mixture 3
       varSNP1=varExplained1S, #variance explained by SNP in mixture 1
       varSNP2=varExplained2S, #variance explained by SNP in mixture 2
       varSNP3=varExplained3S, #variance explained vy snp in mixture 3
       VE95=VE95,#variance explained by probes wiht 95% incl. prob
       gBeta=gBeta,#ggplot with the effects
       VEepi=chainsSigma[,"sigma[phi]^2"]/rowSums(chainsSigma),#VE by probes
       VEgene=chainsSigma[,"sigma[G]^2"]/rowSums(chainsSigma) #VE SNP
  )
}