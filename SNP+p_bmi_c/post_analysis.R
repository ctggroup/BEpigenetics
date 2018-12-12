library(coda)
library(ggmcmc)
library(parallel)
library(bigmemory)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(matrixStats)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
#function to extract parametres from chains, we turn a big.matrix object into normal matrices objects.
extract_parameter <- function(chain,parameter){
  as.mcmc(chain[,grep(parameter,colnames(chain))])
}
changeNamesSigma<-function(x){
  colnames(x)[grep("sigma",colnames(x))]<-c("sigma[epsilon]^2","sigma[G]^2","sigma[phi]^2")
  x
}

process_posterior<- function(chains,probe_names,genexprobe,subsample){
#args <- commandArgs()
# Calculate the number of cores
no_cores <- 4
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
save(list=c("C1","C2","C3","C4"),file="posterior.RData")


#again the bigmatrix objects seem to run into heap issues if I try a lapply over them, thus I will first convert the individual chain parameters to objects and then operate using normal matrix objects
print("extracting sigma parameter")
chainsSigma<-as.mcmc.list(lapply(list(C1[subsample,],C2[subsample,],C3[subsample,],C4[subsample,]),function(x){
  changeNamesSigma(extract_parameter(x,"sigma"))}))
print("extracting beta parameters")
chainsBeta<- as.mcmc.list(lapply(list(C1[subsample,],C2[subsample,],C3[subsample,],C4[subsample,]),function(x){
  tmp<-extract_parameter(x,"beta")
  colnames(tmp)<-probe_names
  tmp
}))
print("extracting mu parameter")
chainsMu <- as.mcmc.list(lapply(list(C1[subsample,],C2[subsample,],C3[subsample,],C4[subsample,]),function(x){extract_parameter(x,"mu")}))
print("extracting components parameters")
chainsComp <- as.mcmc.list(lapply(list(C1[subsample,],C2[subsample,],C3[subsample,],C4[subsample,]),function(x){extract_parameter(x,"comp")}))
#We free up some memory
print("freeing up space")
rm(C1,C2,C3,C4)
gc()

#now we print the traceplots, and the autocorrelation plots
print("printing the auto corrleation plot")


sSigma<-ggs(chainsSigma)

pdf("Autocorrelation.pdf")
g<-ggs_autocorrelation(sSigma,greek=T)
print(g)
dev.off()
sSigma<-ggs(chainsSigma)
print("printing the traceplot")
pdf("Traceplot.pdf")
g<-ggs_traceplot(sSigma,greek=T)
print(g)
dev.off()

print("printing the running means plot")
pdf("Running.pdf")
g<-ggs_running(ggs(chainsSigma),greek=T)
print(g)
dev.off()

print("printing the R_hat plot")
pdf("Rhat.pdf")
g<-ggs_Rhat(sSigma,greek=T)+xlab("R_hat")
print(g)
dev.off()

print("printing the Geweke plot")
pdf("Geweke.pdf")
g<-ggs_geweke(sSigma,greek=T)
print(g)
dev.off()

#Here we should do some thinning if necessary


#here we compute the chains inclussion probability
print("calculating markers in model per sample")
postMarkers<-apply(X = do.call( rbind,chainsComp), MARGIN=2, function(x){as.numeric(x!=0)} )
colnames(postMarkers)<-probe_names

print("calculating proportion of genes mapped to probes per sample")
numberOfGenes<-length(unique(genexprobe))
genesPerSample<-apply(postMarkers,MARGIN=1,FUN=function(x){length(unique(genexprobe[rownames(genexprobe) %in% colnames(postMarkers)[which(abs(x)>0)],1]))/numberOfGenes})



print("Calculating posterior inclusion probability")
postIncl<-apply(X = postMarkers, MARGIN=2, function(x){sum(x)/length(x)} )
names(postIncl)<-probe_names
print("Calculating the variance explained by the probes with 95% IP")
VE95<- apply(X=do.call(rbind,chainsBeta[,which(postIncl>0.95)]),MARGIN=1,FUN=function(x)var(x)*length(which(postIncl>0.95)))/rowSums(do.call(rbind,chainsSigma))

print("Creating grob for the probes with 95% IP")
sBeta<-ggs(chainsBeta[,which(postIncl>0.95),drop=FALSE])
gBeta<- ggs_caterpillar(sBeta)

print("calculating markers in model summaries")
markersInModel<-apply(X=do.call(rbind,chainsComp),MARGIN=1,FUN = function(x){sum(x!=0)})
probesInModel<-apply(X=do.call(rbind,chainsComp),MARGIN=1,FUN = function(x){sum(x[grep("cg",probe_names)]!=0)})
snpsInModel<-apply(X=do.call(rbind,chainsComp),MARGIN=1,FUN = function(x){sum(x[grep("rs",probe_names)]!=0)})

print("calculating snps in mixture per sample")
SNPper1<-apply(X=do.call(rbind,chainsComp),MARGIN = 1,FUN=function(x){sum(as.numeric(x[grep("rs",probe_names)]==1))})
SNPper2<-apply(X=do.call(rbind,chainsComp),MARGIN = 1,FUN=function(x){sum(as.numeric(x[grep("rs",probe_names)]==2))})
SNPper3<-apply(X=do.call(rbind,chainsComp),MARGIN = 1,FUN=function(x){sum(as.numeric(x[grep("rs",probe_names)]==3))})
print("calculating oribes in mixture per sample")
Pper1<-apply(X=do.call(rbind,chainsComp),MARGIN = 1,FUN=function(x){sum(as.numeric(x[grep("cg",probe_names)]==1))})
Pper2<-apply(X=do.call(rbind,chainsComp),MARGIN = 1,FUN=function(x){sum(as.numeric(x[grep("cg",probe_names)]==2))})
Pper3<-apply(X=do.call(rbind,chainsComp),MARGIN = 1,FUN=function(x){sum(as.numeric(x[grep("cg",probe_names)]==3))})


#sBetas<-ggs(chainsBeta)
#sComps<-ggs(chainsComp)

mask1<-apply(X=do.call(rbind,chainsComp),MARGIN = 1,FUN=function(x){as.numeric(x==1)})
mask2<-apply(X=do.call(rbind,chainsComp),MARGIN = 1,FUN=function(x){as.numeric(x==2)})
mask3<-apply(X=do.call(rbind,chainsComp),MARGIN = 1,FUN=function(x){as.numeric(x==3)})

print("calculating explained variance per mixture probes")
varExplained1S=colSums(apply(X=t(mask1)*do.call(rbind,chainsBeta),MARGIN = 1,FUN=function(x){(x[grep("rs",probe_names)]^2)})/(do.call(rbind,chainsSigma)[,"sigma[G]^2"]*0.0001))/SNPper1
varExplained1S=apply(X=t(mask1)*do.call(rbind,chainsBeta),MARGIN=1,FUN=function(x)var(x[grepl("rs",probe_names) & x!=0]))*SNPper1/rowSums(do.call(rbind,chainsSigma))
varExplained2S=apply(X=t(mask2)*do.call(rbind,chainsBeta),MARGIN=1,FUN=function(x)var(x[grepl("rs",probe_names) & x!=0]))*SNPper2/rowSums(do.call(rbind,chainsSigma))
varExplained3S=apply(X=t(mask3)*do.call(rbind,chainsBeta),MARGIN=1,FUN=function(x)var(x[grepl("rs",probe_names) & x!=0]))*SNPper3/rowSums(do.call(rbind,chainsSigma))


#varExplained2S=colSums(apply(X=t(mask2)*do.call(rbind,chainsBeta),MARGIN = 1,FUN=function(x){(x[grep("rs",probe_names)]^2)})/(do.call(rbind,chainsSigma)[,"sigma[G]^2"]*0.001))/SNPper2
#varExplained3S=colSums(apply(X=t(mask3)*do.call(rbind,chainsBeta),MARGIN = 1,FUN=function(x){(x[grep("rs",probe_names)]^2)})/(do.call(rbind,chainsSigma)[,"sigma[G]^2"]*0.01))/SNPper3



print("calculating explained variance per mixture probes")
#varExplained1=colSums(apply(X=t(mask1)*do.call(rbind,chainsBeta),MARGIN = 1,FUN=function(x){(x[grep("cg",probe_names)]^2)})/(do.call(rbind,chainsSigma)[,"sigma[phi]^2"]*0.01))/Pper1
#varExplained2=colSums(apply(X=t(mask2)*do.call(rbind,chainsBeta),MARGIN = 1,FUN=function(x){(x[grep("cg",probe_names)]^2)})/(do.call(rbind,chainsSigma)[,"sigma[phi]^2"]*0.1))/Pper2
#varExplained3=colSums(apply(X=t(mask3)*do.call(rbind,chainsBeta),MARGIN = 1,FUN=function(x){(x[grep("cg",probe_names)]^2)})/(do.call(rbind,chainsSigma)[,"sigma[phi]^2"]*1))/Pper3
varExplained1=apply(X=t(mask1)*do.call(rbind,chainsBeta),MARGIN=1,FUN=function(x)var(x[grepl("cg",probe_names) & x!=0]))*Pper1/rowSums(do.call(rbind,chainsSigma))
varExplained2=apply(X=t(mask2)*do.call(rbind,chainsBeta),MARGIN=1,FUN=function(x)var(x[grepl("cg",probe_names) & x!=0]))*Pper2/rowSums(do.call(rbind,chainsSigma))
varExplained3=apply(X=t(mask3)*do.call(rbind,chainsBeta),MARGIN=1,FUN=function(x)var(x[grepl("cg",probe_names) & x!=0]))*Pper3/rowSums(do.call(rbind,chainsSigma))


varExplained1[is.na(varExplained1)] <- 0
varExplained2[is.na(varExplained2)] <- 0
varExplained3[is.na(varExplained3)] <- 0
varExplained1S[is.na(varExplained1S)] <- 0
varExplained2S[is.na(varExplained2S)] <- 0
varExplained3S[is.na(varExplained3S)] <- 0



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
     VE95=VE95,
     gBeta=gBeta,
     VEepi=do.call(rbind,chainsSigma)[,"sigma[phi]^2"]/rowSums(do.call(rbind,chainsSigma)),
     VEgene=do.call(rbind,chainsSigma)[,"sigma[G]^2"]/rowSums(do.call(rbind,chainsSigma))
     )

}

files<-c("C1bmi_c.csv","C2bmi_c.csv","C3bmi_c.csv","C4bmi_c.csv")
load("bmi_col_names.RData")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(Other)
genexprobe<-as.matrix(Other$UCSC_RefGene_Name)
rownames(genexprobe)<-rownames(Other)

result<-process_posterior(files,bmi_names,genexprobe,750:1000)
save(list='result',file="posteriorSummary.RData")
