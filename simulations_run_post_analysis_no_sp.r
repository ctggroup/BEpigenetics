
lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    p
    }

#function to extract parametres from chains, we turn a big.matrix object into normal matrices objects.
extract_parameter <- function(chain,parameter){
  as.mcmc(chain[,grep(parameter,colnames(chain))])
}
changeNamesSigma<-function(x){
  colnames(x)[grep("sigma",colnames(x))]<-c("sigma[epsilon]^2","sigma[G]^2","sigma[phi]^2")
  x
}

process_posterior<- function(chains,B,y,O,M,R,osca.pve,osca.B){
#args <- commandArgs()
# Calculate the number of cores
#for some reason, if I do this using lapply I run into heap problems, probably the list of bigmemory objects is not well allocated

corR<-sqrt(sum(cor(y,R)^2))
covR<-sqrt(sum(cov(y,R)^2))
file <- paste(chains[1],"/C1.csv",sep="")
print(paste("Reading file: ",file))
C1<-read.big.matrix(file,header=T,type="double")
save(list=c("C1"),file=paste(chains[1],"posterior.RData",sep=""))


#again the bigmatrix objects seem to run into heap issues if I try a lapply over them, thus I will first convert the individual chain parameters to objects and then operate using normal matrix objects
print("extracting sigma parameter")
chainsSigma<-as.mcmc(lapply(list(C1),function(x){
  changeNamesSigma(extract_parameter(x,"sigma"))}))
print("extracting beta parameters")
chainsBeta<- as.mcmc(lapply(list(C1),function(x){
  tmp<-extract_parameter(x,"beta")
  tmp
}))
print("extracting mu parameter")
chainsMu <- as.mcmc(lapply(list(C1),function(x){extract_parameter(x,"mu")}))
print("extracting components parameters")
chainsComp <- as.mcmc(lapply(list(C1),function(x){extract_parameter(x,"comp")}))
#We free up some memory
print("freeing up space")
rm(C1)
gc()

#now we print the traceplots, and the autocorrelation plots
print("printing the auto corrleation plot")


sSigma<-ggs(chainsSigma)

pdf(paste(chains[1],"/Autocorrelation.pdf",sep=""))
g<-ggs_autocorrelation(sSigma,greek=T)
print(g)
dev.off()
sSigma<-ggs(chainsSigma)
print(paste(chains[1],"printing the traceplot",sep=""))
pdf(paste(chains[1],"/Traceplot.pdf",sep=""))
g<-ggs_traceplot(sSigma,greek=T)
print(g)
dev.off()

print("printing the running means plot")
pdf(paste(chains[1],"/Running.pdf",sep=""))
g<-ggs_running(ggs(chainsSigma),greek=T)
print(g)
dev.off()


print("printing the Geweke plot")
pdf(paste(chains[1],"/Geweke.pdf",sep=""))
g<-ggs_geweke(sSigma,greek=T)
print(g)
dev.off()

#Here we should do some thinning if necessary


#here we compute the chains inclussion probability
print("calculating markers in model per sample")
postMarkers<-apply(X = do.call( rbind,chainsComp), MARGIN=2, function(x){as.numeric(x!=0)} )

print("Calculating posterior inclusion probability")
postIncl<-apply(X = postMarkers, MARGIN=2, function(x){sum(x)/length(x)} )
print("Calculating the variance explained by the probes with 95% IP")
VE95<- apply(X=do.call(rbind,chainsBeta[,which(postIncl>0.95)]),MARGIN=1,FUN=function(x)var(x)*length(which(postIncl>0.95)))/rowSums(do.call(rbind,chainsSigma))


print("calculating markers in model summaries")
markersInModel<-apply(X=do.call(rbind,chainsComp),MARGIN=1,FUN = function(x){sum(x!=0)})

meanBetas<-colMeans(do.call(rbind,chainsBeta))
meanMu<-mean(do.call(rbind,chainsMu))
meanSigmaG<-colMeans(do.call(rbind,chainsSigma))[2]
meanSigmaE<-colMeans(do.call(rbind,chainsSigma))[1]
print(meanSigmaG)
print("calculating mse")
mse_mean_B=sqrt(mean((B-meanBetas)^2))
mse_mean_B
print("calculating correlation")
cor_B=cor(B,meanBetas)
lm_summary=summary(lm(data=data.frame(y=scale(y),x = meanMu-scale(O)%*%meanBetas),formula=y~x))

lm_truth=summary(lm(data=data.frame(y=y,x=scale(M)%*%B),formula=y~x))
lambdas <- 10^seq(10, -10, by = -.5)
print("calculating gwas coefficients")
datafm<-data.frame(y=scale(y),o=scale(O))
coefGwas=rep(0,length(B))
for(i in 1:length(B)){
  mod=lm(scale(y) ~ scale(O[,i]))
  if(lmp(mod)<(0.05/length(B)))  
   coefGwas[i]=mod$coeff[2]
}
print("correlation between gwas and true coefficients")
corGwas=cor(B,coefGwas)
print(corGwas)
print("calculating summary of fit")

summarygwas=summary(lm(data=data.frame(y=scale(y),x=scale(O)%*%coefGwas),y~x))
r2gwas=summarygwas$adj.r.squared
print(r2gwas)
varEGwas=var(scale(O)%*%coefGwas)
print(varEGwas)
sigmaEGwas=(summarygwas$sigma)**2
print(sigmaEGwas)
mseGwas=sqrt(mean((B-coefGwas)^2))


library(glmnet)
print("calculating blup coefficients")
lm_blup <-cv.glmnet(scale(O), y, alpha = 1)
coefBlup=coef(lm_blup, s = "lambda.min")[-1]
corBlup=cor(B,coefBlup)
print(corBlup)
summaryBlup=summary(lm(data=data.frame(y=y,x=coef(lm_blup, s = "lambda.min")[1]+scale(O)%*%coefBlup),formula=y~x))
r2Blup=summaryBlup$adj.r.squared
print(r2Blup)
mseBlup=sqrt(mean((B-coefBlup)^2))
varEBlup=var(scale(O)%*%coefBlup)
print(varEBlup)
sigmaEBlup=(summaryBlup$sigma)**2

print("calculating lasso coefficients")
lm_lasso <-cv.glmnet(scale(O), y )
coeflasso=coef(lm_lasso, s = "lambda.min")[-1]
corlasso=cor(B,coeflasso)
print(corlasso)
summarylasso=summary(lm(data=data.frame(y=y,x=coef(lm_blup, s = "lambda.min")[1]+scale(O)%*%coeflasso),formula=y~x))
r2lasso=summarylasso$adj.r.squared
print(r2lasso)
mselasso=sqrt(mean((B-coeflasso)^2))
varElasso=var(scale(O)%*%coeflasso)
print(varElasso)
sigmaElasso=(summarylasso$sigma)**2


corOsca<-cor(B,osca.B)
summaryOsca<-summary(lm(data=data.frame(y=y,x=coef(lm_blup, s = "lambda.min")[1]+scale(O)%*%osca.B),formula=y~x))
r2Osca<-summaryOsca$adj.r.squared
mseOsca<-sqrt(mean((B-osca.B)^2))
varEOsca<-osca.pve
sigmaEOsca=(summaryOsca$sigma)**2


varEM=var(scale(M)%*%B)
varEO=var(scale(O)%*%B)


corGRTruth<-sqrt(sum(cov(scale(M)%*%B,R)^2))
corGRB<- sqrt(sum((cov(scale(M)%*%meanBetas,R) )^2))
corGRBlup<-sqrt(sum((cov(scale(M)%*%coefBlup,R) )^2))
corGRGwas<-sqrt(sum((cov(scale(M)%*%coefGwas,R) )^2))
corGRlasso<-sqrt(sum((cov(scale(M)%*%coeflasso,R)))^2)

PR<-prcomp(scale(M))
corPRB<-cor(B-meanBetas,PR$rotation[,1])
corPRBlup<-cor(B-coefBlup,PR$rotation[,1])
corPRGwas<-cor(B-coefGwas,PR$rotation[,1])
corPRlasso<-cor(B-coeflasso,PR$rotation[,1])

list(
     mse_mean_B=mse_mean_B,
     cor_B=cor_B,
     meanMu=meanMu,
     meanSigmaE=meanSigmaE,
     meanSigmaG=meanSigmaG,
     corGwas=corGwas,
     r2gwas=r2gwas,
     mseGwas= mseGwas,
     varEGwas=varEGwas,
     sigmaEGwas=sigmaEGwas,
     corBlup=corBlup,
     r2Blup= r2Blup,
     mseBlup= mseBlup,
     varEBlup= varEBlup,
     sigmaEBlup= sigmaEBlup,
     corlasso=corlasso,
     r2lasso=r2lasso,
     mselasso=mselasso,
     varElasso= varElasso,
     sigmaElasso= sigmaElasso,
     corR=corR,
     covR=covR,
     varEM=varEM,
     varEO=varEO,
     corGRB=corGRB,
     corGRBlup=corGRBlup,
     corGRGwas=corGRGwas,
     corGRlasso=corGRlasso,
     corGRTruth=corGRTruth,
     corPRB=corPRB,
     corPRGwas=corPRGwas,
     corPRBlup=corPRBlup,
     corPRlasso=corPRlasso,
     corOsca=corOsca,
     r2Osca=r2Osca,
     mseOsca=mseOsca,
     varEOsca=varEOsca,
     sigmaEOsca=sigmaEOsca
     )

}

args = commandArgs(trailingOnly=TRUE)
simulation_root=args[1]
simulations=grep("sim_data",list.dirs(path=simulation_root),value=T)
library(parallel)
mc<-makeCluster(length(simulations),type="FORK")
clusterExport(cl=mc,varlist=c("simulations"))
lapply(1:length(simulations),function(i){
        library(bigmemory)
        library(data.table)
	library(coda)
	library(ggmcmc)
        simulation_id=paste(simulations[i],paste("/",strsplit(simulations[i],split="/")[[1]][4],sep=""),sep="")
        print(simulation_id)
        print(simulations[i])
        y=(as.matrix(fread(paste(simulation_id,"_y.dat",sep="")))[,1])
        X=((as.matrix(fread(paste(simulation_id,"_O.dat",sep="")))))
        B=as.matrix(fread(paste(simulation_id,"_B.dat",sep="")))[,1]
	M=as.matrix(fread(paste(simulation_id,"_M.dat",sep="")))
        R=as.matrix(fread(paste(simulation_id,"_R.dat",sep="")))
        osca.pve<-fread(paste(simulation_id,"_O.hsq",sep=""),fill=T)
        osca.pve<-as.numeric(osca.pve[Source=="V(G)",Variance])
        osca.B<-as.matrix(fread(paste(simulation_id,"_O.mlma",sep=""))$b)
        posteriorSummary<-process_posterior(simulations[i],B,y,X,M,R,osca.pve,osca.B)
        save(list="posteriorSummary",file=paste(simulations[i],"/posteriorSummary.RData",sep=""))
        }
        )


