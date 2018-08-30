
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

process_posterior<- function(chains,B,Bg,y,O,M,X){
#args <- commandArgs()
# Calculate the number of cores
#for some reason, if I do this using lapply I run into heap problems, probably the list of bigmemory objects is not well allocated

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
print("freeing up space")
rm(C1)
gc()

nMethyl<-1:ncol(O)
nG<-(ncol(O)+1):(ncol(O)+ncol(X))
vary<-var(y)
varm<-var(scale(M)%*%B)
varg<-var(scale(X)%*%Bg)

meanSigmaG<-colMeans(do.call(rbind,chainsSigma))[3]
print(meanSigmaG)
meanSigmaPhi<-colMeans(do.call(rbind,chainsSigma))[2]
print(meanSigmaPhi)
PR<-prcomp(t(scale(X)))

datafm<-as.data.frame(cbind(y,PR$rotation[,1:10]))
yadj<-residuals(lm(y~.,data = datafm ))
yadj<-as.matrix(yadj)
print(dim(yadj))
D<-scale(cbind(O,X))

print("calculating gwas coefficients")

coefGwas=rep(0,length(B))

for(i in 1:length(B)){
  mod=lm(scale(yadj) ~ scale(O[,i]))
  if(lmp(mod)<(0.05/length(B)))  
   coefGwas[i]=mod$coeff[2]
}
print("calculating summary of fit")
varEGwasm=var(scale(O)%*%coefGwas)
print(varEGwasm)


library(glmnet)
print("calculating blup coefficients")
lm_blup <-cv.glmnet(D, y, alpha = 1)
coefBlup=coef(lm_blup, s = "lambda.min")[-1]
varEBlupm=var(scale(O)%*%coefBlup[nMethyl])
varEBlupg=var(scale(X)%*%coefBlup[nG])
print(varEBlupm)

print("calculating lasso coefficients")
lm_lasso <-cv.glmnet(D, y )
coeflasso=coef(lm_lasso, s = "lambda.min")[-1]
varElassom=var(scale(O)%*%coeflasso[nMethyl])
varElassog=var(scale(X)%*%coeflasso[nG])
print(varElassom)


list(
     varg=varg,
     varm=varm,
     meanSigmaG=meanSigmaG,
     meanSigmaPhi=meanSigmaPhi,
     varEGwasm=varEGwasm,
     varEBlupm= varEBlupm,
     varEBlupg=varEBlupg,
     varElassom= varElassom, 
     varElassog= varElassog,
     vary=vary
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
        y=(as.matrix(fread(paste(simulation_id,"_y.dat",sep="")))[,-1])
        O=((as.matrix(fread(paste(simulation_id,"_O.dat",sep="")))))
        B=as.matrix(fread(paste(simulation_id,"_B.dat",sep="")))[,-1]
       	M=as.matrix(fread(paste(simulation_id,"_M.dat",sep="")))
	      Bg=as.matrix(fread(paste(simulation_id,"_Bg.dat",sep="")))[,-1]
	      X=((as.matrix(fread(paste(simulation_id,"_X.dat",sep="")))))[,-1]
	      posteriorSummary<-process_posterior(simulations[i],B,Bg,y,O,M,X)
        save(list="posteriorSummary",file=paste(simulations[i],"/posteriorSummary.RData",sep=""))
        }
        )


