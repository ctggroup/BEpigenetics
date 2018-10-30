
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

process_posterior<- function(B,y,O,M,R,DATA_FILE){
 library(lfmm)
  print("calculating lfmm_ridge coefficients")
  #10-fold cross validation
  errs <- lfmm_ridge_CV(Y = scale(y),
                        X = scale(O),
                        n.fold.row = 10,
                        n.fold.col = 10,
                        lambdas = c(1e-10, 1, 1e20),
                        Ks = 5)
  #fit model with lambda for 10-fold cross validation
  mod.lfmm <- lfmm_ridge(Y = scale(y),
                         X = scale(O),
                         K = 5,
                         lambda=errs[errs$err==min(errs$err),"lambda"])
  #association test
  pv <- lfmm_test(Y = scale(y),
                  X = scale(O),
                  lfmm = mod.lfmm,
                  calibrate = "gif")
  
  coeflfmm_ridge=ifelse(pv$calibrated.pvalue<0.05,pv$B,0)
  print("correlation between lfmm_ridge and true coefficients ")
  corlfmm_ridge=cor(B,coeflfmm_ridge)
  print(corlfmm_ridge)
  print("calculating summary of fit ")
  
  summarylfmm_ridge=summary(lm(data=data.frame(y=scale(y),x=scale(O)%*%coeflfmm_ridge),y~x))
  r2lfmm_ridge=summarylfmm_ridge$adj.r.squared
  print("R squared ")
  print( r2lfmm_ridge)
  varElfmm_ridge=r2lfmm_ridge
  print("SigmaE ")
  sigmaElfmm_ridge=(summarylfmm_ridge$sigma)**2
  print(sigmaElfmm_ridge)
  mselfmm_ridge=sqrt(mean((B-coeflfmm_ridge)^2))
  print("MSE")
  print(mselfmm_ridge)
  #now lfmm lasso
  #fit model with lambda for 10-fold cross validation
  print("running lfmm_lasso")
  mod.lfmm <- lfmm_lasso(Y = scale(y),
                         X = scale(O),
                         K = 5)
  #association test
  pv <- lfmm_test(Y = scale(y),
                  X = scale(O),
                  lfmm = mod.lfmm,
                  calibrate = "gif")
  
  coeflfmm_lasso=ifelse(pv$calibrated.pvalue<0.05,pv$B,0)
  print("correlation between gwas and true coefficients")
  corlfmm_lasso=cor(B,coeflfmm_lasso)
  print(corlfmm_lasso)
  print("calculating summary of fit")
  
  summarylfmm_lasso=summary(lm(data=data.frame(y=scale(y),x=scale(O)%*%coeflfmm_lasso),y~x))
  
  r2lfmm_lasso=summarylfmm_lasso$adj.r.squared
  print("R squared lfmm_lasso ")
  print(r2lfmm_lasso)
  varElfmm_lasso=r2lfmm_lasso
  sigmaElfmm_ridge=(summarylfmm_ridge$sigma)**2
  print("sigmaE")
  print(sigmaElfmm_ridge)
  print("MSE")
  mselfmm_lasso=sqrt(mean((B-coeflfmm_lasso)^2))
  print(mselfmm_lasso)
  
  #Now with refactor
  print("Running Refactor")
  source("refactor.R")
  coefRefactor<-rep(0,length(B))
  output <- refactor(DATA_FILE, K, out = "demo_refactor")
  print("calculating coefficients")
   associations_test <- function(O, y, model_append)
  {
    observed_pvalues <- c()
    for (site in 1:ncol(O))
    {
      model <- lm(y ~ O[,site] + model_append)
      
      pvalue <- coef(summary(model))[2,4]
      observed_pvalues[site] = as.numeric(pvalue)
      if(pvalue<0.05/length(B))
        coefRefactor[site]<-model$coeff[2]
    }
    
    return(observed_pvalues) 
   }
   print("correlation between refactor coefficients and true coefficients ")
   corRefactor=cor(B,coefRefactor)
   print(corRefactor)
   print("calculating summary of fit ")
   
   summaryRefactor=summary(lm(data=data.frame(y=y,x=O%*%coefRefactor),y~x))
   print("R squared ")
   r2lRefactor=summaryRefactor$adj.r.squared
   print(r2lRefactor)
   varERefactor=r2Refactor
   print("SigmaE ")
   sigmaERefactor=(summaryRefactor$sigma)**2
   print(sigmaERefactor)
   print("MSE ")
   mseRefactor=sqrt(mean((B-coefRefactor)^2))
   print(mseRefactor)
   
  
  varEM=var(scale(M)%*%B)
  varEO=var(scale(O)%*%B)
  
  
  corGRTruth<-sqrt(sum(cov(scale(M)%*%B,R)^2))
  corGRlfmm_ridge<-sqrt(sum((cov(scale(M)%*%coeflfmm_ridge,R) )^2))
  corGRlfmm_lasso<-sqrt(sum((cov(scale(M)%*%coeflfmm_lasso,R) )^2))
  corGRRefactor<-sqrt(sum((cov(scale(M)%*%coefRefactor,R)))^2)
  
  PR<-prcomp(scale(M))
  corPRlfmm_ridge<-cor(B-coefBlup,PR$rotation[,1])
  corPRlfmm_lasso<-cor(B-coefGwas,PR$rotation[,1])
  corPRRefactor<-cor(B-coeflasso,PR$rotation[,1])
  
  list(
   
    corlfmm_ridge=corlfmm_ridge,
    r2lfmm_ridge=r2lfmm_ridge,
    mselfmm_ridge= mselfmm_ridge,
    varElfmm_ridge=varElfmm_ridge,
    sigmaElfmm_ridge=sigmaElfmm_ridges,
    corlfmm_lasso=corlfmm_lasso,
    r2lfmm_lasso= r2lfmm_lasso,
    mselfmm_lasso= mselfmm_lasso,
    varElfmm_lasso= varElfmm_lasso,
    sigmaElfmm_lasso= sigmaElfmm_lasso,
    corRefactor=corRefactor,
    r2Refactor=r2Refactor,
    mseRefactor=mseRefactor,
    varERefactor= varERefactor,
    sigmaERefactor= sigmaERefactor,
    corGRlfmm_ridge=corGRlfmm_ridge,
    corGRlfmm_lasso=corGRlfmm_lasso,
    corGRRefactor=corGRRefactor,
    corPRlfmm_ridge=corPRlfmm_ridge,
    corPRlfmm_lasso=corPRlfmm_lasso,
    corPRRefactor=corPRRefactor
  )
  
}

args = commandArgs(trailingOnly=TRUE)
simulation_root=args[1]
simulations=grep("sim_data",list.dirs(path=simulation_root),value=T)
library(parallel)
mc<-makeCluster(length(simulations),type="FORK")
clusterExport(cl=mc,varlist=c("simulations"))
parLapply(mc,1:length(simulations),function(i){
  library(bigmemory)
  library(data.table)
  library(coda)
  library(ggmcmc)
  simulation_id=paste(simulations[i],paste("/",strsplit(simulations[i],split="/")[[1]][4],sep=""),sep="")
  print(simulation_id)
  print(simulations[i])
  y=(as.matrix(fread(paste(simulation_id,"_y.dat",sep="")))[,1])
  X=(t(as.matrix(fread(paste(simulation_id,"_O.dat",sep="")))))
  B=as.matrix(fread(paste(simulation_id,"_B.dat",sep="")))[,1]
  M=as.matrix(fread(paste(simulation_id,"_M.dat",sep="")))
  R=as.matrix(fread(paste(simulation_id,"_R.dat",sep="")))
  osca.pve<-fread(paste(simulation_id,"_O.hsq"),fill=T)
  osca.pve<-as.numeric(osca.pve[Source=="V(G)",Variance])
  osca.B<-as.matrix(fread(paste(simulation_id,"_O.mlma"))$b)
  posteriorSummary<-process_posterior(simulations[i],B,y,X,M,R,paste(simulation_id,"_O.dat",sep=""))
  save(list="posteriorSummary",file=paste(simulations[i],"/posteriorSummary_lfmm_ref.RData",sep=""))
}
)
