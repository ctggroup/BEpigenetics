args = commandArgs(trailingOnly=TRUE)
simulation_root=args[1]
simulations=grep("sim_data",list.dirs(path=simulation_root),value=T)
print(simulations)
library(parallel)
mc<-makeCluster(min(length(simulations),20))

sim_list=parLapply(mc,simulations,function(x){res<-try(load(paste(x,"/posteriorSummary.RData",sep="")))
  if(inherits(res, "try-error")){return(list(mse=NA,
	          corB=NA,
        	  adjR=NA,
                  sigmaG=NA,
                  sigmaE=NA,
                  mseGwas=NA,
                  corGwas=NA,
                  adjRGwas=NA,
                  varEGwas=NA,
                  sigmaEGwas=NA,
                  mseBlup=NA,
                  corBlup=NA,
                  adjRBlup=NA,
                  varEBlup=NA,
                  sigmaEBlup=NA,
                  mseLasso=NA,
                  corLasso=NA,
                  adjRLasso=NA,
                  varELasso=NA,
                  sigmaELasso=NA,
                  corR=NA,
                  covR=NA,
                  varEM=NA,
                  corGRB=NA,
                  corGRBlup=NA,
                  corGRGwas=NA,
                  corGRlasso=NA,
                   corGRTruth=NA,
                   corPRBlup=NA,
		   corPRB=NA,
		   corPRGwas=NA,
	           corPRlasso=NA,
		   corOsca<-NA,
		   summaryOsca<-NA,
		   r2Osca<-NA,
		   mseOsca<-NA,
		   varEOsca<-NA,
		   sigmaEOsca=NA
				))}
     list(mse=posteriorSummary$mse_mean_B,
	  corB=posteriorSummary$cor_B,
          adjR=(posteriorSummary$lm_summary)$adj.r.squared,
          sigmaG=posteriorSummary$meanSigmaG,
          sigmaE=posteriorSummary$meanSigmaE,
          mseGwas=posteriorSummary$mseGwas,
	  corGwas=posteriorSummary$corGwas,
	  adjRGwas=posteriorSummary$r2gwas,
          varEGwas=posteriorSummary$varEGwas,
          sigmaEGwas=posteriorSummary$sigmaEGwas,
          mseBlup=posteriorSummary$mseBlup,
	  corBlup=posteriorSummary$varEBlup,	
	  adjRBlup=posteriorSummary$r2Blup,
	  varEBlup=posteriorSummary$varEBlup,
	  sigmaEBlup=posteriorSummary$sigmaEBlup,
          mseLasso=posteriorSummary$mselasso,
	  corLasso=posteriorSummary$corlasso,
	  adjRLasso=posteriorSummary$r2lasso,
	  varELasso=posteriorSummary$varElasso,
	  sigmaELasso=posteriorSummary$sigmaElasso,
	  corR=posteriorSummary$corR,
	  covR=posteriorSummary$covR,   
	  varEM=posteriorSummary$varEM,
          corGRB=posteriorSummary$corGRB,
          corGRBlup=posteriorSummary$corGRBlup,
     	  corGRGwas=posteriorSummary$corGRGwas,     
          corGRlasso=posteriorSummary$corGRlasso,
	  corGRTruth=posteriorSummary$corGRTruth,
          corPRBlup=posteriorSummary$corPRBlup,
                   corPRB=posteriorSummary$corPRB,
                   corPRGwas=posteriorSummary$corPRGwas,
                   corPRlasso=posteriorSummary$corPRlasso,
	  corOsca<-posteriorSummary$corOsca,
	  summaryOsca<-posteriorSummary$summaryOsca,
	  r2Osca<-posteriorSummary$r2Osca,
	  mseOsca<-posteriorSummary$mseOsca,
	  varEOsca<-posteriorSummary$varEOsca,
	  sigmaEOsca=posteriorSummary$sigmaEOsca
	)})
sim_df<-do.call(rbind,sim_list)
rownames(sim_df)<-simulations
save(list="sim_df",file=paste(simulation_root,"/simulationPlot.RData",sep=""))