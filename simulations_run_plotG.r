args = commandArgs(trailingOnly=TRUE)
simulation_root=args[1]
simulations=grep("sim_data",list.dirs(path=simulation_root),value=T)
print(simulations)
library(parallel)
mc<-makeCluster(min(length(simulations),20))

sim_list=parLapply(mc,simulations,function(x){res<-try(load(paste(x,"/posteriorSummary.RData",sep="")))
  if(inherits(res, "try-error")){return(list(meanSigmaG=NA,
     meanSigmaPhi=NA,
     varEGwasm=NA,
     varEBlupm=NA,
     varEBlupg=NA,
     varElassom=NA,
     varElassog=NA,
     varg=NA,
     varm=NA,
     vary=NA

				))}
     list(meanSigmaG=posteriorSummary$meanSigmaG,
     meanSigmaPhi=posteriorSummary$meanSigmaPhi,
     varEGwasm=posteriorSummary$varEGwasm,
     varEBlupm= posteriorSummary$varEBlupm,
     varEBlupg=posteriorSummary$varEBlupg,
     varElassom= posteriorSummary$varElassom,
     varElassog= posteriorSummary$varElassog,
     vary=posteriorSummary$vary,
     varm=posteriorSummary$varm,
     varg=posteriorSummary$varg
	)})
sim_df<-do.call(rbind,sim_list)
rownames(sim_df)<-simulations
save(list="sim_df",file=paste(simulation_root,"/simulationPlot.RData",sep=""))
