args = commandArgs(trailingOnly=TRUE)
simulation_root=args[1]
simulations=grep("sim_data",list.dirs(path=simulation_root),value=T)
print(simulations)
library(parallel)
mc<-makeCluster(min(length(simulations),20))

sim_list=parLapply(mc,simulations,function(x){res<-try(load(paste(x,"/oscaSummary.RData",sep="")))
  if(inherits(res, "try-error")){return(list(
		   corOsca=NA,
		   mseOsca=NA,
		   varEOsca=NA
				))}
     list(
	  corOsca=oscaSummary$corOsca,
	  mseOsca=oscaSummary$mseOsca,
	  varEOsca=oscaSummary$varEOsca
	)})
sim_df<-do.call(rbind,sim_list)
rownames(sim_df)<-simulations
save(list="sim_df",file=paste(simulation_root,"/oscaPlot.RData",sep=""))
