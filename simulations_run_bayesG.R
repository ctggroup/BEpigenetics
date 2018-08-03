args = commandArgs(trailingOnly=TRUE)
simulation_root=args[1]
simulations=grep("sim_data",list.dirs(path=simulation_root),value=T)
library(BayesRRcpp)
library(parallel)
mc<-makeCluster(length(simulations),type="FORK")
clusterExport(cl=mc,varlist=c("simulations"))
parLapply(mc,1:length(simulations),function(i){
  library(data.table)
  library(BayesRRcpp)
  simulation_id=paste(simulations[i],paste("/",strsplit(simulations[i],split="/")[[1]][4],sep=""),sep="")
  chain_name=paste(paste(simulations[i],"/C1",sep=""),".csv",sep="")
  print(simulation_id)
  y=scale(as.matrix(fread(paste(simulation_id,"_y.dat",sep="")))[,1])
  O=scale((as.matrix(fread(paste(simulation_id,"_O.dat",sep="")))))
  X=scale(as.matrix(fread(paste(simulation_id,"_X.dat",sep=""))))
  cva=matrix(data=c(c(0.1,0.01,0.001,0.0001),c(0.01,0.001,0.0001,0.00001)),nrow=2,byrow = T)
  gAssign=c(rep(0,ncol(O)),rep(1,ncol(X)))
  BayesRSamplerV2Groups(chain_name, 1, 20000,10000 , 10,cbind(O,X), y, 0.001, 0.001,0.001, 0.001, 0.001, cva=cva,2,gAssign = gAssign)
  
  
})
