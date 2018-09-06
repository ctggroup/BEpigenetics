
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

process_posterior<- function(B,osca.pve,osca.B){


corOsca<-cor(B,osca.B)
mseOsca<-sqrt(mean((B-osca.B)^2))
varEOsca<-osca.pve



list(
     corOsca=corOsca,
     mseOsca=mseOsca,
     varEOsca=varEOsca
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
        B=as.matrix(fread(paste(simulation_id,"_B.dat",sep="")))[,1]
        osca.pve<-fread(paste(simulation_id,"_O.hsq",sep=""),fill=T)
        osca.pve<-as.numeric(osca.pve[osca.pve$Source=="V(O)/Vp","Variance"])
        osca.B<-as.matrix(fread(paste(simulation_id,"_O.mlma",sep=""))$b)
        oscaSummary<-process_posterior(B,osca.pve,osca.B)
        save(list="oscaSummary",file=paste(simulations[i],"/oscaSummary.RData",sep=""))
        }
        )


