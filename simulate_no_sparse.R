args = commandArgs(trailingOnly=TRUE)
simulation_root<-args[1]
simulations<-unique(dirname(list.files(path=simulation_root,pattern=".*dat",full.names=T,recursive=T)))
library(stats)
library(data.table)
lapply(simulations,function(simulation_dir){
   print(simulation_dir) 
   m_file=list.files(path=simulation_dir,pattern=".*_M.dat",full.names=T,recursive=T)
   print(m_file)
   print("reading methylation matrix")
   M <- as.matrix(fread(m_file))
   b_file=list.files(path=simulation_dir,pattern=".*_B.dat",full.names=T,recursive=T)
   print(b_file) 
   B <- as.matrix(fread(b_file))
   y_file=list.files(path=simulation_dir,pattern=".*_y.dat",full.names=T,recursive=T)
   print(y_file)
   
   markers<-dim(M)[2]
   print(markers)
   samples<-dim(M)[1]
   print(samples)
   print("simulating effects")
   nonZmarker<-sample(1:markers,size=1000)
   B<-rep(0,markers)
   B[nonZmarker]<-rnorm(1000,0,sqrt(0.6/1000))
   
   print("simulating phenotype")
   y=scale(M)%*%as.matrix(B)+rnorm(samples,sd=sqrt(0.4))
   print(head(y))
   print("saving simulations phenotype")
   write.table(y,file=y_file,col.names=F,sep="\t")
    print("saving simulations m. effects")
   write.table(B,file=b_file,col.names=F,sep="\t")
    print("saving simulations genetic effects")
}
)
