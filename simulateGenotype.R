load( methylation matrix)
load(B)
library(stats)
markers<-dim(M)[2]
samples<-dim(M)[1]
X<- matrix(rbinom(markers*samples,2,0.5),nrow=samples,ncol=markers)
X[,B!=0]<-apply(M[,B!=0],MARGIN = 2,function(x){
  cofcor<-runif(1,20,25)
  tmp<-rbinom(samples,2,1/(1+exp(10-cofcor*x)))
  if(sd(tmp)==0)
    tmp<-rbinom(samples,2,0.5)
  tmp
})

B[B!=0]<-rnorm(100,0,sqrt(0.5/100))
Bg<-rep(0,markers)
Bg[B!=0]<-rnorm(100,0,sqrt(0.2/100))

Bg[sample(setdiff(1:markers,which(B!=0)),900)]<-rnorm(900,0,sqrt(0.1/900))
y=scale(M)%*%as.matrix(B)+scale(X)%*%Bg+rnorm(samples,sd=sqrt(0.2))
save(list=c("X","y","B","Bg"),file=)