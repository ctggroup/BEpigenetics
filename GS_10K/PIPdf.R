PIPdf <- function(posteriorSummary, flatAnno, filename){
  result <- posteriorSummary
  genes <- merge(data.frame(cpg=names(result$PIP),PIP=as.numeric(result$PIP)),flatAnno,by.x="cpg", by.y="cpg")
  genes <- genes[order(genes$PIP,decreasing = T),]
  genes <- genes[genes$PIP>0.5,]
  print(xtable(genes),file = filename)
  
  data.frame("0<5%"=length(which(result$PIP>0 & result$PIP<0.05)),'5%<50%'=length(which(result$PIP>0.05 & result$PIP<0.5)),'50%<95%'=length(which(result$PIP>0.5 & result$PIP<0.95)),'95%>'=length(which(result$PIP>0.95)))
  
  
}