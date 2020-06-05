library(readr)
HDIofMCMC = function( sampleVec , credMass=0.95 ) {
  # Computes highest density interval from a sample of representative values,
  #   estimated as shortest credible interval.
  # Arguments:
  #   sampleVec
  #     is a vector of representative values from a probability distribution.
  #   credMass
  #     is a scalar between 0 and 1, indicating the mass within the credible
  #     interval that is to be estimated.
  # Value:
  #   HDIlim is a vector containing the limits of the HDI
  sortedPts = sort( sampleVec )
  ciIdxInc = ceiling( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( "CI(l)"=HDImin , "CI(h)"=HDImax )
  return( HDIlim )
}

plot_GO<-function(enrich.df, title, outputFile, PIP, trait){
  enrich.df.agg<-aggregate(enrichment ~ go_id,data=enrich.df,function(x)c(unclass(summary(x)),HDIofMCMC(x),"ROPE"=sum(x<1.5 & x>0.5)/400) )
  #tmp<-as.data.frame(enrich.df.agg$enrichment)
  if(is.matrix(enrich.df.agg$enrichment)){
    tmp<-as.data.frame(enrich.df.agg$enrichment)
  }else{
    tmp<-as.data.frame(do.call(rbind,enrich.df.agg$enrichment))
  }
  enrich.df.smry<-cbind(go_id=enrich.df.agg$go_id,tmp)
  enrich.df$freq<-NULL
  enrich.df.prob<-plyr::count(enrich.df,"go_id")
  enrich.df.prob$prob<-enrich.df.prob$freq/400
  enrich.df.prob$freq<-NULL
  df<-merge(enrich.df.smry,enrich.df.prob)
  
  topGO<-df
  #topGO$group<-ifelse(df$prob>0.95 & df$Mean>2.0,"Enriched","Non-enriched")
  #with ROPE
  topGO$group<-ifelse(df$prob>0.95 & df$ROPE<0.05,"Enriched","Non-enriched")
  ribbon.df<-data.frame(ymin=rep(0.5,100),ymax=rep(1.5,100),x=seq(0,1,length.out = 100))
  p<-ggplot()+  geom_point(data=topGO, aes(x=prob, y=Mean,colour=group)) +ylab("Enrichment") + xlab("Posterior inclusion probability") + ggtitle(paste("Enrichment of terms according to PIP for", toupper(trait)))
  topGO$prob=as.numeric(topGO$prob)
  topGO[,"CI(l)"]<-as.numeric(topGO[,"CI(l)"])
  topGO[,"CI(h)"]<-as.numeric(topGO[,"CI(h)"])
  p<-p+geom_errorbar(data=topGO,aes(x=prob,ymin=topGO[,"CI(l)"],ymax=topGO[,"CI(h)"],colour=group))
  p<-p+geom_ribbon(data=ribbon.df,aes(ymin=ymin, ymax=ymax, x=x, fill = "ROPE"),linetype="solid",alpha=1.0)+scale_fill_manual(values=c( "grey80"), name="")
  ggsave(plot = p,filename = paste(c("Enrichment_",trait,".pdf"), collapse = ""))
  
  
  
  topGO<-df[df$prob>0.95 & df$ROPE<0.05,]
 
  topGO <- as_tibble(topGO) %>% arrange(desc(Mean))
  write_delim(x = topGO,delim = " ", path = paste(c("GO_",trait,".csv"),collapse = ""))#list of go terms for REVIGO
  topGO <- as.data.frame(topGO)
  print(xtable(topGO),file = paste(c("GO_",trait,".tex"),collapse = ""))
  
  enrich.df.Pprobes<-aggregate(P.probes ~ go_id,data=enrich.df,function(x)c(unclass(summary(x)),HDIofMCMC(x)) )
  if(is.matrix(enrich.df.agg$enrichment)){
    tmp<-as.data.frame(enrich.df.Pprobes$P.probes)
  }else{
    tmp<-as.data.frame(do.call(rbind,enrich.df.Pprobes$P.probes))
  }
  enrich.df.smry.Pprobes<-cbind(go_id=enrich.df.Pprobes$go_id,tmp)
  
  enrich.df.PVE<-aggregate(P.VE ~ go_id,data=enrich.df,function(x)c(unclass(summary(x)),HDIofMCMC(x)) )
  if(is.matrix(enrich.df.PVE$P.VE)){
    tmp<-as.data.frame(enrich.df.PVE$P.VE)
  }else{
    tmp<-as.data.frame(do.call(rbind,enrich.df.PVE$P.VE))
  }
  enrich.df.smry.PVE<-cbind(go_id=enrich.df.PVE$go_id,tmp)
  
  
  enrich.df.smry.pxv<-merge(data.frame(go_id=enrich.df.smry.Pprobes$go_id,mean_pprobes=enrich.df.smry.Pprobes$Mean),data.frame(go_id=enrich.df.PVE$go_id,mean_pve=enrich.df.smry.PVE$Mean))
  enrich.df.smry.pxv<-merge(enrich.df.smry.pxv,df)
  enrich.df.smry.pxv$group<-ifelse(enrich.df.smry.pxv$prob>0.95 & enrich.df.smry.pxv$ROPE<0.05,"Enriched","Non-enriched")
  
  plot_repr<-ggplot(enrich.df.smry.pxv, aes(x=mean_pprobes,y=mean_pve,colour=group)) +geom_point()+geom_line(data=enrich.df.smry.pxv,aes(x=mean_pprobes,y=mean_pprobes,color="Expected")) +ggtitle(paste("Enrichment of terms according to representation in the model"),toupper(trait)) +ylab("Mean proportion of PVE per term") +xlab("Mean proportion of probes mapping to term")+theme_light()
  ggsave(paste(c("enrichmentxproportion_", trait, ".pdf"), collapse = ""), plot_repr)
  plot_repr
  
}
