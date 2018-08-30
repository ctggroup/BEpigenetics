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

plot_GO_bmi<-function(){
library(readr)
load("~/Documents/Methylation/SNP+p_bmi_c/enrichmentGO.RData")
enrich.bmi <- enrich 
enrich.df<- do.call(rbind,enrich.bmi)
enrich.df.agg<-aggregate(enrichment ~ go_id,data=enrich.df,function(x)c(unclass(summary(x)),HDIofMCMC(x),"ROPE"=sum(x<1.5 & x>0.5)/1004) )
#tmp<-as.data.frame(enrich.df.agg$enrichment)
if(is.matrix(enrich.df.agg$enrichment)){
  tmp<-as.data.frame(enrich.df.agg$enrichment)
}else{
  tmp<-as.data.frame(do.call(rbind,enrich.df.agg$enrichment))
}
enrich.df.smry<-cbind(go_id=enrich.df.agg$go_id,tmp)
enrich.df$freq<-NULL
enrich.df.prob<-plyr::count(enrich.df,"go_id")
enrich.df.prob$prob<-enrich.df.prob$freq/1004
enrich.df.prob$freq<-NULL
df<-merge(enrich.df.smry,enrich.df.prob)

topGO<-df
#topGO$group<-ifelse(df$prob>0.95 & df$Mean>2.0,"Enriched","Non-enriched")
#with ROPE
topGO$group<-ifelse(df$prob>0.95 & df$ROPE<0.05,"Enriched","Non-enriched")
ribbon.df<-data.frame(ymin=rep(0.5,100),ymax=rep(1.5,100),x=seq(0,1,length.out = 100))
p<-ggplot()+  geom_point(data=topGO, aes(x=prob, y=Mean,colour=group)) +ylab("Enrichment") + xlab("Posterior inclusion probability") + ggtitle("Enrichment of terms according to PIP for BMI")
topGO$prob=as.numeric(topGO$prob)
topGO[,"CI(l)"]<-as.numeric(topGO[,"CI(l)"])
topGO[,"CI(h)"]<-as.numeric(topGO[,"CI(h)"])
p<-p+geom_errorbar(data=topGO,aes(x=prob,ymin=topGO[,"CI(l)"],ymax=topGO[,"CI(h)"],colour=group))
p<-p+geom_ribbon(data=ribbon.df,aes(ymin=ymin, ymax=ymax, x=x, fill = "ROPE"),linetype="solid",alpha=1.0)+scale_fill_manual(values=c( "grey80"), name="")
ggsave(plot=p,filename="Enrichment_bmi.pdf")
#topGO<-df[df$prob>0.95 & df$Mean>2.0,]
topGO<-df[df$prob>0.95 & df$ROPE<0.05,]
#topGO.output<-topGO[,c("go_id","Mean")]
#write.csv(topGO.output,"GObmi.csv",row.names = F,sep = " ",quote = F)
write_delim(x = topGO,delim = " ",path="GObmi.csv")#list of go terms for REVIGO
print(xtable(topGO),file="GObmi.tex")
enrich.df.Pprobes<-aggregate(P.probes ~ go_id,data=enrich.df,function(x)c(unclass(summary(x)),HDIofMCMC(x)) )
#tmp<-as.data.frame(enrich.df.Pprobes$P.probes)
if(is.matrix(enrich.df.agg$enrichment)){
  tmp<-as.data.frame(enrich.df.Pprobes$P.probes)
}else{
  tmp<-as.data.frame(do.call(rbind,enrich.df.Pprobes$P.probes))
}
enrich.df.smry.Pprobes<-cbind(go_id=enrich.df.Pprobes$go_id,tmp)

enrich.df.PVE<-aggregate(P.VE ~ go_id,data=enrich.df,function(x)c(unclass(summary(x)),HDIofMCMC(x)) )
#tmp<-as.data.frame(enrich.df.PVE$P.VE)
if(is.matrix(enrich.df.PVE$P.VE)){
  tmp<-as.data.frame(enrich.df.PVE$P.VE)
}else{
  tmp<-as.data.frame(do.call(rbind,enrich.df.PVE$P.VE))
}
enrich.df.smry.PVE<-cbind(go_id=enrich.df.PVE$go_id,tmp)


enrich.df.smry.pxv<-merge(data.frame(go_id=enrich.df.smry.Pprobes$go_id,mean_pprobes=enrich.df.smry.Pprobes$Mean),data.frame(go_id=enrich.df.PVE$go_id,mean_pve=enrich.df.smry.PVE$Mean))
enrich.df.smry.pxv<-merge(enrich.df.smry.pxv,df)
#enrich.df.smry.pxv$group<-ifelse(enrich.df.smry.pxv$prob>0.95 & enrich.df.smry.pxv$Mean>2.0,"Enriched","Non-enriched")
enrich.df.smry.pxv$group<-ifelse(enrich.df.smry.pxv$prob>0.95 & enrich.df.smry.pxv$ROPE<0.05,"Enriched","Non-enriched")

plot_repr_bmi<-ggplot(enrich.df.smry.pxv, aes(x=mean_pprobes,y=mean_pve,colour=group)) +geom_point()+geom_line(data=enrich.df.smry.pxv,aes(x=mean_pprobes,y=mean_pprobes,color="Expected")) +ggtitle("Enrichment of terms according to representation in the model BMI") +ylab("Mean proportion of PVE per term") +xlab("Mean proportion of probes mapping to term")+theme_light()
ggsave("enrichmentxproportion_bmi.pdf",plot_repr_bmi)
plot_repr_bmi
}


####################################################################
plot_GO_smk<-function(){
load("~/Documents/Methylation/SNP+p_smk/enrichmentGO.RData")
enrich.smk <- enrich 
enrich.df<- do.call(rbind,enrich.smk)
enrich.df.agg<-aggregate(enrichment ~ go_id,data=enrich.df,function(x)c(unclass(summary(x)),HDIofMCMC(x),"ROPE"=sum(x<1.5 & x>0.5)/1004) )
#tmp<-as.data.frame(enrich.df.agg$enrichment)
if(is.matrix(enrich.df.agg$enrichment)){
  tmp<-as.data.frame(enrich.df.agg$enrichment)
}else{
  tmp<-as.data.frame(do.call(rbind,enrich.df.agg$enrichment))
}
enrich.df.smry<-cbind(go_id=enrich.df.agg$go_id,tmp)
head(enrich.df.smry[order(enrich.df.smry$Mean,decreasing = T),])
dim(enrich.df.smry[enrich.df.smry$Mean>1,])
dim(enrich.df.smry[enrich.df.smry$Mean<1,])
enrich.df$freq<-NULL
enrich.df.prob<-plyr::count(enrich.df,"go_id")
enrich.df.prob$prob<-enrich.df.prob$freq/1004
enrich.df.prob$freq<-NULL
df<-merge(enrich.df.smry,enrich.df.prob)


topGO<-df
#topGO$group<-ifelse(df$prob>0.95 & df$Mean>2.0,"Enriched","Non-enriched")
topGO$group<-ifelse(df$prob>0.95 & df$ROPE<0.05,"Enriched","Non-enriched")
ribbon.df<-data.frame(ymin=rep(0.5,100),ymax=rep(1.5,100),x=seq(0,1,length.out = 100))
p<-ggplot()+  geom_point(data=topGO, aes(x=prob, y=Mean,colour=group)) +ylab("Enrichment") + xlab("Posterior inclusion probability") + ggtitle("Enrichment of terms according to PIP for BMI")
topGO$prob=as.numeric(topGO$prob)
topGO[,"CI(l)"]<-as.numeric(topGO[,"CI(l)"])
topGO[,"CI(h)"]<-as.numeric(topGO[,"CI(h)"])
p<-p+geom_errorbar(data=topGO,aes(x=prob,ymin=topGO[,"CI(l)"],ymax=topGO[,"CI(h)"],colour=group))
p<-p+geom_ribbon(data=ribbon.df,aes(ymin=ymin, ymax=ymax, x=x, fill = "ROPE"),linetype="solid",alpha=1.0)+scale_fill_manual(values=c( "grey80"), name="")

ggsave(filename="Enrichment_smk.pdf")

#topGO<-df[df$prob>0.95 & df$Mean>2.0,]
topGO<-df[df$prob>0.95 & df$ROPE<0.05,]

#topGO.output<-topGO[,c("go_id","Mean")]
#write.csv(topGO.output,"GObmi.csv",row.names = F,sep = " ",quote = F)
write_delim(x = topGO,delim = " ",path="GOsmk.csv") #list of go terms for REVIGO
print(xtable(topGO),file="GOsmk.tex")


enrich.df.Pprobes<-aggregate(P.probes ~ go_id,data=enrich.df,function(x)c(unclass(summary(x)),HDIofMCMC(x)))
#tmp<-as.data.frame(enrich.df.Pprobes$P.probes)
if(is.matrix(enrich.df.agg$enrichment)){
  tmp<-as.data.frame(enrich.df.Pprobes$P.probes)
}else{
  tmp<-as.data.frame(do.call(rbind,enrich.df.Pprobes$P.probes))
}
enrich.df.smry.Pprobes<-cbind(go_id=enrich.df.Pprobes$go_id,tmp)

enrich.df.PVE<-aggregate(P.VE ~ go_id,data=enrich.df,function(x)c(unclass(summary(x)),HDIofMCMC(x)) )
#tmp<-as.data.frame(enrich.df.PVE$P.VE)
if(is.matrix(enrich.df.agg$enrichment)){
  tmp<-as.data.frame(enrich.df.PVE$P.VE)
}else{
  tmp<-as.data.frame(do.call(rbind,enrich.df.PVE$P.VE))
}
enrich.df.smry.PVE<-cbind(go_id=enrich.df.PVE$go_id,tmp)


enrich.df.smry.pxv<-merge(data.frame(go_id=enrich.df.smry.Pprobes$go_id,mean_pprobes=enrich.df.smry.Pprobes$Mean),data.frame(go_id=enrich.df.PVE$go_id,mean_pve=enrich.df.smry.PVE$Mean))
enrich.df.smry.pxv<-merge(enrich.df.smry.pxv,df)
#enrich.df.smry.pxv$group<-ifelse(enrich.df.smry.pxv$prob>0.95 & enrich.df.smry.pxv$Mean>2.0,"Enriched","Non-enriched")
enrich.df.smry.pxv$group<-ifelse(enrich.df.smry.pxv$prob>0.95 & enrich.df.smry.pxv$ROPE<0.05,"Enriched","Non-enriched")

plot_repr_smk<-ggplot(enrich.df.smry.pxv, aes(x=mean_pprobes,y=mean_pve,colour=group)) +geom_point()+geom_line(data=enrich.df.smry.pxv,aes(x=mean_pprobes,y=mean_pprobes,color="Expected"))+ggtitle("Enrichment of terms according to representation in the model smoking") +ylab("Mean proportion of PVE per term") +xlab("Mean proportion of probes mapping to term")+theme_light()
ggsave("enrichmentxproportion_smk.pdf",plot=plot_repr_smk)
plot_repr_smk
}




plot_tissue_enrichment<-function(input,output,title){
  load(input)
  enR.df<- do.call(rbind,enR)
  enR.df$freq<-NULL
  enR.df.prob<-plyr::count(enR.df,"variable") 
  enR.df.prob$prob<-enR.df.prob$freq/1004
  enR.df.prob$freq<-NULL



  enR.df.agg<-aggregate(enrichment ~ variable,data=enR.df,function(x)c(unclass(summary(x)),HDIofMCMC(x),"ROPE"=sum(x<1.5 & x>0.5)/1004) )
  if(is.matrix(enR.df.agg$enrichment)){
    tmp<-as.data.frame(enR.df.agg$enrichment)
  }else{
    tmp<-as.data.frame(do.call(rbind,enR.df.agg$enrichment))
  }
  enrich.df.smry<-cbind(tissue=enR.df.agg$variable,tmp)
  

 enrich.df.Pprobes<-aggregate(P.probes ~ variable,data=enR.df,function(x)c(unclass(summary(x)),HDIofMCMC(x)))
 #tmp<-as.data.frame(enrich.df.Pprobes$P.probes)
 if(is.matrix(enrich.df.Pprobes$P.probes)){
   tmp<-as.data.frame(enrich.df.Pprobes$P.probes)
 }else{
   tmp<-as.data.frame(do.call(rbind,enrich.df.Pprobes$P.probes))
 }
 
 enrich.df.smry.Pprobes<-cbind(tissue=enrich.df.Pprobes$variable,tmp)
  
  enrich.df.PVE<-aggregate(P.VE ~ variable,data=enR.df,function(x)c(unclass(summary(x)),HDIofMCMC(x)) )
 # tmp<-as.data.frame(enrich.df.PVE$P.VE)
  if(is.matrix(enrich.df.PVE$P.VE)){
    tmp<-as.data.frame(enrich.df.PVE$P.VE)
  }else{
    tmp<-as.data.frame(do.call(rbind,enrich.df.PVE$P.VE))
  }
  enrich.df.smry.PVE<-cbind(tissue= enrich.df.PVE$variable,tmp)
  
  
  enrich.df.smry.pxv<-merge(data.frame(
    tissue=enrich.df.smry.Pprobes$tissue,mean_pprobes=enrich.df.smry.Pprobes$Mean),data.frame(
      tissue=enrich.df.smry.PVE$tissue,mean_pve=enrich.df.smry.PVE$Mean))
  df<-merge(enrich.df.smry,enR.df.prob,by.x="tissue",by.y="variable")
  topGO<-df
  #topGO$group<-ifelse(df$prob>0.95 & df$Mean>2.0,1,2)
  #lets see with the ROPE
  topGO$group<-ifelse(df$prob>0.95 & df$ROPE<0.05,"Enriched","Non-enriched")
  ribbon.df<-data.frame(ymin=rep(0.5,100),ymax=rep(1.5,100),x=seq(0,1,length.out = 100))
  p<-ggplot()
  p<-p+ geom_point(data=topGO, aes(x=prob, y=Mean,colour=group)) +ylab("Enrichment") + xlab("Posterior inclusion probability") + ggtitle(title)+theme_light()
 
  topGO$prob=as.numeric(topGO$prob)
  topGO[,"CI(l)"]<-as.numeric(topGO[,"CI(l)"])
  topGO[,"CI(h)"]<-as.numeric(topGO[,"CI(h)"])
  p<-p+geom_errorbar(data=topGO,aes(x=prob,ymin=topGO[,"CI(l)"],ymax=topGO[,"CI(h)"],colour=group))
  p<-p+geom_ribbon(data=ribbon.df,aes(ymin=ymin, ymax=ymax, x=x, fill = "ROPE"), alpha = 0.3)+scale_fill_manual(values=c( "grey80"), name="")
  ggsave(plot = p,filename = output[1]) 
  enrich.df.smry.pxv<-merge(enrich.df.smry.pxv,df)
  #with rope
  #enrich.df.smry.pxv$group<-ifelse(enrich.df.smry.pxv$prob>0.95 & (enrich.df.smry.pxv$Mean>2.0|enrich.df.smry.pxv$Mean<0.5),1,2)
  enrich.df.smry.pxv$group<-ifelse(enrich.df.smry.pxv$prob>0.95 & enrich.df.smry.pxv$ROPE<0.05,"Enriched","Non-enriched")
  library(ggplot2)
  q<-ggplot(enrich.df.smry.pxv, aes(x=mean_pprobes,y=mean_pve,colour=group)) +geom_point()+geom_line(data=enrich.df.smry.pxv,aes(x=mean_pprobes,y=mean_pprobes,color="Expected"))+ylab("Mean proportion of PVE per term") +xlab("Mean proportion of probes mapping to term")+theme_light() + ggtitle(title)
  ggsave(plot=q,filename = output[2])
  write_delim(x = enrich.df.smry.pxv,delim = " ",path=output[3])
  print(xtable(enrich.df.smry.pxv),file=output[4])
  list(table=enrich.df.smry.pxv, plot.prob=p, plot.prop=q)
}



en_bmi_o<-plot_tissue_enrichment("../SNP+p_bmi_c/enrichmentGTEover.RData",c("probxenrich_GTE_o_bmi.pdf","propxpve_GTE_o_bmi.pdf","GTE_o_bmi.csv","GTE_o_bmi.tex"), "Enrichment of GTEX over expressed genes for BMI")

en_bmi_u<-plot_tissue_enrichment("../SNP+p_bmi_c/enrichmentGTEunder.RData",c("probxenrich_GTE_u_bmi.pdf","propxpve_GTE_u_bmi.pdf","GTE_u_bmi.csv","GTE_u_bmi.tex"), "Enrichment of GTEX repressed genes for BMI")


en_smk_o<-plot_tissue_enrichment("../SNP+p_smk/enrichmentGTEover.RData",c("probxenrich_GTE_o_smk.pdf","propxpve_GTE_o_smk.pdf","GTE_o_smk.csv","GTE_o_smk.tex"), "Enrichment of GTEX over expressed genes for smoking")

en_smk_u<-plot_tissue_enrichment("../SNP+p_smk/enrichmentGTEunder.RData",c("probxenrich_GTE_u_smk.pdf","propxpve_GTE_u_smk.pdf","GTE_u_bmi.csv","GTE_u_smk.tex"),"Enrichment of GTEX repressed genes for BMI")

order.table<-function(enrichment.table){enrichment.table[order(enrichment.table$Mean,decreasing = T),]}

order.table(en_bmi_o$table)
order.table(en_bmi_u$table)
order.table(en_smk_u$table)
order.table(en_smk_o$table)



en_bmi_o_depict<-plot_tissue_enrichment("../SNP+p_bmi_c/enrichmentDEPICTover.RData",c("probxenrich_DEPICT_o_bmi.pdf","propxpve_DEPICT_o_bmi.pdf","DEPICT_o_bmi.csv","DEPICT_o_bmi.tex"), "Enrichment of DEPICT over expressed genes for BMI")
order.table(en_bmi_o_depict$table)

en_bmi_u_depict<-plot_tissue_enrichment("../SNP+p_bmi_c/enrichmentDEPICTunder.RData",c("probxenrich_DEPICT_u_bmi.pdf","propxpve_DEPICT_u_bmi.pdf","DEPICT_u_bmi.csv","DEPICT_u_bmi.tex"), "Enrichment of DEPICT repressed genes for BMI")
order.table(en_bmi_u_depict$table)


en_smk_o_depict<-plot_tissue_enrichment("../SNP+p_smk/enrichmentDEPICTover.RData",c("probxenrich_DEPICT_o_smk.pdf","propxpve_DEPICT_o_smk.pdf","DEPICT_o_smk.csv","DEPICT_o_smk.tex"), "Enrichment of DEPICT over expressed genes for smoking")
order.table(en_smk_o_depict$table)


en_smk_u_depict<-plot_tissue_enrichment("../SNP+p_smk/enrichmentDEPICTunder.RData",c("probxenrich_DEPICT_u_smk.pdf","propxpve_DEPICT_u_smk.pdf","DEPICT_u_smk.csv","DEPICT_u_smk.tex"), "Enrichment of DEPICT repressed genes for smoking")
order.table(en_smk_u_depict$table)


#bmi.plots<-grid.arrange(en_bmi_o$plot.prob,en_bmi_o$plot.prop,en_bmi_u$plot.prob,en_bmi_u$plot.prop,en_bmi_o_depict$plot.prob,en_bmi_o_depict$plot.prop,en_bmi_u_depict$plot.prob,en_bmi_u_depict$plot.prop,ncol=2)
#
bmi.plots<-list(en_bmi_o$plot.prob,en_bmi_o$plot.prop,en_bmi_u$plot.prob,en_bmi_u$plot.prop,en_bmi_o_depict$plot.prob,en_bmi_o_depict$plot.prop,en_bmi_u_depict$plot.prob,en_bmi_u_depict$plot.prop)
bmi.plots.aux<-lapply(1:length(bmi.plots),function(i){ arrangeGrob(bmi.plots[[i]], top = textGrob(letters[i], x = unit(0, "npc")
                                                                                                   , y   = unit(1, "npc"), just=c("left","top"),
                                                                                                   gp=gpar(col="black", fontsize=18, fontfamily="Helvetica")))})
bmi.plots<-grid.arrange(grobs=bmi.plots.aux,ncol=2)

ggsave("bmi_gtex_depict.pdf",plot=bmi.plots,scale=2)




smk.plots<-list(en_smk_o$plot.prob,en_smk_o$plot.prop,en_smk_u$plot.prob,en_smk_u$plot.prop,en_smk_o_depict$plot.prob,en_smk_o_depict$plot.prop,en_smk_u_depict$plot.prob,en_smk_u_depict$plot.prop)
smk.plots.aux<-lapply(1:length(smk.plots),function(i){ arrangeGrob(smk.plots[[i]], top = textGrob(letters[i], x = unit(0, "npc")
                                                                                                  , y   = unit(1, "npc"), just=c("left","top"),
                                                                                                  gp=gpar(col="black", fontsize=18, fontfamily="Helvetica")))})
smk.plots<-grid.arrange(grobs=smk.plots.aux,ncol=2)
ggsave("smk_gtex_depict.pdf",plot=smk.plots,scale=2)


