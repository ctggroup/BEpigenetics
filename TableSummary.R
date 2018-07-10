table.summary<-function(result){
  
  table.summary<-rbind(
  "Markers in model"=unclass(summary(result$perSampleMarkers)),
  "cpg probes in model"=unclass(summary(result$perSampleProbes)),
  "SNP in model"=unclass(summary(result$perSampleSNP)),
  "SNP in mixture 1"=unclass(summary(result$SNPper1)),
  "SNP in mixture 2"=unclass(summary(result$SNPper2)),
  "SNP in mixture 3"=unclass(summary(result$SNPper3)),
  "cpg probes in mixture 1"=unclass(summary(result$Pper1)),
  "cpg probes in mixture 2"=unclass(summary(result$Pper2)),
  "cpg probes in mixture 3"=unclass(summary(result$Pper3)),
  "VE(%) SNP"=unclass(summary(100*result$VEgen)),
  "VE(%) SNP in mixture 1"=unclass(summary(100*result$varSNP1)),
  "VE(%) SNP in mixture 2"=unclass(summary(100*result$varSNP2)),
  "VE(%) SNP in mixture 3"=unclass(summary(100*result$varSNP3)),
  "VE(%) cpg"=unclass(summary(100*result$VEepi)),
  "VE(%) cpg in mixture 1"=unclass(summary(100*result$varP1)),
  "VE(%) cpg in mixture 2"=unclass(summary(100*result$varP2)),
  "VE(%) cpg in mixture 3"=unclass(summary(100*result$varP3)),
  "VE(%) by probes with 95% IP"=unclass(summary(100*result$VE95)),
  "(%) genome mapping to probes in model"=unclass(summary(100*result$genesPerSample)))
   as.data.frame(table.summary)
}
