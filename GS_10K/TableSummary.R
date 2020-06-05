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
  HDIlim = c( "95% HDI(low)"=HDImin , "95% HDI(high)"=HDImax )
  return( HDIlim )
}

hdi.summary<-function(x){
  c(unclass(summary(x)),HDIofMCMC(x))
}

table.summary<-function(result){
  
  table.summary<-rbind(
    "Markers in model"=hdi.summary(result$perSampleMarkers),
    "cpg probes in model"=hdi.summary(result$perSampleProbes),
    "SNP in model"=hdi.summary(result$perSampleSNP),
    "SNP in mixture 1"=hdi.summary(result$SNPper1),
    "SNP in mixture 2"=hdi.summary(result$SNPper2),
    "SNP in mixture 3"=hdi.summary(result$SNPper3),
    "cpg probes in mixture 1"=hdi.summary(result$Pper1),
    "cpg probes in mixture 2"=hdi.summary(result$Pper2),
    "cpg probes in mixture 3"=hdi.summary(result$Pper3),
    "VE(%) SNP"=hdi.summary(100*result$VEgen),
    "VE(%) SNP in mixture 1"=hdi.summary(100*result$varSNP1),
    "VE(%) SNP in mixture 2"=hdi.summary(100*result$varSNP2),
    "VE(%) SNP in mixture 3"=hdi.summary(100*result$varSNP3),
    "VE(%) cpg"=hdi.summary(100*result$VEepi),
    "VE(%) cpg in mixture 1"=hdi.summary(100*result$varP1),
    "VE(%) cpg in mixture 2"=hdi.summary(100*result$varP2),
    "VE(%) cpg in mixture 3"=hdi.summary(100*result$varP3),
    "VE(%) by probes with 95% IP"=hdi.summary(100*result$VE95),
    "(%) genome mapping to probes in model"=hdi.summary(100*result$genesPerSample))
  as.data.frame(table.summary)
}
