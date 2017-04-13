calH <- function(betas){
  betas. <- betas
  pars. <- pars.cal.t(betas., thetaZs, gammaZs, etaZs, muZs)
  t.t <- matrix(pars.[,1],ncol=1)
  t.g <-matrix(pars.[,2],ncol=1)
  t.e<- matrix(pars.[,3],ncol=1)
  t.m<- matrix(pars.[,4],ncol=1)
  
  H.01 <- H.02 <- rep(0.001,len.et)
  H.diff1 <- H.diff2 <- 1; iterH <- 1
  
  while ((H.diff1>tolerance | H.diff2>tolerance) & iterH <= maxit){
    dem1 =colSums(sweep(risk.1,1, t.t, multiply))
    H.1 = ifelse(dem1==0,0,dN1r/dem1)
    
    dem2 = colSums(sweep(risk.1,1, t.g, multiply))+
      colSums(sweep(risk.1r2d,1, t.e, multiply))+
      colSums(sweep(risk.1nr2d,1, t.m, multiply))
    H.2 = ifelse(dem2==0,0,dNdeath/dem2)
    
    H.diff1 <- sqrt(sum((H.1-H.01)^2))/sqrt(length(H.1))
    H.diff2 <- sqrt(sum((H.2-H.02)^2))/sqrt(length(H.2))
    H.01 <- H.1; H.02<-H.2; iterH <- iterH+1
    if( iterH %% 40 == 0 ) cat(paste(iterH))
  }
  result <- list(H.01, H.02, iterH); names(result) <- c('H1','H2','iterH')
  result
}