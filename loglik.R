# Function to calculate -loglikelihood: 
loglik <- function(betas){
  #cat(betas);cat("\n")
  betas. <- c(betas)
  pars. <- pars.cal.t(betas., thetaZs, gammaZs, etaZs, muZs)
  t.t <- matrix(pars.[,1],ncol=1)
  t.g <-matrix(pars.[,2],ncol=1)
  t.e<- matrix(pars.[,3],ncol=1)
  t.m<- matrix(pars.[,4],ncol=1)
  
  #Find baseline hazard
  Hfit=calH(betas.)
  H.w1 <- Hfit$H1
  H.w2 <- Hfit$H2    
  
  logl=0
  for (i in 1:len.et){
    THETA_1r = t.t
    THETA_1d = t.g
    THETA_1r2d = t.e
    THETA_1nr2d = t.m
    
    S1t = sum(multiply(THETA_1r,risk.1[,i]))
    S2t = sum(multiply(THETA_1d,risk.1[,i])+multiply(THETA_1r2d,risk.1r2d[,i])+multiply(THETA_1nr2d,risk.1nr2d[,i]))
    
    sum1=ifelse(S1t==0,0,sum(ifelse(event.1r[,i]==0, 0, log.(THETA_1r/S1t)*event.1r[,i])))
    sum2=ifelse(S2t==0,0,sum(ifelse(event.1d[,i]==0, 0, log.(THETA_1d/S2t)*event.1d[,i]))+
                  sum(ifelse(event.1r2d[,i]==0, 0, log.(THETA_1r2d/S2t)*event.1r2d[,i]))+
                  sum(ifelse(event.1nr2d[,i]==0, 0, log.(THETA_1nr2d/S2t)*event.1nr2d[,i])))
    
    logl=logl+sum1+sum2
  }
  return(-logl) #-logl --> minimize --> optim
}