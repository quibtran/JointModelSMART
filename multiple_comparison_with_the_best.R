########################################################## 
# Multiple comparison with the best 
# Modified from : Ertefaie et al (2015)
# Last modified: April 13, 2017
########################################################## 
rm(list = ls())
library(survival)
library(DTR)
library(MASS)
source('helper.R') # Misc. helper functions
source('simSMART.R') # simulate SMART with survival outcomes
source('calH.R') # Estimate 2 non-parametric baseline hazards for response and death
source('loglik.R') # Estimate log-likelihood

########################################################## Set parameters ####
bootstrapN = 5 # number of bootstrap replication; set to 500 in paper
seed=1234; n=200; tolerance=1e-5; maxit=50; pi.x <- 0.5; pi.z1.R <- .4; pi.z1.NR <- .6; pi.z2.R=.6; pi.z2.NR <- .7; Tcheck=5; adt=c(3,20); decimal=1; bshape1=3/4; bscale1=4; bshape2=3; bscale2=4
thetaZ=c('X'); gammaZ=c('X'); etaZ=c('X','Zr'); muZ=NULL
np.t=length(thetaZ); np.g=length(gammaZ); np.e=length(etaZ); np.m=length(muZ); n.p=np.t+np.g+np.e+np.m 
TimeRDu = round(seq(0, adt[2], by=10^(-decimal)),decimal)
tseq = c(7,15) # time at which we want to compare S(t) between regimens
delta=1.3 #results in paper show delta ranging between 0 and 2
tbetas=c(delta, -delta, -delta/4, -delta/2) 

########################################################## Get bootstrapped SE estimate of S(t) ####
savefit = list(NULL)
for (bootiter in 1:bootstrapN){
  if( bootiter %% 5 == 0 ) cat(paste('boostrap',bootiter),'\n')
  simudata <- simudata <- simSMART(seed+bootiter, n, tbetas, pi.x, pi.z1.R, pi.z1.NR, pi.z2.R, pi.z2.NR, 
                                   Tcheck, bshape1, bscale1, bshape2, bscale2, adt, decimal)
  Tind <- order(simudata$T1,simudata$T2)
  sorted <- simudata[Tind,]
  len.et=length(TimeRDu)
  T1match= match(sorted$T1,TimeRDu)
  Tnr = which(TimeRDu==Tcheck); if(length(Tnr)==0){Tnr=len.et}
  
  thetaZs <- as.matrix(subset(sorted,select=thetaZ)) # covariate set for theta
  gammaZs <- as.matrix(subset(sorted,select=gammaZ)) # covariate set for gamma
  etaZs <- as.matrix(subset(sorted,select=etaZ)) # covariate set for gamma
  muZs <- as.matrix(subset(sorted,select=muZ)) # covariate set for gamma
  
  ## Change to counting process
  dN1r=dN1d=dN1r2d=dN1nr2d=array(0,len.et)
  for (t in 1:len.et){
    dN1r[t]=sum((sorted$T1==TimeRDu[t] & sorted$R==1))
    dN1d[t]=sum((sorted$T1==TimeRDu[t] & sorted$D==1))  
    dN1r2d[t]=sum((sorted$T2==TimeRDu[t] & sorted$R==1 & sorted$C2==1))  
    dN1nr2d[t]=sum((sorted$T2==TimeRDu[t] & sorted$NR==1 & sorted$C2==1))  
  }
  dNdeath=dN1d+dN1r2d+dN1nr2d
  
  risk.1 <- risk.1r2d <- risk.1nr2d <-array(0,dim=c(n,len.et))
  event.1r <- event.1d <- event.1r2d<- event.1nr2d <-array(0,dim=c(n,len.et))
  for (i in 1:nrow(sorted)){
    subjecti = sorted[i,]
    risk.1[i,] <- (subjecti$T1>=TimeRDu)
    risk.1r2d[i,] <- (subjecti$T1<TimeRDu)*(subjecti$R==1)*(subjecti$T2>=TimeRDu) 
    risk.1nr2d[i,] <- (subjecti$T1<TimeRDu)*(subjecti$NR==1)*(subjecti$T2>=TimeRDu)
    event.1r[i,] <- (subjecti$T1==TimeRDu)*(subjecti$R==1)
    event.1d[i,] <- (subjecti$T1==TimeRDu)*(subjecti$D==1)
    event.1r2d[i,] <- (subjecti$T2==TimeRDu)*(subjecti$R==1)*(subjecti$C2==1)
    event.1nr2d[i,] <- (subjecti$T2==TimeRDu)*(subjecti$NR==1)*(subjecti$C2==1)
  }   
  
  ## RUN MLE (may take a couple minutes for n=500)####
  test <- optim(par=rep(0.5,n.p),fn=loglik,gr=NULL, method = "L-BFGS-B", 
                lower=rep(-10,n.p), upper=rep(10,n.p), hessian=T, control=list(maxit=1000,trace=1))
  para.b <- test$p
  hessian.mat=test$hessian
  p.se <- sqrt(abs(diag(solve(-hessian.mat))))
  
  ## SURVIVAL PREDICTIONS ####
  H.fit=calH(para.b)
  H.fit.1=H.fit$H1 
  H.fit.2=H.fit$H2 
  Zt.V=0; Zt.X=para.b[1]
  Zg.V=0; Zg.X=para.b[2]
  Ze.1=0; Ze.X=para.b[3]; Ze.Zr=para.b[4]; Ze.X1Zr1=0; Ze.X1Zr0=0; Ze.X1Zr1V1=0
  Zm.1=Zm.X=Zm.Znr=Zm.XZnr=0
  
  fit.A1B1=sapply(TimeRDu, FUN=s.km2.pred, theta=exp(Zt.V), gamma=exp(Zg.V), eta=exp(Ze.1), mu=exp(Zm.1))
  fit.A1B2=sapply(TimeRDu,FUN=s.km2.pred, theta=exp(Zt.V), gamma=exp(Zg.V), eta=exp(Ze.1+Ze.Zr), mu=exp(Zm.1))
  fit.A2B1=sapply(TimeRDu,FUN=s.km2.pred, theta=exp(Zt.V+Zt.X), gamma=exp(Zg.V+Zg.X), eta=exp(Ze.1+Ze.X+Ze.X1Zr0), mu=exp(Zm.1+Zm.X))
  fit.A2B2=sapply(TimeRDu,FUN=s.km2.pred, theta=exp(Zt.V+Zt.X), gamma=exp(Zg.V+Zg.X), eta=exp(Ze.1+Ze.X+Ze.Zr+Ze.X1Zr1+Ze.X1Zr1V1), mu=exp(Zm.1+Zm.X))
  savefit[[bootiter]] = cbind(TimeRDu, fit.A1B1, fit.A1B2, fit.A2B1, fit.A2B2)
}

########################################################## Multiple comparison with the best ####
getC<-function(V,alpha, nrep){
  x<-mvrnorm(nrep,rep(0,4),V)
  for(i in 1:4){
    temp<-matrix(0,nrep,4)
    for(j in 1:4){
      if(j!=i){temp[,j]<-(x[,j]-x[,i])/sqrt(V[i,i]+V[j,j]-2*V[i,j])}
    }	
    temp<-temp[,-i]
    temp<-apply(temp,1,sort)
    temp<-apply(temp,2,max)
    c[i]<-quantile(temp,1-alpha)
  }
  c
} 

seqt = which(TimeRDu %in% tseq)
for (t in 1:length(tseq)){
  St=do.call("rbind",lapply(savefit, function(x) x[seqt[t],]))[,-1]
  V=cov(St)
  theta = St[1,]
  c<-rep(0,4) 
  for(i in 1:10){c<-c+getC(V,alpha=0.05, nrep=500)} 
  c<-c/10
  S.pos<-rep(1,4)
  for(i in 1:4){
    for(j in 1:4){
      S.pos[i]= S.pos[i]*(theta[i]>=theta[j]-c[i]*sqrt(V[i,i]+V[j,j]-2*V[i,j]))
    }
  }
  if(t==1) cat(' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~','\n')
  cat(' At t = ',tseq[t],', best regimen(s) are', c('A1B1','A1B2','A2B1','A2B2')[S.pos==1],'\n')
  if(t==length(tseq)) cat(' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~','\n')
}


