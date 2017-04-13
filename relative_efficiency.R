##############################################################
# Compare Variance of S(t) from joint model vs. WRSE vs. IPWE
# Last modified: April 13, 2017
############################################################## 
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
seed=1234; n=200; tolerance=1e-5; maxit=50; pi.x <- 0.5; pi.z1.R <- .5; pi.z1.NR <- .5; pi.z2.R=.5; pi.z2.NR <- .5
Tcheck=5; adt=c(3,20); decimal=1; bshape1=3/4; bscale1=4; bshape2=3; bscale2=4
tbetas <- c(0.6, -1.3, -0.7, -1, 1.8)#eta;
thetaZ=c('X'); gammaZ=c('X'); etaZ=c('X','Zr','X1Zr1'); muZ=NULL
np.t=length(thetaZ); np.g=length(gammaZ); np.e=length(etaZ); np.m=length(muZ); n.p=np.t+np.g+np.e+np.m 
TimeRDu = round(seq(0, adt[2], by=10^(-decimal)),decimal)

########################################################## Get bootstrapped SE estimate of S(t) ####
savefit = guo_tsiatis = lunceford = list(NULL)
for (bootiter in 1:bootstrapN){
  if( bootiter %% 1 == 0 ) cat(paste('boostrap',bootiter),'\n')
  simudata <- simSMART(seed=seed+bootiter, n, tbetas, pi.x, pi.z1.R, pi.z1.NR, pi.z2.R, pi.z2.NR, Tcheck, bshape1, bscale1, bshape2, bscale2, adt, decimal)
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

  ##Guo & Tsiatis (2005) Weighted-risk set estimator:
  simudata$TR = 0; simudata$TR[simudata$R==1]=simudata$T1[simudata$R==1]
  simudata$Z=simudata$Zr; simudata$Z[is.na(simudata$Z)]=0
  simudata$U=simudata$T2
  simudata$delta=simudata$C2
  simudata2=subset(simudata, select=c(X,Z,TR,R,U,delta))
  est01 = WRSEestimate(data=simudata2); est01
  guo_tsiatis[[bootiter]] = cbind(est01$time, est01$SURV11, est01$SURV12, est01$SURV21, est01$SURV22, est01$SE11, est01$SE12, est01$SE21, est01$SE22, est01$COV1112, est01$COV2122)
  colnames(guo_tsiatis[[bootiter]]) = c('time','SURV11', 'SURV12', 'SURV21', 'SURV22', 'SE11', 'SE12', 'SE21', 'SE22', 'COV1112', 'COV2122')
  ##Lunceford et al (2002)
  est01LDT = LDTestimate(simudata2)
  lunceford[[bootiter]] = cbind(est01LDT$time, est01LDT$SURV11, est01LDT$SURV12, est01LDT$SURV21, est01LDT$SURV22, est01LDT$SE11, est01LDT$SE12, est01LDT$SE21, est01LDT$SE22, est01LDT$COV1112, est01LDT$COV2122)
  colnames(lunceford[[bootiter]]) = c('time','SURV11', 'SURV12', 'SURV21', 'SURV22', 'SE11', 'SE12', 'SE21', 'SE22', 'COV1112', 'COV2122')
}


########################################################## RELATIVE EFFICIENCY ####
### Bootstrapped SE from joint model
timeseq= c(1,2,3,4,5,6,7)
TimeRDu =  savefit[[1]][,1]; len.et=length(TimeRDu)
seqt = which(TimeRDu %in% timeseq)
sigma =matrix(NA,nrow=length(seqt), ncol=4)
for (t in 1:length(seqt)){
  St=do.call("rbind",lapply(savefit, function(x) x[seqt[t],]))
  StV0 = St[,2:5]
  sigma[t,] = apply(StV0,2,sd)
}
joint.sigma=cbind(timeseq,sigma)

### Bootstrapped SE from WRSE
gt.sigma =matrix(NA,nrow=length(seqt), ncol=4)
for (t in 1:length(seqt)){
  St=do.call("rbind",lapply(guo_tsiatis, function(x) x[x[,1]==timeseq[t]]))
  StV0 = St[,2:5]
  gt.sigma[t,] = apply(StV0,2,sd)
}
gt.sigma=cbind(timeseq,gt.sigma)

### Bootstrapped SE from IPWE
lf.sigma=matrix(NA,nrow=length(seqt), ncol=4)
for (t in 1:length(seqt)){
  St=do.call("rbind",lapply(lunceford, function(x) x[x[,1]==timeseq[t]]))
  StV0 = St[,2:5]
  lf.sigma[t,] = apply(StV0,2,sd)
}
lf.sigma=cbind(timeseq,lf.sigma)

#### Relative efficiency: 
relative.eff.gt=relative.eff.lf = matrix(NA,nrow=nrow(gt.sigma), ncol=ncol(gt.sigma))
for (i in 1:nrow(relative.eff.gt)){
  #Var(Joint model) / Var (WRSE):
  relative.eff.gt[i,]=c(gt.sigma[i,1], (joint.sigma[i,-1]^2)/(gt.sigma[i,-1]^2))
  #Var(Joint model) / Var (IPWE):
  relative.eff.lf[i,]=c(lf.sigma[i,1], (joint.sigma[i,-1]^2)/(lf.sigma[i,-1]^2))
}
round(relative.eff.gt,2)
round(relative.eff.lf,2)

#### Instead of using bootstrapped SE, WRSE and IPWE also have estimator for SE 
library(plyr)
St.gt=data.frame(do.call("rbind",lapply(guo_tsiatis, function(x) {x[x[,1] %in% timeseq,]})))
gt.plyr <- ddply( St.gt, .(time), function(x) colMeans(x) )
gt.sigma=gt.plyr[,c('time','SE11','SE12','SE21','SE22')]
St.lf=data.frame(do.call("rbind",lapply(lunceford, function(x) {x[x[,1] %in% timeseq,]})))
lf.plyr <- ddply( St.lf, .(time), function(x) colMeans(x) )
lf.sigma=lf.plyr[,c('time','SE11','SE12','SE21','SE22')]

relative.eff.gt=relative.eff.lf = matrix(NA,nrow=nrow(gt.plyr), ncol=ncol(gt.sigma))
for (i in 1:nrow(relative.eff.gt)){
  #Var(Joint model) / Var (WRSE):
  relative.eff.gt[i,]=c(gt.sigma[i,1], (joint.sigma[i,-1]^2)/(gt.sigma[i,-1]^2))
  #Var(Joint model) / Var (WRSE):
  relative.eff.lf[i,]=c(lf.sigma[i,1], (joint.sigma[i,-1]^2)/(lf.sigma[i,-1]^2))
}
round(relative.eff.gt,2)
round(relative.eff.lf,2)