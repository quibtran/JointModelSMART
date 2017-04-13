##################################################################
# Joint modeling of time-to-response and time-to-death in SMART 
# Data: Simulated data from function simSMART()
# Estimators: 
#   - Beta coefficients from NPMLE using optim()
#   - SE from hessian matrix
# Survival diagnostic plots:
#   - Time-to-response: Kaplan-Meier curve vs. joint model prediction 
#   - Time-to-death: Joint model vs. WRSE (Guo&Tsiatis,2005) vs. 
#                     IPWE (Lunceford et al, 2002) vs. Tang Wahed (2015)
# Last modified: April 13, 2017
##################################################################
rm(list = ls())
library(survival)
library(numDeriv)
library(plyr)

source('helper.R') # Misc. helper functions
source('simSMART.R') # simulate SMART with survival outcomes
source('calH.R') # Estimate 2 non-parametric baseline hazards for response and death
source('loglik.R') # Estimate log-likelihood

########################################################## Set parameters ####
seed=1234; n=500; tolerance=1e-5; maxit=50; pi.x <- 0.5; pi.z1.R <- .5; pi.z1.NR <- .5; pi.z2.R=.5; pi.z2.NR <- .5
Tcheck=5; adt=c(3,20); decimal=1; bshape1=3/4; bscale1=4; bshape2=3; bscale2=4
tbetas <- c(-1.5, 0.6, -0.5, -1.3, -1.2, -0.7, -1, 1.8, -2.4)
thetaZ=c('V','X'); gammaZ=c('V','X'); etaZ=c('1','X','Zr','X1Zr1','X1Zr1V1'); muZ=NULL
np.t=length(thetaZ); np.g=length(gammaZ); np.e=length(etaZ); np.m=length(muZ); n.p=np.t+np.g+np.e+np.m 
TimeRDu = round(seq(0, adt[2], by=10^(-decimal)),decimal)

########################################################## Simulate data ####
simudata <- simSMART(seed, n, tbetas, pi.x, pi.z1.R, pi.z1.NR, pi.z2.R, pi.z2.NR, 
                     Tcheck, bshape1, bscale1, bshape2, bscale2, adt, decimal)
summary(simudata)
Tind <- order(simudata$T1,simudata$T2)
sorted <- simudata[Tind,]
TimeRDu = sort(unique(c(simudata$T1,simudata$T2)))
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

########################################################## RUN MLE (may take a couple minutes for n=500)####
test <- optim(par=rep(0.5,n.p),fn=loglik,gr=NULL, method = "L-BFGS-B", 
              lower=rep(-10,n.p), upper=rep(10,n.p), hessian=T, control=list(maxit=1000,trace=1))
hessian.mat=test$hessian
p.se <- sqrt(abs(diag(solve(-hessian.mat))))
para.b <- test$p
par.var.b <- (p.se)^2 
para.true.b <- rep(0,length(para.b)) 
p.value.b <- 2*pnorm(abs(para.b - para.true.b)/sqrt(abs(par.var.b)),lower.tail=F) 
result <- round(data.frame(para.b,sqrt(abs(par.var.b)),p.value.b),3)
names(result) <- c('est','se','pvalue')
result

########################################################## DIAGNOSTIC PLOTS ####
hex=c("#88CCEE","#ff7f00","#4daf4a","#e41a1c","#984ea3","#ffff33","#a65628","#f781bf", "#BEBEBE", "#666666")
H.fit=calH(para.b)
H.fit.1=H.fit$H1 
H.fit.2=H.fit$H2 
Zt.V=para.b[1]; Zt.X=para.b[2]
Zg.V=para.b[3]; Zg.X=para.b[4]
Ze.1=para.b[5]; Ze.X=para.b[6]; Ze.Zr=para.b[7]; Ze.X1Zr1=para.b[8]; Ze.X1Zr1V1=para.b[9]
Zm.1=Zm.X=Zm.Zr=0
########################### Time-to-response: Kaplan-Meier vs joint model  ####
#Kaplan-Meier:
dN1.V1X0=dN1.V1X1=dN1.V0X0=dN1.V0X1=array(0,len.et)
N1.risk.V1X0=N1.risk.V1X1=N1.risk.V0X0=N1.risk.V0X1=array(0,len.et)
for (t in 1:len.et){
  dN1.V1X0[t]=sum((sorted$T1==TimeRDu[t] & sorted$R==1 & sorted$V==1 & sorted$X==0))
  dN1.V1X1[t]=sum((sorted$T1==TimeRDu[t] & sorted$R==1 & sorted$V==1 & sorted$X==1))
  dN1.V0X0[t]=sum((sorted$T1==TimeRDu[t] & sorted$R==1 & sorted$V==0 & sorted$X==0))
  dN1.V0X1[t]=sum((sorted$T1==TimeRDu[t] & sorted$R==1 & sorted$V==0 & sorted$X==1))
  N1.risk.V1X0[t]=sum((sorted$T1>=TimeRDu[t] & sorted$V==1 & sorted$X==0))
  N1.risk.V1X1[t]=sum((sorted$T1>=TimeRDu[t] & sorted$V==1 & sorted$X==1))
  N1.risk.V0X0[t]=sum((sorted$T1>=TimeRDu[t] & sorted$V==0 & sorted$X==0))
  N1.risk.V0X1[t]=sum((sorted$T1>=TimeRDu[t] & sorted$V==0 & sorted$X==1))
}
S.KM.1.V1X0 <- cumprod(1-(ifelse(N1.risk.V1X0==0,0,dN1.V1X0/N1.risk.V1X0)))
S.KM.1.V1X1 <- cumprod(1-(ifelse(N1.risk.V1X1==0,0,dN1.V1X1/N1.risk.V1X1)))
S.KM.1.V0X0 <- cumprod(1-(ifelse(N1.risk.V0X0==0,0,dN1.V0X0/N1.risk.V0X0)))
S.KM.1.V0X1 <- cumprod(1-(ifelse(N1.risk.V0X1==0,0,dN1.V0X1/N1.risk.V0X1)))
#Joint model:
s.km1.pred <- function(t,h1,h2,theta){
  t=ifelse(t<=Tcheck,t,Tcheck)
  cumH1t = cumsum(h1)[TimeRDu==t]
  cumH2t = cumsum(h2)[TimeRDu==t]
  exp(-cumH1t*theta)
}
S.KM.1.V1X0.hat <- sapply(TimeRDu,FUN=s.km1.pred, h1=H.fit.1,h2=H.fit.2, theta=exp(Zt.V))
S.KM.1.V1X1.hat <- sapply(TimeRDu,FUN=s.km1.pred, h1=H.fit.1,h2=H.fit.2, theta=exp(Zt.V+Zt.X))
S.KM.1.V0X0.hat <- sapply(TimeRDu,FUN=s.km1.pred, h1=H.fit.1,h2=H.fit.2, theta=exp(0))
S.KM.1.V0X1.hat <- sapply(TimeRDu,FUN=s.km1.pred, h1=H.fit.1,h2=H.fit.2, theta=exp(0+Zt.X))
#PLOT:
plot(TimeRDu, S.KM.1.V1X0, col=hex[1], ylim=range(S.KM.1.V1X0.hat,S.KM.1.V1X1.hat,S.KM.1.V0X0.hat,S.KM.1.V0X1.hat), xlim=c(0,Tcheck+2),
     main="Survival curve for event of response \n stratified by baseline covariates and first treatment",
     xlab="Age", ylab="Fraction Surviving",  lwd=2, type="s")
lines(TimeRDu, S.KM.1.V1X1, col=hex[2], lty=1, lwd=2, type="s")
lines(TimeRDu, S.KM.1.V0X0, col=hex[3], lty=1, lwd=2, type="s")
lines(TimeRDu, S.KM.1.V0X1, col=hex[4], lty=1, lwd=2, type="s")
lines(TimeRDu, S.KM.1.V1X0.hat, col=hex[1], lty=2, type='l',lwd=3)
lines(TimeRDu, S.KM.1.V1X1.hat, col=hex[2], lty=2, type='l',lwd=3)
lines(TimeRDu, S.KM.1.V0X0.hat, col=hex[3], lty=2, type='l',lwd=3)
lines(TimeRDu, S.KM.1.V0X1.hat, col=hex[4], lty=2, type='l',lwd=3)
legend('bottomleft',col=hex[c(1,2,3,4)], lwd=3, legend=c("V=1, A1", "V=1, A2", "V=0, A1","V=0, A2"))

########################### Time-to-death: Joint model vs. WRSE (Guo&Tsiatis,2005) vs. IPWE (Lunceford et al, 2002) vs. Tang Wahed (2015)  ####
## Joint model 
fit.V0A1B1=sapply(TimeRDu, FUN=s.km2.pred, theta=exp(0), gamma=exp(0), eta=exp(Ze.1), mu=exp(Zm.1))
fit.V0A1B2=sapply(TimeRDu,FUN=s.km2.pred, theta=exp(0), gamma=exp(0), eta=exp(Ze.1+Ze.Zr), mu=exp(Zm.1))
fit.V0A2B1=sapply(TimeRDu,FUN=s.km2.pred, theta=exp(0+Zt.X), gamma=exp(0+Zg.X), eta=exp(Ze.1+Ze.X), mu=exp(Zm.1+Zm.X))
fit.V0A2B2=sapply(TimeRDu,FUN=s.km2.pred, theta=exp(0+Zt.X), gamma=exp(0+Zg.X), eta=exp(Ze.1+Ze.X+Ze.Zr+Ze.X1Zr1), mu=exp(Zm.1+Zm.X))
fit.V1A1B1=sapply(TimeRDu, FUN=s.km2.pred, theta=exp(Zt.V), gamma=exp(Zg.V), eta=exp(Ze.1), mu=exp(Zm.1))
fit.V1A1B2=sapply(TimeRDu,FUN=s.km2.pred, theta=exp(Zt.V), gamma=exp(Zg.V), eta=exp(Ze.1+Ze.Zr), mu=exp(Zm.1))
fit.V1A2B1=sapply(TimeRDu,FUN=s.km2.pred, theta=exp(Zt.V+Zt.X), gamma=exp(Zg.V+Zg.X), eta=exp(Ze.1+Ze.X), mu=exp(Zm.1+Zm.X))
fit.V1A2B2=sapply(TimeRDu,FUN=s.km2.pred, theta=exp(Zt.V+Zt.X), gamma=exp(Zg.V+Zg.X), eta=exp(Ze.1+Ze.X+Ze.Zr+Ze.X1Zr1+Ze.X1Zr1V1), mu=exp(Zm.1+Zm.X))
savefit = cbind(TimeRDu, fit.V0A1B1, fit.V0A1B2, fit.V0A2B1, fit.V0A2B2,fit.V1A1B1, fit.V1A1B2, fit.V1A2B1, fit.V1A2B2)

## Survival estimator (Guo&Tsiatis) from package 'DTR':
library(DTR)
simudata$TR = 0; simudata$TR[simudata$R==1]=simudata$T1[simudata$R==1]
simudata$Z=simudata$Zr; simudata$Z[is.na(simudata$Z)]=0
simudata$U=simudata$T2
simudata$delta=simudata$C2
simudata2=subset(simudata, select=c(V,X,Z,TR,R,U,delta))
est01 = WRSEestimate(data=simudata2); est01
est0 = WRSEestimate(data=simudata2[simudata2$V==0,]); est0
est1 = WRSEestimate(data=simudata2[simudata2$V==1,]); est1
guo_tsiatis0 = cbind(est0$time, est0$SURV11, est0$SURV12, est0$SURV21, est0$SURV22, est0$SE11, est0$SE12, est0$SE21, est0$SE22, est0$COV1112, est0$COV2122)
guo_tsiatis1 = cbind(est1$time, est1$SURV11, est1$SURV12, est1$SURV21, est1$SURV22, est1$SE11, est1$SE12, est1$SE21, est1$SE22, est1$COV1112, est1$COV2122)
colnames(guo_tsiatis0) = colnames(guo_tsiatis1) = c('time','SURV11', 'SURV12', 'SURV21', 'SURV22', 'SE11', 'SE12', 'SE21', 'SE22', 'COV1112', 'COV2122')

## Survival estimator (Lunceford) from package 'DTR':
est01LDT = LDTestimate(simudata2)
est0LDT = LDTestimate(simudata2[simudata2$V==0,])
est1LDT = LDTestimate(simudata2[simudata2$V==1,])
lunceford0 = cbind(est0LDT$time, est0LDT$SURV11, est0LDT$SURV12, est0LDT$SURV21, est0LDT$SURV22, est0LDT$SE11, est0LDT$SE12, est0LDT$SE21, est0LDT$SE22, est0LDT$COV1112, est0LDT$COV2122)
lunceford1 = cbind(est1LDT$time, est1LDT$SURV11, est1LDT$SURV12, est1LDT$SURV21, est1LDT$SURV22, est1LDT$SE11, est1LDT$SE12, est1LDT$SE21, est1LDT$SE22, est1LDT$COV1112, est1LDT$COV2122)
colnames(lunceford0) = c('time','SURV11', 'SURV12', 'SURV21', 'SURV22', 'SE11', 'SE12', 'SE21', 'SE22', 'COV1112', 'COV2122')
colnames(lunceford1) = c('time','SURV11', 'SURV12', 'SURV21', 'SURV22', 'SE11', 'SE12', 'SE21', 'SE22', 'COV1112', 'COV2122')

## Survival estimator (Tang & Wahed, 2015), derived from cummulative hazard from package 'DTR':
est0TW =  CHRestimate2(data=simudata2, covar=c("V")); est0TW
est1TW = est0TW
est1TW$S11=est0TW$S11^(exp(est0TW$coefficients))
est1TW$S12=est0TW$S12^(exp(est0TW$coefficients))
est1TW$S21=est0TW$S21^(exp(est0TW$coefficients))
est1TW$S22=est0TW$S22^(exp(est0TW$coefficients))

#PLOT:
fit=savefit[,2:5]; est=est0; estLDT=est0LDT; estTW=est0TW
#fit=savefit[,6:9]; est=est1; estLDT=est1LDT; estTW = est1TW
plot(estLDT$time, estLDT$SURV11, type='l', lwd=2, lty=2, col=hex[10],
     main="Predicted survival rates for group with baseline V=0",
     ylim=c(0,1), xlim=c(0,max(TimeRDu)), ylab='Fraction Surviving', xlab='t')
lines(estLDT$time, estLDT$SURV12, type='l',col=hex[10], lwd=2, lty=2)
lines(estLDT$time, estLDT$SURV21, type='l',col=hex[10], lwd=2, lty=2)
lines(estLDT$time, estLDT$SURV22, type='l',col=hex[10], lwd=2, lty=2)
lines(est$time, est$SURV11, type='l', col=hex[10],lwd=2, lty=3)
lines(est$time, est$SURV12, type='l', col=hex[10],lwd=2, lty=3)
lines(est$time, est$SURV21, type='l', col=hex[10],lwd=2, lty=3)
lines(est$time, est$SURV22, type='l', col=hex[10],lwd=2, lty=3)
lines(estTW$t, estTW$S11, type='l', col=hex[1],lwd=2, lty=5)
lines(estTW$t, estTW$S12, type='l', col=hex[2],lwd=2, lty=5)
lines(estTW$t, estTW$S21, type='l', col=hex[3],lwd=2, lty=5)
lines(estTW$t, estTW$S22, type='l', col=hex[4],lwd=2, lty=5)
lines(TimeRDu, fit[,1], type='l',col=hex[1],  lwd=3)  
lines(TimeRDu, fit[,2], type='l',col=hex[2],  lwd=3)  
lines(TimeRDu, fit[,3], type='l',col=hex[3],  lwd=3)  
lines(TimeRDu, fit[,4], type='l',col=hex[4],  lwd=3)  
legend('topright', col=c(hex[c(10,10,10,1:4)]), lwd=2,lty=c(2,3,5,1,1,1,1), legend=c("IPWE", "WRSE", "Tang & Wahed (2013)","Joint model: A1B1", "Joint model: A1B2", "Joint model: A2B1","Joint model: A2B2"))
