# Function to simulate SMART design II with survival outcomes
simSMART = function (seed, n, tbetas, pi.x, pi.z1.R, pi.z1.NR, pi.z2.R, pi.z2.NR, 
                     Tcheck, bshape1, bscale1, bshape2, bscale2, adt, decimal) 
{
  set.seed(seed)
  b.theta <- matrix(tbetas[(1):(np.t)],ncol=1)
  b.gamma <- matrix(tbetas[(np.t+1):(np.t+np.g)],ncol=1)
  b.eta <- matrix(tbetas[(np.t+np.g+1):(np.t+np.g+np.e)],ncol=1)
  b.mu <- matrix(tbetas[(np.t+np.g+np.e+1):(np.t+np.g+np.e+np.m)],ncol=1)
  V <- rbinom(n, 1, 0.75)#rnorm(n,2,1)
  X <- rbinom(n, 1, pi.x)
  dat=cbind('1'=1, V,X)
  thetaZs = as.matrix(subset(dat,select=thetaZ))
  gammaZs = as.matrix(subset(dat,select=gammaZ))
  
  #Simulate Tr, time to response
  theta=exp(thetaZs%*%b.theta) 
  u <- runif(n)
  Tr <- (-log(u)/(theta))^(1/bshape1)*bscale1 #time diagnosis
  
  #Death time before response or check-up time
  gamma <- exp(gammaZs%*%b.gamma) 
  v <- runif(n)
  Td <- (-log(v)/(gamma))^(1/bshape2)*bscale2
  
  #Response status
  minT = pmin(Tr, Td, Tcheck)
  R = 1*(Tr==minT)
  NR = 1*(Tcheck==minT)
  D = 1*(Td==minT)
  
  #Second treatment 
  Zr <- Znr <- rep(NA, n) #NA for subjects that die before response OR no-response by time of check-up
  Zr[which(X==0 & R == 1)] <- rbinom(length(which(X==0 & R == 1)), 1, pi.z1.R)
  Zr[which(X==1 & R == 1)] <- rbinom(length(which(X==1 & R == 1)), 1, pi.z2.R)
  Znr[which(X==0 & NR == 1)] <- rbinom(length(which(X==0 & NR == 1)), 1, pi.z1.NR)
  Znr[which(X==1 & NR == 1)] <- rbinom(length(which(X==1 & NR == 1)), 1, pi.z2.NR)
  #Z[which(R==0)]=0
  X1Zr1=1*(Zr==1 & X==1)
  X1Zr0=1*(Zr==0 & X==1)
  X1Zr1V1=1*(Zr==1 & X==1 & V==1)
  #XZnr=Znr*X
  dat = cbind(dat,Zr,Znr,X1Zr1,X1Zr0, X1Zr1V1)
  
  #Simulate new Td for subject survive after response/check-up time
  etaZs = as.matrix(subset(dat,select=etaZ))
  muZs = as.matrix(subset(dat,select=muZ))
  eta <- gamma*exp(etaZs%*%b.eta) 
  if(ncol(muZs)==0) {mu <- gamma}else{mu=gamma*exp(muZs%*%b.mu)}  # gamma
  Td.after.response <- ((-log(v)-(gamma-eta)*(Tr/bscale2)^bshape2)/(eta))^(1/bshape2)*bscale2
  #Tr[R==1] < Td.after.response)[R==1] #check to make sure responded subject's death time is after response time
  Td.after.noresponse <- ((-log(v)-(gamma-mu)*(Tcheck/bscale2)^bshape2)/(mu))^(1/bshape2)*bscale2
  #Tcheck < Td.after.noresponse[NR==1] #check to make sure non-responded subject's death time is after check-up time
  Td[which(R == 1)] <- Td.after.response[which(R == 1)]
  Td[which(NR == 1)] <- Td.after.noresponse[which(NR == 1)]
  #cbind(minT,Td,R,NR,D)
  
  #Simulate censoring time
  tc = runif(n,adt[1],adt[2])
  T1 = pmin(minT,tc)
  R = R*(minT <= tc)
  NR = NR*(minT <= tc)
  D = D*(minT <= tc)
  T2 = pmin(Td, tc)
  C2 = (Td <= tc)
  simudata = data.frame(cbind(dat,T1, R, NR, D, T2, C2))
  names(simudata)=c(colnames(dat),'T1','R','NR','D','T2','C2')
  row.names(simudata) <- NULL
  ceiling_dec <- function(x, level) round(x + 5*10^(-level-1), level)
  simudata$T1 <- ceiling_dec(simudata$T1, decimal)
  simudata$T2 <- ceiling_dec(simudata$T2, decimal)
  simudata
}
