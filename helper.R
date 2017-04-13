## Generic helper functions
pars.cal.t <-function(betas, tz, gz, ez, mz){
  b.theta <- matrix(betas[1:(np.t)],ncol=1)
  b.gamma <- matrix(betas[(np.t+1):(np.t+np.g)],ncol=1)
  b.eta <- matrix(betas[(np.t+np.g+1):(np.t+np.g+np.e)],ncol=1)
  b.mu <- matrix(betas[(np.t+np.g+np.e+1):(np.t+np.g+np.e+np.m)],ncol=1)
  theta <- exp(tz%*%b.theta) # theta
  gamma <- exp(gz%*%b.gamma)   # gamma
  eta <- gamma*exp(ez%*%b.eta)   # gamma
  if(ncol(mz)==0) {mu <- gamma}else{mu=gamma*exp(mz%*%b.mu)}  # gamma
  pars <- cbind(theta,gamma,eta,mu)
  return(pars)
}

multiply=function(a,b){ifelse(a==0|b==0,0,a*b)}
divide = function(a,b){ifelse(b==0, 0, a/b)}
log. <- function(x) {b<-ifelse(x>0,log(x),1e-10); return(b)}

hex=c("#88CCEE","#ff7f00","#4daf4a","#e41a1c","#984ea3","#ffff33","#a65628","#f781bf", "#BEBEBE", "#666666")

## Survival prediction
s.km2.pred <- function(t, theta, gamma, eta, mu){
  eta=gamma*eta; mu_g=gamma*mu
  cumH1t = cumsum(H.fit.1)[TimeRDu==t]
  cumH2t = cumsum(H.fit.2)[TimeRDu==t]
  cumH1tnr = cumsum(H.fit.1)[TimeRDu==Tcheck]
  cumH2tnr = cumsum(H.fit.2)[TimeRDu==Tcheck]
  if(t<=Tcheck){
    temp1=exp(-(cumH1t*theta+cumH2t*gamma))
    trRange = TimeRDu[TimeRDu<=t]  
    dtr=diff(trRange)
    integaltr=0
    if(length(trRange)>1){
      for (i in 2:length(trRange)){
        tr = trRange[i]
        h1r=H.fit.1[TimeRDu==tr]
        cumH1r=cumsum(H.fit.1)[TimeRDu==tr]
        cumH2r=cumsum(H.fit.2)[TimeRDu==tr]
        temp2 = theta*h1r*exp(-(cumH1r*theta+cumH2r*(gamma-eta)+cumH2t*eta))
        integaltr = integaltr + temp2#*dtr[i-1]
      }
    }
    surv = temp1+integaltr
  }
  if(t>Tcheck){
    temp1=exp(-(cumH1tnr*theta+cumH2tnr*(gamma-mu_g)+cumH2t*mu_g))
    trRange = TimeRDu[TimeRDu<=Tcheck]  
    dtr=diff(trRange)
    integaltr=0
    if(length(trRange)>1){
      for (i in 2:length(trRange)){
        tr = trRange[i]
        h1r=H.fit.1[TimeRDu==tr]
        cumH1r=cumsum(H.fit.1)[TimeRDu==tr]
        cumH2r=cumsum(H.fit.2)[TimeRDu==tr]
        temp2 = theta*h1r*exp(-(cumH1r*theta+cumH2r*(gamma-eta)+cumH2t*eta))
        integaltr = integaltr + temp2#*dtr[i-1]
      }
    }
    surv = temp1+integaltr
  }
  return(surv)
}

## Survival estimator (Tang & Wahed, 2015), derived from cummulative hazard from package 'DTR':
CHRestimate2 = function (data, covar = names(data)[!names(data) %in% c("X", "R", "Z", "U", "delta")]){
  n = nrow(data)
  X = data$X #1st treatment
  R = data$R #response status
  Z = data$Z #2nd treatment
  U = data$U #survival/censor time
  delta = data$delta #censor indicator
  t = unique(U[which(delta == 1)])
  t = t[order(t)]
  n.risk = apply(as.array(t), 1, function(x) sum(as.numeric(U >= x)))
  n.event = apply(as.array(t), 1, function(x) length(which(U == x & delta == 1)))
  pi.x = sum(X)/n
  pi.z1 = sum((1 - X) * R * Z)/sum((1 - X) * R)
  pi.z2 = sum(X * R * Z)/sum(X * R)
  w11 = (1 - X) * (1 - R)/(1 - pi.x) + (1 - X) * R * (1 - Z)/((1 - pi.x) * (1 - pi.z1))
  w12 = (1 - X) * (1 - R)/(1 - pi.x) + (1 - X) * R * Z/((1 - pi.x) * pi.z1)
  w21 = X * (1 - R)/pi.x + X * R * (1 - Z)/(pi.x * (1 - pi.z1))
  w22 = X * (1 - R)/pi.x + X * R * Z/(pi.x * pi.z1)
  if (length(covar) == 0) {
    stop("Covariate(s) can not be empty")
  }else {
    if (FALSE %in% (covar %in% names(data))) {
      stop("Covariate(s) can not be found in the data")
    }else {
      V = as.matrix(data[, names(data) %in% covar])
    }
  }
  est = lest = NULL
  if (NCOL(V) == 1) {
    beta = as.numeric(coxph(Surv(U, delta) ~ ., data = data[, names(data) %in% c("U", "delta", covar)])$coef)
    cat("Calling for updateBeta() function to solve for coefficients... \n")
    for (p in 1:1000) {
      ebeta = updateBeta(beta, V, U, delta, w11, w12, w21, w22)
      index = max(abs(ebeta - beta))
      if (index <= 10^(-6)) 
        break
      else {
        beta = ebeta
        p = p + 1
      }
    }
    cat("Calculating survival rates... \n")
    V = as.numeric(V)
    e = exp(ebeta * V)
    for (j in 1:length(t)) { 
      s00 = s01 = s02 = s03 = rep(0, n)
      v10_bar = v11_bar = v12_bar = v13_bar = rep(0,n)
      lambda11 = lambda12 = lambda21 = lambda22 = 0
      omega = 0
      h11 = h12 = h21 = h22 = 0
      for (i in 1:n) {
        ind = as.numeric(U >= U[i])
        s00[i] = sum(ind * w11 * e)/n
        s01[i] = sum(ind * w12 * e)/n
        s02[i] = sum(ind * w21 * e)/n
        s03[i] = sum(ind * w22 * e)/n
        s10 = sum(ind * w11 * e * V)/n
        s11 = sum(ind * w12 * e * V)/n
        s12 = sum(ind * w21 * e * V)/n
        s13 = sum(ind * w22 * e * V)/n
        s20 = sum(ind * w11 * e * V * V)/n
        s21 = sum(ind * w12 * e * V * V)/n
        s22 = sum(ind * w21 * e * V * V)/n
        s23 = sum(ind * w22 * e * V * V)/n
        if (s00[i] != 0) 
          v10_bar[i] = s10/s00[i]
        else v10_bar[i] = 0
        if (s01[i] != 0) 
          v11_bar[i] = s11/s01[i]
        else v11_bar[i] = 0
        if (s02[i] != 0) 
          v12_bar[i] = s12/s02[i]
        else v12_bar[i] = 0
        if (s03[i] != 0) 
          v13_bar[i] = s13/s03[i]
        else v13_bar[i] = 0
        if (s00[i] != 0) 
          v20_bar = s20/s00[i]
        else v20_bar = 0
        if (s01[i] != 0) 
          v21_bar = s21/s01[i]
        else v21_bar = 0
        if (s02[i] != 0) 
          v22_bar = s22/s02[i]
        else v22_bar = 0
        if (s03[i] != 0) 
          v23_bar = s23/s03[i]
        else v23_bar = 0
        if (s00[i] != 0) lambda11 = lambda11 + delta[i] * w11[i] * as.numeric(U[i] <= t[j])/(n * s00[i])
        if (s01[i] != 0) lambda12 = lambda12 + delta[i] * w12[i] * as.numeric(U[i] <= t[j])/(n * s01[i])
        if (s02[i] != 0) lambda21 = lambda21 + delta[i] * w21[i] * as.numeric(U[i] <= t[j])/(n * s02[i])
        if (s03[i] != 0) lambda22 = lambda22 + delta[i] * w22[i] * as.numeric(U[i] <= t[j])/(n * s03[i])
        if (s00[i] != 0) h11 = h11 - delta[i] * w11[i] * as.numeric(U[i] <= t[j]) * v10_bar[i]/(n * s00[i])
        if (s01[i] != 0) h12 = h12 - delta[i] * w12[i] * as.numeric(U[i] <= t[j]) * v11_bar[i]/(n * s01[i])
        if (s02[i] != 0) h21 = h21 - delta[i] * w21[i] * as.numeric(U[i] <= t[j]) * v12_bar[i]/(n * s02[i])
        if (s03[i] != 0) h22 = h22 - delta[i] * w22[i] * as.numeric(U[i] <= t[j]) * v13_bar[i]/(n * s03[i])
        tao11i = v20_bar - v10_bar[i] * v10_bar[i]
        tao12i = v21_bar - v11_bar[i] * v11_bar[i]
        tao21i = v22_bar - v12_bar[i] * v12_bar[i]
        tao22i = v23_bar - v13_bar[i] * v13_bar[i]
        omega = omega + (delta[i] * w11[i] * tao11i + 
                           delta[i] * w12[i] * tao12i + 
                           delta[i] * w21[i] * tao21i + 
                           delta[i] * w22[i] * tao22i)/n
      }
      temp = c(exp(-lambda11),exp(-lambda12),exp(-lambda21),exp(-lambda22))
      #temp[which(temp == 0 | is.na(temp) == TRUE | temp == Inf | temp == -Inf)] = NA
      est = rbind(est, temp)
      rownames(est) = NULL
    }
  }
  results = list(Call = match.call(), coefficients = ebeta, t=t, S11 = est[, 1], S12 = est[, 2], S21 = est[,3], S22 = est[, 4])
  return(results)
}