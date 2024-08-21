library(dplyr)

expit = function(x) {return(exp(x)/(1+exp(x)))}

sim_trt_func = function(X, D) {
  mu0.true = -3 + X%*%c(-1,0.7,0.6,-0.6)
  muD.true = 1 - 0.4/10*D + 3*X%*%c(0.2,0.2,0.3,-0.1) + (X*matrix(rep(D,4),ncol=4))%*%c(-0.1,0,0.1,0) - 0.03/10*D^3
  return(list(muD=muD.true, mu0=mu0.true))
}

gen_sim = function(seed, n) {
  set.seed(seed)
  X = matrix(rnorm(n*4),ncol=4)
  A.mean = expit(cbind(1,X) %*% (2*c(-0.1, 0.05, 0.05, -0.05, 0.15))) #c(-0.3, 0.25, 0.25, -0.25, 0.45)) #c(-0.3, 0.15, 0.15, -0.15, 0.25))
  A = rbinom(n=n, size=1, prob=A.mean)
  D.mean = ifelse(A == 0, NA, cbind(1,X) %*% c(3,0.2, 0.25,-0.3,0.5)) #c(3,0.45,0.4,-0.45,0.9))
  D <- rnorm(n, mean=D.mean,sd=2)
  
  trt_func_data = sim_trt_func(X, D)
  mu0.true <- trt_func_data$mu0
  muD.true <- trt_func_data$muD
  alpha.i = rnorm(n, mean=)
  Y0 = rnorm(n, mean=10+0.4*X[,1]-X[,2]+0.4*X[,3]+0.3*X[,4], sd=0.3)
  deltaY = rnorm(n,mean=ifelse(A == 0, mu0.true, muD.true), sd=0.7) #sd=0.1)
  Y1 = Y0 + deltaY
  W = cbind( exp(X[,1]/2) , 10+X[,2]/(1+exp(X[,1])) , (X[,1]*X[,3]/25 + 0.6)^3 , (X[,2]+X[,4]+20)^2)
  W = apply(W, 2, function(x) (x-mean(x))/sd(x))
  ids = 1:n
  dfk = data.frame(ids,X,W,A,D,Y0,Y1,deltaY)
  head(dfk)
  names(dfk) = c("ID", paste0("X", 1:4), paste0("W", 1:4), "A", "D", "Y0", "Y1", "deltaY")
  dfk$D2 = dfk$D^2
  dfk$D3 = dfk$D^3
  dfk$boot_wts = 1
  D.min <- min(dfk[A==1,]$D)
  D.max <- max(dfk[A==1,]$D)
  D.vals <- seq(D.min, D.max, length.out=50) #sort(dfk[A==1,]$D)[round(seq(1,sum(A==1),length.out=100))]
  return(list(data=dfk, D.vals=D.vals))
}

get_true_params = function() {
  n=1000000
  sim_info = gen_sim(1, n)
  trt_data = sim_info$data[sim_info$data$A == 1,]
  X = as.matrix(trt_data[,paste0("X", 1:4)])
  D.quants = quantile(trt_data$D, probs=c(0.10, 0.90))
  D.vals = seq(D.quants[1], D.quants[2], length.out=50)
  thetaD = sapply(D.vals, function(d.val) {
    d = rep(d.val, nrow(trt_data))
    mean(sim_trt_func(X, d)$muD)})
  theta0 = mean(sim_trt_func(X, rep(0, nrow(trt_data)))$mu0)
  D.dens = approx(density(trt_data$D), xout=D.vals)$y
  return(list(thetaD=thetaD, theta0=theta0, trt_effect=(thetaD - theta0), D.vals=D.vals, D.dens=D.dens))
}

get_true_params_setD = function(D.vals) {
  n=1000000
  sim_info = gen_sim(1, n)
  trt_data = sim_info$data[sim_info$data$A == 1,]
  X = as.matrix(trt_data[,paste0("X", 1:4)])
  thetaD = sapply(D.vals, function(d.val) {
    d = rep(d.val, nrow(trt_data))
    mean(sim_trt_func(X, d)$muD)})
  theta0 = mean(sim_trt_func(X, rep(0, nrow(trt_data)))$mu0)
  D.dens = approx(density(trt_data$D), xout=D.vals)$y
  return(list(thetaD=thetaD, theta0=theta0, trt_effect=(thetaD - theta0), D.vals=D.vals, D.dens=D.dens))
}

get_true_params_attD = function() {
  n=1000000
  sim_info = gen_sim(1, n)
  trt_data = sim_info$data[sim_info$data$A == 1,]
  X = as.matrix(trt_data[,paste0("X", 1:4)])
  D.quants = quantile(trt_data$D, probs=c(0.10, 0.90))
  D.vals = seq(D.quants[1], D.quants[2], length.out=50)
  thetaD = approx_fxn(x=trt_data$D, y=trt_data$deltaY, newx=D.vals)
  theta0 = sapply(D.vals, function(d.val) {
    d = rep(d.val, nrow(trt_data))
    mean(sim_trt_func(X, d)$mu0)})
  D.dens = approx(density(trt_data$D), xout=D.vals)$y
  return(list(thetaD=thetaD, theta0=theta0, trt_effect=(thetaD - theta0), D.vals=D.vals, D.dens=D.dens))
}