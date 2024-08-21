library(rlist)
get_sand_var_info_dose = function(data, psD_info, orD_info, wts, pseudo_info, bw,
                                  muhat_mat, mhat_mat, gps_mat) {
  var_data = cbind(ID=data$ID, psD_wts=wts, pseudo_out=pseudo_info$pseudo_dr, 
                   gps=psD_info$gps, predD=psD_info$predD, predMuD=orD_info$predMuD)
  return(list(data=var_data, bw=bw, muhat_mat=muhat_mat, mhat_mat=mhat_mat, gps_mat=gps_mat,
              orD_form=orD_info$orD_form, psD_form=psD_info$psD_form, sigma=psD_info$sigma))
}

get_sand_var_info_ctl = function(data, ps0_info, or0_info, wts, est) {
  var_data = cbind(data[,c("ID", "deltaY", "A", "D", "D2", "D3", paste0("X", 1:4), paste0("W", 1:4))],
                   ps0_wts=wts, predMu0=or0_info$predMu0, ps=ps0_info$ps)
  return(list(data=var_data, theta0ests=est, ps0_form=ps0_info$ps0_form, or0_form=or0_info$or0_form))
}

merge_sand_var_info = function(svid, svic, testD, thetaDests, D.vals=D.vals) {
  var_data = merge(svic$data, svid$data, by="ID", all=T)
  return(list(data=var_data, bw=svid$bw, theta0est=svic$theta0est,
              testD=testD, thetaDests=thetaDests, D.vals=D.vals,
              sigma=svid$sigma, orD_form=svid$orD_form, psD_form=svid$psD_form,
              or0_form=svic$or0_form, ps0_form=svic$ps0_form,
              muhat_mat=svid$muhat_mat, mhat_mat=svid$mhat_mat, gps_mat=svid$gps_mat))
}

join_EEs = function(EE, A, A_vals) {
  if (length(A_vals) > 1) return(EE)
  EE_df = data.frame(cbind(A=A, EE=matrix(0, nrow=length(A), ncol=if (is.null(ncol(EE))) 1 else ncol(EE))))
  EE_df[EE_df$A == A_vals,-1] = EE
  return(EE_df[,-1])
}

get_sand_vars = function(sand_var_info) {
  var_info1 = var_info2 = var_info3 = var_info4 = var_info5 = var_info6 = list(psiD=c(), thetaD=c(), theta0=c())
  
  data = sand_var_info$data
  trt_data = data[data$A == 1,]
  ctl_data = data[data$A == 0,]
  n = nrow(data)
  n1 = nrow(trt_data)
  p1 = mean(data$A)
  
  X_psD = create_X(data=trt_data, form=sand_var_info$psD_form)
  X_ps0 = create_X(data=data, form=sand_var_info$ps0_form)
  X_ps0_trt = X_ps0[data$A == 1,]
  X_ps0_ctl = X_ps0[data$A == 0,]
  X_or0 = create_X(data=ctl_data, form=sand_var_info$or0_form)
  trt_Xstar = create_X(data=trt_data, form=sand_var_info$orD_form)
  trt_XstarD = create_XstarD(data=trt_data, n1=n1, form=sand_var_info$orD_form)
  
  theta00_info = EE_theta00_info(data=ctl_data, n=n, p1=p1, X_psD=X_psD, X_ps0=X_ps0_ctl, X_or0=X_or0, Xstar=trt_Xstar)
  theta01_info = EE_theta01_info(data=trt_data, n=n, p1=p1, X_psD=X_psD, X_ps0=X_ps0_trt, X_or0=X_or0, Xstar=trt_Xstar)
  alpha1_info = EE_alpha1_info(data=trt_data, X_psD=X_psD, X_ps0=X_ps0_trt, X_or0=X_or0, Xstar=trt_Xstar)
  sigma_info = EE_sigma_info(data=trt_data, sigma=sand_var_info$sigma, X_psD=X_psD, X_ps0=X_ps0_trt, X_or0=X_or0, Xstar=trt_Xstar)
  lambda1_info = EE_lambda1_info(data=trt_data, X_psD=X_psD, X_ps0=X_ps0_trt, X_or0=X_or0, Xstar=trt_Xstar)
  alpha0_info = EE_alpha0_info(data=data, X_psD=X_psD, X_ps0=X_ps0, X_or0=X_or0, Xstar=trt_Xstar)
  lambda0_info = EE_lambda0_info(data=ctl_data, X_psD=X_psD, X_ps0=X_ps0_ctl, X_or0=X_or0, Xstar=trt_Xstar)
  for (i in 1:length(sand_var_info$testD)) {
    delta = sand_var_info$testD[i]
    thetaDest = sand_var_info$thetaDests[i]
    kenD_info = calc_ken_thetaD_eq_info(full_data=data, dose_val=delta, D.vals=sand_var_info$D.vals, 
                                        bw=sand_var_info$bw, muhat_mat=sand_var_info$muhat_mat, 
                                        mhat_mat=sand_var_info$mhat_mat, gps_mat=sand_var_info$gps_mat, thetaDest=thetaDest)
    thetaD1_info = EE_thetaD1_info(data=trt_data, delta=delta, thetaD=thetaDest, bw=sand_var_info$bw, 
                                   n1=n1, sigma=sand_var_info$sigma, X_psD=X_psD, X_ps0=X_ps0, X_or0=X_or0, 
                                   Xstar=trt_Xstar, XstarD=trt_XstarD)
    thetaD2_info = EE_thetaD2_info(data=trt_data, delta=delta, bw=sand_var_info$bw, n1=n1, sigma=sand_var_info$sigma,  
                                   X_psD=X_psD, X_ps0=X_ps0, X_or0=X_or0, Xstar=trt_Xstar, XstarD=trt_XstarD)
    thetaD3_info = thetaD2_info
    thetaD3_info$EE$beta1 = kenD_info$EE$beta1
    thetaD3_info$EE$beta2 = kenD_info$EE$beta2
    thetaD3_info$dEEs$beta1$dthetaD2 = list(beta1 = kenD_info$dEEs$beta1[1], beta2 = kenD_info$dEEs$beta1[2])
    thetaD3_info$dEEs$beta2$dthetaD2 = list(beta1 = kenD_info$dEEs$beta2[1], beta2 = kenD_info$dEEs$beta2[2])
    dose_var_info1 = calc_sand_var(EEs=rbind(join_EEs(thetaD1_info$EE, data$A, 1), 
                                             join_EEs(theta00_info$EE, data$A, 0),
                                             join_EEs(theta01_info$EE, data$A, 1),
                                             t(join_EEs(alpha1_info$EE, data$A, 1)),
                                             join_EEs(sigma_info$EE, data$A, 1), 
                                             t(join_EEs(lambda1_info$EE, data$A, 1)), 
                                             t(join_EEs(alpha0_info$EE, data$A, c(0,1))), 
                                             t(join_EEs(lambda0_info$EE, data$A, 0))),
                                   dEEs=rbind(set_dEEs(thetaD1_info$dEEs, v2=F), set_dEEs(theta00_info$dEEs, v2=F), 
                                              set_dEEs(theta01_info$dEEs, v2=F), set_dEEs(alpha1_info$dEEs, v2=F), 
                                              set_dEEs(sigma_info$dEEs, v2=F), set_dEEs(lambda1_info$dEEs, v2=F), 
                                              set_dEEs(alpha0_info$dEEs, v2=F), set_dEEs(lambda0_info$dEEs, v2=F)))
    dose_var_info2 = calc_sand_var(EEs=rbind(join_EEs(thetaD2_info$EE$beta1, data$A, 1),
                                             join_EEs(theta00_info$EE, data$A, 0), 
                                             join_EEs(theta01_info$EE, data$A, 1), 
                                             join_EEs(thetaD2_info$EE$beta2, data$A, 1),
                                             t(join_EEs(alpha1_info$EE, data$A, 1)), 
                                             join_EEs(sigma_info$EE, data$A, 1), 
                                             t(join_EEs(lambda1_info$EE, data$A, 1)), 
                                             t(join_EEs(alpha0_info$EE, data$A, c(0,1))), 
                                             t(join_EEs(lambda0_info$EE, data$A, 0))),
                                   dEEs=rbind(set_dEEs(thetaD2_info$dEEs$beta1, v2=T), set_dEEs(theta00_info$dEEs, v2=T), 
                                              set_dEEs(theta01_info$dEEs, v2=T), set_dEEs(thetaD2_info$dEEs$beta2, v2=T),
                                              set_dEEs(alpha1_info$dEEs, v2=T), set_dEEs(sigma_info$dEEs, v2=T), 
                                              set_dEEs(lambda1_info$dEEs, v2=T), set_dEEs(alpha0_info$dEEs, v2=T), 
                                              set_dEEs(lambda0_info$dEEs, v2=T)))
    dose_var_info3 = calc_sand_var(EEs=rbind(join_EEs(thetaD3_info$EE$beta1, data$A, 1),
                                             join_EEs(theta00_info$EE, data$A, 0), 
                                             join_EEs(theta01_info$EE, data$A, 1), 
                                             join_EEs(thetaD3_info$EE$beta2, data$A, 1),
                                             t(join_EEs(alpha1_info$EE, data$A, 1)), 
                                             join_EEs(sigma_info$EE, data$A, 1), 
                                             t(join_EEs(lambda1_info$EE, data$A, 1)), 
                                             t(join_EEs(alpha0_info$EE, data$A, c(0,1))), 
                                             t(join_EEs(lambda0_info$EE, data$A, 0))),
                                   dEEs=rbind(set_dEEs(thetaD3_info$dEEs$beta1, v2=T), set_dEEs(theta00_info$dEEs, v2=T), 
                                              set_dEEs(theta01_info$dEEs, v2=T), set_dEEs(thetaD3_info$dEEs$beta2, v2=T),
                                              set_dEEs(alpha1_info$dEEs, v2=T), set_dEEs(sigma_info$dEEs, v2=T), 
                                              set_dEEs(lambda1_info$dEEs, v2=T), set_dEEs(alpha0_info$dEEs, v2=T), 
                                              set_dEEs(lambda0_info$dEEs, v2=T)))
    dose_var_info4 = calc_sand_var(EEs=rbind(join_EEs(thetaD1_info$EE, data$A, 1), 
                                             join_EEs(theta00_info$EE, data$A, 0),
                                             join_EEs(theta01_info$EE, data$A, 1)),
                                   dEEs=rbind(set_dEEs2(thetaD1_info$dEEs, v2=F), set_dEEs2(theta00_info$dEEs, v2=F), 
                                              set_dEEs2(theta01_info$dEEs, v2=F)))
    dose_var_info5 = calc_sand_var(EEs=rbind(join_EEs(thetaD2_info$EE$beta1, data$A, 1),
                                             join_EEs(theta00_info$EE, data$A, 0), 
                                             join_EEs(theta01_info$EE, data$A, 1), 
                                             join_EEs(thetaD2_info$EE$beta2, data$A, 1)),
                                   dEEs=rbind(set_dEEs2(thetaD2_info$dEEs$beta1, v2=T), set_dEEs2(theta00_info$dEEs, v2=T), 
                                              set_dEEs2(theta01_info$dEEs, v2=T), set_dEEs2(thetaD2_info$dEEs$beta2, v2=T)))
    dose_var_info6 = calc_sand_var(EEs=rbind(join_EEs(thetaD3_info$EE$beta1, data$A, 1),
                                             join_EEs(theta00_info$EE, data$A, 0), 
                                             join_EEs(theta01_info$EE, data$A, 1), 
                                             join_EEs(thetaD3_info$EE$beta2, data$A, 1)),
                                   dEEs=rbind(set_dEEs2(thetaD3_info$dEEs$beta1, v2=T), set_dEEs2(theta00_info$dEEs, v2=T), 
                                              set_dEEs2(theta01_info$dEEs, v2=T), set_dEEs2(thetaD3_info$dEEs$beta2, v2=T)))
    for (lbl in names(var_info1)) var_info1[[lbl]] = c(var_info1[[lbl]], dose_var_info1[[lbl]])
    for (lbl in names(var_info2)) var_info2[[lbl]] = c(var_info2[[lbl]], dose_var_info2[[lbl]])
    for (lbl in names(var_info3)) var_info3[[lbl]] = c(var_info3[[lbl]], dose_var_info3[[lbl]])
    for (lbl in names(var_info4)) var_info4[[lbl]] = c(var_info4[[lbl]], dose_var_info4[[lbl]])
    for (lbl in names(var_info5)) var_info5[[lbl]] = c(var_info5[[lbl]], dose_var_info5[[lbl]])
    for (lbl in names(var_info6)) var_info6[[lbl]] = c(var_info6[[lbl]], dose_var_info6[[lbl]])
  }
  for (lbl in names(var_info1)) var_info1[[lbl]] = t(var_info1[[lbl]])
  for (lbl in names(var_info2)) var_info2[[lbl]] = t(var_info2[[lbl]])
  for (lbl in names(var_info3)) var_info3[[lbl]] = t(var_info3[[lbl]])
  for (lbl in names(var_info4)) var_info4[[lbl]] = t(var_info4[[lbl]])
  for (lbl in names(var_info5)) var_info5[[lbl]] = t(var_info5[[lbl]])
  for (lbl in names(var_info6)) var_info6[[lbl]] = t(var_info6[[lbl]])
  return(list(v1=var_info1, v2=var_info2, v3=var_info3, v4=var_info4, v5=var_info5, v6=var_info6))
}

calc_ken_thetaD_eq_info = function(full_data, dose_val, D.vals, bw, muhat_mat, mhat_mat, gps_mat, thetaDest) {
  data = full_data[full_data$A == 1,]
  d.std = (data$D - dose_val)/bw
  kern.std = dnorm(d.std)/bw
  beta = coef(lm(data$pseudo_out ~ d.std, weights=kern.std))
  kern.mat = matrix(rep(dnorm((D.vals-dose_val)/bw)/bw, nrow(data)), byrow=T, nrow=nrow(data))
  g2 = matrix(rep((D.vals-dose_val)/bw, nrow(data)), byrow=T, nrow=nrow(data))
  intfn1.mat = kern.mat * (muhat_mat - mhat_mat) * gps_mat
  intfn2.mat = g2 * kern.mat * (muhat_mat - mhat_mat) * gps_mat
  # int1 = apply(matrix(rep((D.vals[-1] - D.vals[-length(D.vals)])/2, nrow(data)), 
  #                     byrow=T, nrow=nrow(data))*intfn1.mat[,-1]+intfn1.mat[,-length(D.vals)], 1, sum)
  # int2 = apply(matrix(rep((D.vals[-1]-D.vals[-length(D.vals)])/2, nrow(data)),
  #                     byrow=T, nrow=nrow(data))*intfn2.mat[,-1]+intfn2.mat[,-length(D.vals)], 1, sum)
  int1 = apply(matrix(rep((D.vals[-1] - D.vals[-length(D.vals)]), nrow(data)), 
                      byrow=T, nrow=nrow(data))*intfn1.mat[,-1], 1, sum)
  int2 = apply(matrix(rep((D.vals[-1]-D.vals[-length(D.vals)]), nrow(data)),
                      byrow=T, nrow=nrow(data))*intfn2.mat[,-1], 1, sum)
  EE1 = kern.std * (data$pseudo_out - beta[1] - beta[2]*d.std) + int1
  EE2 = d.std * kern.std * (data$pseudo_out - beta[1] - beta[2]*d.std) + int2
  return(list(EE=list(beta1 = EE1, beta2 = EE2), dEEs=list(beta1=c(-sum(kern.std), -sum(kern.std * d.std)), 
                                                           beta2=c(-sum(kern.std * d.std), -sum(kern.std * (d.std^2))))))
}

calc_sand_var = function(EEs, dEEs) {
  mat1 = EEs %*% t(EEs)
  varmat = solve(dEEs) %*% mat1 %*% solve(t(dEEs))
  return(list(psiD=varmat[1,1]+varmat[2,2]+varmat[3,3]+2*varmat[2,3]-2*varmat[1,2]-2*varmat[1,3],
              thetaD=varmat[1,1], theta0=(varmat[2,2]+varmat[3,3]+2*(varmat[2,3]))))
}

get_lis = function(D, delta, bw) {
  ki = dnorm((D-delta)/bw)/bw
  sn1 = sum(ki*(D-delta))
  sn2 = sum(ki*(D-delta)^2)
  bi = ki*(sn2 - (D - delta)*sn1)
  li = bi/sum(bi)
  return(li)
}

get_dalpha1 = function(pref, deltaY, predMuD, n1, sigma, gps, X, fij) {
  V = pref*(deltaY-predMuD)/(n1*sigma^2*gps)
  W = matrix(nrow=nrow(X), ncol=ncol(X))
  for (i in 1:nrow(X)) {
    Ximinusj = t(t(X)-X[i,])
    W[i,] = t(fij[i,]) %*% Ximinusj
  }
  dalpha1 = t(V) %*% W
  return(dalpha1)
}

get_dsigma = function(pref, deltaY, predMuD, n1, sigma, gps, D, predD, fij) {
  V  = pref*(deltaY-predMuD)/(n1*sigma^3*gps)
  Dmat = matrix(rep(D, n1), nrow=n1, ncol=n1, byrow=F)
  predDmat = matrix(rep(predD, n1), nrow=n1, ncol=n1, byrow=T)
  Dminus = (Dmat - predDmat)^2 - (Dmat - predDmat)^2
  W = colSums(t(fij) * Dminus)
  dsigma = t(V) %*% W
  return(dsigma)
}

get_dlambda1 = function(pref, psD_wts, Xstar, n1, XstarD) {
  dlambda1 = -t(pref) %*% (psD_wts * Xstar) - 1/n1*t(colSums(pref*XstarD))
  return(dlambda1)
}

set_dEEs = function(dlist, v2=F) {
  dEEs = cbind(if (v2) dlist$dthetaD2$beta1 else dlist$dthetaD1, dlist$dtheta00, dlist$dtheta01, if (v2) dlist$dthetaD2$beta2, 
               dlist$dalpha1, dlist$dsigma, dlist$dlambda1, dlist$dalpha0, dlist$dlambda0)
  return(dEEs)
}

set_dEEs2 = function(dlist, v2=F) {
  dEEs = cbind(if (v2) dlist$dthetaD2$beta1 else dlist$dthetaD1, dlist$dtheta00, dlist$dtheta01, if (v2) dlist$dthetaD2$beta2)
  return(dEEs)
}

create_X = function(data, form) {
  X = model.matrix(form, data=data)
  return(X)
}

create_XstarD = function(data, n1, form) {
  data2 = cbind(data[rep(1:n1, n1),!(names(data) %in% c("D", "D2", "D3"))], D=rep(data$D, each=n1), 
                D2=rep(data$D^2, each=n1), D3=rep(data$D^3, each=n1))
  XstarD = model.matrix(form, data=data2)
  return(XstarD)
}

create_fij = function(D, predD, n1, sigma) {
  return(matrix(dnorm(x=rep(D, each=n1), mean=rep(D, n1), sd=sigma), nrow=n1, ncol=n1, byrow=T))
}

EE_thetaD1_info = function(data, delta, thetaD, bw, n1, sigma, X_psD, X_ps0, X_or0, Xstar, XstarD) {
  lis = get_lis(D=data$D, delta=delta, bw=bw)
  EE = thetaD/n1 - lis*data$pseudo_out
  pref = -1*lis
  fij = create_fij(data$D, data$predD, n1, sigma)
  dEEs = list(dthetaD1=1, dthetaD2=NA, dtheta00=0, dtheta01=0, 
              dalpha1=get_dalpha1(pref=pref, deltaY=data$deltaY, predMuD=data$predMuD, 
                                  n1=n1, sigma=sigma, gps=data$gps, X=X_psD, fij=fij), 
              dsigma=get_dsigma(pref=pref, deltaY=data$deltaY, predMuD=data$predMuD, n1=n1, 
                                sigma=sigma, gps=data$gps, D=data$D, predD=data$predD, fij=fij),
              dlambda1=get_dlambda1(pref=pref, psD_wts=data$psD_wts, Xstar=Xstar, n1=n1, XstarD=XstarD),
              dalpha0=matrix(0, nrow=1, ncol=ncol(X_ps0)), dlambda0=matrix(0, nrow=1, ncol=ncol(X_or0)))
  return(list(EE = EE, dEEs = dEEs))
}

get_khas = function(D.std, bw) {
  return(dnorm(D.std)/bw)
}

get_betas_kern = function(pseudo_out, D.std, khas) {
  return(coef(lm(pseudo_out ~ D.std, weights=khas)))
}

EE_thetaD2_info = function(data, delta, bw, n1, sigma, X_psD, X_ps0, X_or0, Xstar, XstarD) {
  D.std = (data$D-delta)/bw
  khas = get_khas(D.std=D.std, bw=bw)
  betas = get_betas_kern(pseudo_out=data$pseudo_out, D.std=D.std, khas=khas)
  
  EE1 = khas*(data$pseudo_out - betas[1] - D.std*betas[2])
  EE2 = EE1*D.std
  
  pref1 = khas
  fij = create_fij(data$D, data$predD, n1, sigma)
  dEEs1 = list(dthetaD1=NA, dthetaD2=list(beta1=sum(-1*khas), beta2=sum(-1*khas*D.std)), 
               dtheta00=0, dtheta01=0, 
               dalpha1=get_dalpha1(pref=pref1, deltaY=data$deltaY, predMuD=data$predMuD, 
                                   n1=n1, sigma=sigma, gps=data$gps, X=X_psD, fij=fij), 
               dsigma=get_dsigma(pref=pref1, deltaY=data$deltaY, predMuD=data$predMuD, n1=n1, 
                                 sigma=sigma, gps=data$gps, D=data$D, predD=data$predD, fij=fij),
               dlambda1=get_dlambda1(pref=pref1, psD_wts=data$psD_wts, Xstar=Xstar, n1=n1, XstarD=XstarD),
               dalpha0=matrix(0, nrow=1, ncol=ncol(X_ps0)), dlambda0=matrix(0, nrow=1, ncol=ncol(X_or0)))
  
  pref2 = khas * D.std
  dEEs2 = list(dthetaD1=NA, dthetaD2=list(beta1=sum(-1*khas*D.std), beta2=sum(-1*khas*D.std^2)), 
               dtheta00=0, dtheta01=0, 
               dalpha1=get_dalpha1(pref=pref2, deltaY=data$deltaY, predMuD=data$predMuD, 
                                   n1=n1, sigma=sigma, gps=data$gps, X=X_psD, fij=fij), 
               dsigma=get_dsigma(pref=pref2, deltaY=data$deltaY, predMuD=data$predMuD, n1=n1, 
                                 sigma=sigma, gps=data$gps, D=data$D, predD=data$predD, fij=fij),
               dlambda1=get_dlambda1(pref=pref2, psD_wts=data$psD_wts, Xstar=Xstar, n1=n1, XstarD=XstarD),
               dalpha0=matrix(0, nrow=1, ncol=ncol(X_ps0)), dlambda0=matrix(0, nrow=1, ncol=ncol(X_or0)))
  return(list(EE = list(beta1=EE1, beta2=EE2), dEEs = list(beta1=dEEs1, beta2=dEEs2)))
}

EE_theta00_info = function(data, n, p1, X_psD, X_ps0, X_or0, Xstar) {
  theta00i = data$ps0_wt/(n*p1)*(data$deltaY-data$predMu0)
  EE = mean(theta00i) - theta00i
  dEEs = list(dthetaD1=0, dthetaD2=list(beta1=0, beta2=0), dtheta00=1, dtheta01=0, 
              dalpha1=matrix(0, nrow=1, ncol=ncol(X_psD)), dsigma=0, dlambda1=matrix(0, nrow=1, ncol=ncol(Xstar)),
              dalpha0=t(colSums(-1*data$ps0_wts*(data$deltaY-data$predMu0)/(n*p1)*X_ps0)), 
              dlambda0=t(colSums(data$ps0_wts/(n*p1)*X_or0)))
  return(list(EE = EE, dEEs = dEEs))
  
}

EE_theta01_info = function(data, n, p1, X_psD, X_ps0, X_or0, Xstar) {
  theta01i = data$predMu0/(n*p1)
  EE = mean(theta01i) - theta01i
  # EE = mean(data$predMu0)/p1 - data$predMu0/p1
  dEEs = list(dthetaD1=0, dthetaD2=list(beta1=0, beta2=0), dtheta00=0, dtheta01=1, 
              dalpha1=matrix(0, nrow=1, ncol=ncol(X_psD)), dsigma=0, dlambda1=matrix(0, nrow=1, ncol=ncol(Xstar)),
              dalpha0=matrix(0, nrow=1, ncol=ncol(X_ps0)), dlambda0=t(colSums(-1*X_or0/(n*p1))))
  return(list(EE = EE, dEEs = dEEs))
}

EE_alpha1_info = function(data, X_psD, X_ps0, X_or0, Xstar) {
  EE = X_psD*(data$D - data$predD)
  nr = ncol(X_psD)
  dEEs = list(dthetaD1=rep(0, nr), dthetaD2=list(beta1=rep(0, nr), beta2=rep(0, nr)), dtheta00=rep(0, nr), dtheta01=rep(0, nr),
              dalpha1=-1*t(X_psD)%*%X_psD, dsigma=rep(0, nr),  dlambda1=matrix(0, nrow=nr, ncol=ncol(Xstar)),
              dalpha0=matrix(0, nrow=nr, ncol=ncol(X_ps0)), dlambda0=matrix(0, nrow=nr, ncol=ncol(X_or0)))
  return(list(EE = EE, dEEs = dEEs))
}

EE_sigma_info = function(data, sigma, X_psD, X_ps0, X_or0, Xstar) {
  EE = (data$D-data$predD)^2/sigma^2 - 1
  dEEs = list(dthetaD1=0, dthetaD2=list(beta1=0, beta2=0), dtheta00=0, dtheta01=0, 
              dalpha1=t(colSums(2*(data$D-data$predD)/sigma^2*X_psD)), 
              dsigma=sum(-2*(data$D-data$predD)^2/(sigma^3)), 
              dlambda1=matrix(0, nrow=1, ncol=ncol(Xstar)),
              dalpha0=matrix(0, nrow=1, ncol=ncol(X_ps0)), dlambda0=matrix(0, nrow=1, ncol=ncol(X_or0)))
  return(list(EE = EE, dEEs = dEEs))
}

EE_lambda1_info = function(data, X_psD, X_ps0, X_or0, Xstar) {
  EE = Xstar*(data$deltaY-data$predMuD)
  nr = ncol(Xstar)
  dEEs = list(dthetaD1=rep(0, nr), dthetaD2=list(beta1=rep(0, nr), beta2=rep(0, nr)), dtheta00=rep(0, nr), dtheta01=rep(0, nr), 
              dalpha1=matrix(0, nrow=nr, ncol=ncol(X_psD)), dsigma=rep(0, nr), dlambda1=-1*t(Xstar)%*%Xstar,
              dalpha0=matrix(0, nrow=nr, ncol=ncol(X_ps0)), dlambda0=matrix(0, nrow=nr, ncol=ncol(X_or0)))
  return(list(EE = EE, dEEs = dEEs))
}

EE_alpha0_info = function(data, X_psD, X_ps0, X_or0, Xstar) {
  EE = X_ps0*(data$A - data$ps)
  nr = ncol(X_ps0)
  dEEs = list(dthetaD1=rep(0, nr), dthetaD2=list(beta1=rep(0, nr), beta2=rep(0, nr)), dtheta00=rep(0, nr), dtheta01=rep(0, nr), 
              dalpha1=matrix(0, nrow=nr, ncol=ncol(X_psD)), dsigma=rep(0, nr), dlambda1=matrix(0, nrow=nr, ncol=ncol(Xstar)),
              dalpha0=t(-1*data$ps*(1-data$ps)*X_ps0) %*% X_ps0, dlambda0=matrix(0, nrow=nr, ncol=ncol(X_or0)))
  return(list(EE = EE, dEEs = dEEs))
}

EE_lambda0_info = function(data, X_psD, X_ps0, X_or0, Xstar) {
  EE = X_or0*(data$deltaY-data$predMu0)
  nr = ncol(X_or0)
  dEEs = list(dthetaD1=rep(0, nr), dthetaD2=list(beta1=rep(0, nr), beta2=rep(0, nr)), dtheta00=rep(0, nr), dtheta01=rep(0, nr), 
              dalpha1=matrix(0, nrow=nr, ncol=ncol(X_psD)), dsigma=rep(0, nr), dlambda1=matrix(0, nrow=nr, ncol=ncol(Xstar)),
              dalpha0=matrix(0, nrow=nr, ncol=ncol(X_ps0)), dlambda0=-1*t(X_or0)%*%X_or0)
  return(list(EE = EE, dEEs = dEEs))
}
