library(rlist)
# Function to get sandwich variance information for the dose-component of the estimation
get_sand_var_info_dose = function(data, psD_info, orD_info, wts, pseudo_info, bw,
                                  muhat_mat, mhat_mat, gps_mat) {
  var_data = cbind(ID=data$ID, psD_wts=wts, pseudo_out=pseudo_info$pseudo_dr, 
                   gps=psD_info$gps, predD=psD_info$predD, predMuD=orD_info$predMuD)
  return(list(data=var_data, bw=bw, muhat_mat=muhat_mat, mhat_mat=mhat_mat, gps_mat=gps_mat,
              orD_form=orD_info$orD_form, psD_form=psD_info$psD_form, sigma=psD_info$sigma))
}

# Function to get sandwich variance information for the binary control-component of the estimation
get_sand_var_info_ctl = function(data, ps0_info, or0_info, wts, est) {
  var_data = cbind(data[,c("ID", "deltaY", "A", "D", "D2", "D3", paste0("X", 1:4), paste0("W", 1:4))],
                   ps0_wts=wts, predMu0=or0_info$predMu0, ps=ps0_info$ps)
  return(list(data=var_data, theta0ests=est, ps0_form=ps0_info$ps0_form, or0_form=or0_info$or0_form))
}

# Helper function to create data structure with relevant information sandwich variance
merge_sand_var_info = function(svid, svic, testD, thetaDests, D.vals=D.vals) {
  var_data = merge(svic$data, svid$data, by="ID", all=T)
  return(list(data=var_data, bw=svid$bw, theta0est=svic$theta0est,
              testD=testD, thetaDests=thetaDests, D.vals=D.vals,
              sigma=svid$sigma, orD_form=svid$orD_form, psD_form=svid$psD_form,
              or0_form=svic$or0_form, ps0_form=svic$ps0_form,
              muhat_mat=svid$muhat_mat, mhat_mat=svid$mhat_mat, gps_mat=svid$gps_mat))
}

# Helper function to structure estimating equation information
join_EEs = function(EE, A, A_vals) {
  if (length(A_vals) > 1) return(EE)
  EE_df = data.frame(cbind(A=A, EE=matrix(0, nrow=length(A), ncol=if (is.null(ncol(EE))) 1 else ncol(EE))))
  EE_df[EE_df$A == A_vals,-1] = EE
  return(EE_df[,-1])
}

get_sand_vars = function(sand_var_info) {
  var_info = list(psiD=c(), thetaD=c(), theta0=c())
  
  data = sand_var_info$data
  trt_data = data[data$A == 1,]
  ctl_data = data[data$A == 0,]
  n = nrow(data)
  p1 = mean(data$A)
  
  theta00_info = EE_theta00_info(data=ctl_data, n=n, p1=p1)
  theta01_info = EE_theta01_info(data=trt_data, n=n, p1=p1)
  for (i in 1:length(sand_var_info$testD)) {
    delta = sand_var_info$testD[i]
    thetaDest = sand_var_info$thetaDests[i]
    thetaD_info = EE_thetaD_info(data=trt_data, dose_val=delta, D.vals=sand_var_info$D.vals, 
                                 bw=sand_var_info$bw, muhat_mat=sand_var_info$muhat_mat, 
                                 mhat_mat=sand_var_info$mhat_mat, gps_mat=sand_var_info$gps_mat)
    dose_var_info = calc_sand_var(EEs=rbind(join_EEs(thetaD_info$EE$beta1, data$A, 1),
                                             join_EEs(theta00_info$EE, data$A, 0), 
                                             join_EEs(theta01_info$EE, data$A, 1), 
                                             join_EEs(thetaD_info$EE$beta2, data$A, 1)),
                                   dEEs=rbind(set_dEEs(thetaD_info$dEEs$beta1), set_dEEs(theta00_info$dEEs), 
                                              set_dEEs(theta01_info$dEEs), set_dEEs(thetaD_info$dEEs$beta2)))
    for (lbl in names(var_info)) var_info[[lbl]] = c(var_info[[lbl]], dose_var_info[[lbl]])
  }
  for (lbl in names(var_info)) var_info[[lbl]] = t(var_info[[lbl]])
  return(list(var=var_info))
}

calc_sand_var = function(EEs, dEEs) {
  mat1 = EEs %*% t(EEs)
  varmat = solve(dEEs) %*% mat1 %*% solve(t(dEEs))
  return(list(psiD=varmat[1,1]+varmat[2,2]+varmat[3,3]+2*varmat[2,3]-2*varmat[1,2]-2*varmat[1,3],
              thetaD=varmat[1,1], theta0=(varmat[2,2]+varmat[3,3]+2*(varmat[2,3]))))
}


set_dEEs = function(dlist) {
  dEEs = cbind(dlist$dbeta1, dlist$dtheta00, dlist$dtheta01, dlist$dbeta2)
  return(dEEs)
}

create_X = function(data, form) {
  X = model.matrix(form, data=data)
  return(X)
}


# Function to calculate sandwich variance for estimating equations corresponding to dose-specific pseudo-outcomes
EE_thetaD_info = function(data, dose_val, D.vals, bw, muhat_mat, mhat_mat, gps_mat) {
  d.std = (data$D - dose_val)/bw
  kern.std = dnorm(d.std)/bw
  beta = coef(lm(data$pseudo_out ~ d.std, weights=kern.std))
  kern.mat = matrix(rep(dnorm((D.vals-dose_val)/bw)/bw, nrow(data)), byrow=T, nrow=nrow(data))
  g2 = matrix(rep((D.vals-dose_val)/bw, nrow(data)), byrow=T, nrow=nrow(data))
  intfn1.mat = kern.mat * (muhat_mat - mhat_mat) * gps_mat
  intfn2.mat = g2 * kern.mat * (muhat_mat - mhat_mat) * gps_mat
  int1 = apply(matrix(rep((D.vals[-1] - D.vals[-length(D.vals)]), nrow(data)), 
                      byrow=T, nrow=nrow(data))*intfn1.mat[,-1], 1, sum)
  int2 = apply(matrix(rep((D.vals[-1]-D.vals[-length(D.vals)]), nrow(data)),
                      byrow=T, nrow=nrow(data))*intfn2.mat[,-1], 1, sum)
  EE1 = kern.std * (data$pseudo_out - beta[1] - beta[2]*d.std) + int1
  EE2 = d.std * kern.std * (data$pseudo_out - beta[1] - beta[2]*d.std) + int2
  EE=list(beta1 = EE1, beta2 = EE2)
  dEEs=list(beta1=list(dbeta1=-sum(kern.std), dbeta2=-sum(kern.std * d.std), dtheta00=0, dtheta01=0),
            beta2=list(dbeta1=-sum(kern.std * d.std), dbeta2=-sum(kern.std * (d.std^2)), dtheta00=0, dtheta01=0))
  return(list(EE = EE, dEEs = dEEs))
}

EE_theta00_info = function(data, n, p1) {
  theta00i = data$ps0_wt/(n*p1)*(data$deltaY-data$predMu0)
  EE = mean(theta00i) - theta00i
  dEEs = list(dbeta1=0, dbeta2=0, dtheta00=1, dtheta01=0)
  return(list(EE = EE, dEEs = dEEs))
  
}

EE_theta01_info = function(data, n, p1) {
  theta01i = data$predMu0/(n*p1)
  EE = mean(theta01i) - theta01i
  dEEs = list(dbeta1=0, dbeta2=0, dtheta00=0, dtheta01=1)
  return(list(EE = EE, dEEs = dEEs))
}
