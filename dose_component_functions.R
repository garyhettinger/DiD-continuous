approx_fxn = function(x,y,newx) { predict(smooth.spline(x,y),x=newx)$y}

get_psD_info = function(data, ps_good, D.vals, dose_col="D", good_covs=paste0("X", 1:4), 
                        bad_covs = paste0("W", 1:4), is_boot=F, wt_col="boot_wts") {
  data$boot_wts = if (!is_boot) rep(1, nrow(data)) else data[[wt_col]]
  psD_form = as.formula(paste0(dose_col, " ~ ", paste0(if (ps_good) good_covs else bad_covs, collapse=" + ")))
  psD_model = lm(psD_form, data=data, weights=boot_wts)
  psD_sigma = summary(psD_model)$sigma
  piD_mean = predict(psD_model, newdata=data, type="response")
  piD_denom = dnorm(x=data[[dose_col]], mean=piD_mean, sd=psD_sigma)
  piD_num = approx_fxn(x=D.vals, y=sapply(D.vals, function(d) 
    stats::weighted.mean(x=dnorm(x=rep(d, nrow(data)), mean=piD_mean, sd=psD_sigma), w=data$boot_wts)), newx=data[[dose_col]])
  gps_mat1 = matrix(dnorm(x=rep(D.vals, each=nrow(data)), mean=rep(piD_mean, length(D.vals)), sd=psD_sigma), nrow=nrow(data), ncol=length(D.vals))
  gps_mat2 = matrix(rep(apply(gps_mat1, 2, mean), nrow(data)), byrow=T, nrow=nrow(data))
  return(list(psD_form=psD_form, sigma=psD_sigma, psD_model=psD_model, predD=piD_mean, gps=piD_denom, gps_num=piD_num,
              gps_mat=gps_mat2))
}

get_psD_wts = function(gps, gps_num, doses=NULL, normalize=T, trim=F) {
  wts = gps_num / gps
  if (trim) wts = WeightIt::trim(x=wts, at=0.95, lower=F, treat=doses)
  if (normalize) wts = wts/mean(wts)
  return(wts)
}

get_orD_info = function(data, or_good, D.vals, dose_col="D", outcome_col="deltaY", 
                        good_cov_form="X1 + X2 + X3 + X4 + D + D:X1 + D:X3 + D3", 
                        bad_cov_form="W1 + W2 + W3 + W4 + D + D:W1 + D:W3",
                        is_boot=F, wt_col="boot_wts") {
  data$boot_wts = if (!is_boot) rep(1, nrow(data)) else data[[wt_col]]
  orD_form = as.formula(paste0(outcome_col, " ~ ", if (or_good) good_cov_form else bad_cov_form))
  orD_model = glm(orD_form, data=data, family=gaussian(link="identity"), weights=boot_wts)
  muD = predict(orD_model, newdata=data, type="response")
  mD = approx_fxn(x=D.vals, y=sapply(D.vals, function(d) stats::weighted.mean(x=predict(
    orD_model, newdata=data.frame(data %>% mutate(!!dose_col := d, D3=d^3)), type="response"), 
    w=data$boot_wts)), newx=data[[dose_col]])
  test_data = data.frame(cbind( data[rep(1:nrow(data),length(D.vals)), !(names(data) %in% c("D", "D2", "D3"))],
                                D=rep(D.vals,rep(nrow(data),length(D.vals)))))
  test_data$D2 = test_data$D^2
  test_data$D3 = test_data$D^3
  muD_mat = matrix(predict(orD_model, newdata=test_data), nrow=nrow(data), ncol=length(D.vals))
  mD_mat = matrix(rep(apply(muD_mat, 2, mean), nrow(data)), byrow=T, nrow=nrow(data))
  return(list(orD_form=orD_form, orD_model=orD_model, predMuD=muD, predMD=mD, muhat_mat=muD_mat, mhat_mat=mD_mat))
}

get_pseudo_out_info = function(data, wts, predMuD, predMD, outcome_col="deltaY") {
  pseudo.dr = wts*(data[[outcome_col]] - predMuD) + predMD
  pseudo.ipw = wts*data[[outcome_col]]
  return(list(pseudo_dr=pseudo.dr, pseudo_ipw=pseudo.ipw))
}

get_thetaD_est_twfe = function(data, D.vals, or_good, good_covs=paste0("X", 1:4), 
                               bad_covs=paste0("W", 1:4), outcome_col0="Y0", outcome_col1="Y1", 
                               trt_col="A", dose_col="D", is_boot=F, wt_col="boot_wts") {
  data$boot_wts = if (!is_boot) rep(1, nrow(data)) else data[[wt_col]]
  form = as.formula(paste0("Y ~ ", paste0(c(if (or_good) good_covs else bad_covs, "Time", 
                                            trt_col, paste0(trt_col, ":Time"), paste0(trt_col, ":", dose_col, ":Time")), collapse=" + ")))
  stacked_data1 = data %>% select(-!!sym(outcome_col0)) %>% rename(Y=outcome_col1) %>% mutate(Time=1)
  stacked_data0 = data %>% select(-!!sym(outcome_col1)) %>% rename(Y=outcome_col0) %>% mutate(Time=0)
  stacked_data = rbind(stacked_data1, stacked_data0)
  stacked_data[[dose_col]] = ifelse(stacked_data[[trt_col]] == 1, stacked_data[[dose_col]], 0)
  m1 = lm(form, data=stacked_data, weights=boot_wts)
  est = m1$coefficients[[paste0("Time:", trt_col, ":", dose_col)]]*D.vals
  return(est)
}

get_thetaD_est_naive = function(data, D.vals, outcome_col="deltaY", dose_col="D", is_boot=F, wt_col="boot_wts") {
  kern_wts = if (!is_boot) rep(1, nrow(data)) else data[[wt_col]]
  return(wt_kernel_est(data=data, dose_vals=D.vals, outcome_col=outcome_col, 
                       dose_col=dose_col, kern_wts=kern_wts)$est)
}

get_thetaD_est_or = function(data, D.vals, orD_model, dose_col="D", is_boot=F, wt_col="boot_wts") {
  data$boot_wts = if (!is_boot) rep(1, nrow(data)) else data[[wt_col]]
  return(sapply(D.vals, function(d) stats::weighted.mean(x=predict(orD_model, newdata=data.frame(data %>% mutate(!!dose_col := d, D3=d^3)), 
                                                                   type="response"), w=data$boot_wts)))
}

get_thetaD_est_ipw = function(data, D.vals, pseudo, dose_col="D", is_boot=F, wt_col="boot_wts") {
  kern_wts = if (!is_boot) rep(1, nrow(data)) else data[[wt_col]]
  data$pseudo_ipw = pseudo
  return(wt_kernel_est(data=data, dose_vals=D.vals, outcome_col="pseudo_ipw", 
                       dose_col=dose_col, kern_wts=kern_wts)$est)
}

get_thetaD_est_dr_info = function(data, D.vals, pseudo, dose_col="D", is_boot=F, wt_col="boot_wts") {
  kern_wts = if (!is_boot) rep(1, nrow(data)) else data[[wt_col]]
  data$pseudo_out = pseudo
  return(wt_kernel_est(data=data, dose_vals=D.vals, outcome_col="pseudo_out", 
                       dose_col=dose_col, kern_wts=kern_wts))
}

get_thetaD_est_parametric = function(data, D.vals, pseudo, form, dose_col="D", is_boot=F, wt_col="boot_wts") {
  data$boot_wts = if (!is_boot) rep(1, nrow(data)) else data[[wt_col]]
  data$pseudo_dr = pseudo
  form = as.formula(paste0("pseudo_dr ~ ", form))
  m1 = lm(form, data=data, weights=boot_wts)
  newdata = data.frame(dose_holder=D.vals, D2=D.vals^2, D3=D.vals^3)
  names(newdata)[1] = dose_col
  est = predict(m1, newdata=newdata, type="response")
  return(est)
}

get_thetaD_ests = function(data, D.vals, or_good, ps_good, is_boot=F, wt_col="boot_wts", wt_norm=T) {
  psD_info = get_psD_info(data=data, ps_good=ps_good, D.vals=D.vals, is_boot=is_boot, wt_col=wt_col)
  orD_info = get_orD_info(data=data, or_good=or_good, D.vals=D.vals, is_boot=is_boot, wt_col=wt_col)
  psD_wts = get_psD_wts(gps=psD_info$gps, gps_num=psD_info$gps_num, normalize=wt_norm)
  pseudo_info = get_pseudo_out_info(data=data, wts=psD_wts, predMuD=orD_info$predMuD, predMD=orD_info$predMD)
  
  twfe_est = get_thetaD_est_twfe(data=data, D.vals=D.vals, or_good=or_good, is_boot=is_boot, wt_col=wt_col)
  naive_est = get_thetaD_est_naive(data=data, D.vals=D.vals, is_boot=is_boot, wt_col=wt_col)
  or_est = get_thetaD_est_or(data=data, D.vals=D.vals, orD_model=orD_info$orD_model, is_boot=is_boot, wt_col=wt_col)
  ipw_est = get_thetaD_est_ipw(data=data, D.vals=D.vals, pseudo=pseudo_info$pseudo_ipw, is_boot=is_boot, wt_col=wt_col)
  dr_est_info = get_thetaD_est_dr_info(data=data, D.vals=D.vals, pseudo=pseudo_info$pseudo_dr, is_boot=is_boot, wt_col=wt_col)
  dr_est = dr_est_info$est
  dr_paramc_est = get_thetaD_est_parametric(data=data, D.vals=D.vals, pseudo=pseudo_info$pseudo_dr, 
                                            form="D + D3", is_boot=is_boot, wt_col=wt_col)
  dr_parami_est = get_thetaD_est_parametric(data=data, D.vals=D.vals, pseudo=pseudo_info$pseudo_dr, 
                                            form="D + D2", is_boot=is_boot, wt_col=wt_col)
  ests = list(twfe=twfe_est, naive=naive_est, or=or_est, ipw=ipw_est, dr=dr_est, dr_paramc=dr_paramc_est, dr_parami=dr_parami_est)
  sand_var_info = get_sand_var_info_dose(data=data, psD_info=psD_info, orD_info=orD_info, wts=psD_wts,
                                         pseudo_info=pseudo_info, bw=dr_est_info$bw,
                                         muhat_mat=orD_info$muhat_mat, mhat_mat=orD_info$mhat_mat, gps_mat=psD_info$gps_mat)
  return(list(ests=ests, sand_var_info=sand_var_info))
}