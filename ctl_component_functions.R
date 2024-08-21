get_ps0_info = function(data, ps_good, trt_col="A", good_covs=paste0("X", 1:4), bad_covs=paste0("W", 1:4),
                        is_boot=F, wt_col="boot_wts") {
  data$boot_wts = if (!is_boot) rep(1, nrow(data)) else data[[wt_col]]
  ps0_form = as.formula(paste0(trt_col, " ~ ", paste0(if (ps_good) good_covs else bad_covs, collapse=" + ")))
  ps0_model = glm(ps0_form, data=data, family=binomial(link="logit"), weights=boot_wts)
  piA_mean = predict(ps0_model, newdata=data, type="response")
  return(list(ps0_form=ps0_form, ps0_model=ps0_model, ps=piA_mean))
}

get_ps0_wts = function(data, ps, trt_col="A", normalize=T, trim=F) {
  wts = ifelse(data[[trt_col]] == 1, 1, ps/(1-ps))
  if (trim) wts = WeightIt::trim(x=wts, at=0.95, lower=T, treat=data[[trt_col]])
  if (normalize) wts[data[[trt_col]]==0] = wts[data[[trt_col]]==0]/mean(wts[data[[trt_col]]==0])
  return(wts)
}

get_or0_info = function(data, or_good, outcome_col="deltaY", trt_col="A",
                        good_covs=paste0("X", 1:4), bad_covs=paste0("W", 1:4),
                        is_boot=F, wt_col="boot_wts") {
  data$boot_wts = if (!is_boot) rep(1, nrow(data)) else data[[wt_col]]
  or0_form = as.formula(paste0(outcome_col, " ~ ", paste0(if (or_good) good_covs else bad_covs, collapse=" + ")))
  or0_model = glm(or0_form, data=data[data[[trt_col]]==0,], family=gaussian(link="identity"), weights=boot_wts)
  mu0 = predict(or0_model, newdata=data, type="response")
  return(list(or0_form=or0_form, or0_model=or0_model, predMu0=mu0))
}

get_theta0_est_twfe = function(data, or_good, good_covs=paste0("X", 1:4), 
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
  est = -1*m1$coefficients[[paste0("Time:", trt_col)]]
  return(est)
}

get_theta0_est_naive = function(data, outcome_col="deltaY", trt_col="A", is_boot=F, wt_col="boot_wts") {
  data$boot_wts = if (!is_boot) rep(1, nrow(data)) else data[[wt_col]]
  ctl_data = data[data[[trt_col]]==0,]
  return(stats::weighted.mean(x=ctl_data[[outcome_col]], w=ctl_data$boot_wts))
}

get_theta0_est_or = function(data, predMu0, outcome_col="deltaY", trt_col="A", is_boot=F, wt_col="boot_wts") {
  data$boot_wts = if (!is_boot) rep(1, nrow(data)) else data[[wt_col]]
  return(stats::weighted.mean(x=predMu0[data[[trt_col]]==1], w=data$boot_wts[data[[trt_col]]==1]))
}

get_theta0_est_ipw = function(data, wts, outcome_col="deltaY", trt_col="A", is_boot=F, wt_col="boot_wts") {
  data$boot_wts = if (!is_boot) rep(1, nrow(data)) else data[[wt_col]]
  return(stats::weighted.mean(x=wts[data[[trt_col]]==0]*data[[outcome_col]][data[[trt_col]]==0], w=data$boot_wts[data[[trt_col]]==0]))
}

get_theta0_est_dr = function(data, wts, predMu0, outcome_col="deltaY", trt_col="A", is_boot=F, wt_col="boot_wts") {
  data$boot_wts = if (!is_boot) rep(1, nrow(data)) else data[[wt_col]]
  return(stats::weighted.mean(x=ifelse(data[[trt_col]]==1, wts*predMu0, wts*(data[[outcome_col]] - predMu0)), w=data$boot_wts)/
           stats::weighted.mean(x=data[[trt_col]]==1, w=data$boot_wts))
}

get_theta0_ests = function(data, or_good, ps_good, is_boot=F, wt_col="boot_wts", wt_norm=T) {
  ps0_info = get_ps0_info(data=data, ps_good=ps_good, is_boot=is_boot, wt_col=wt_col)
  or0_info = get_or0_info(data=data, or_good=or_good, is_boot=is_boot, wt_col=wt_col)
  ps0_wts = get_ps0_wts(data=data, ps=ps0_info$ps, normalize=wt_norm)
  
  twfe_est = get_theta0_est_twfe(data=data, or_good=or_good, is_boot=is_boot, wt_col=wt_col) 
  naive_est = get_theta0_est_naive(data=data, is_boot=is_boot, wt_col=wt_col)
  or_est = get_theta0_est_or(data=data, predMu0=or0_info$predMu0, is_boot=is_boot, wt_col=wt_col)
  ipw_est = get_theta0_est_ipw(data=data, wts=ps0_wts, is_boot=is_boot, wt_col=wt_col)
  dr_est = get_theta0_est_dr(data=data, wts=ps0_wts, predMu0=or0_info$predMu0, is_boot=is_boot, wt_col=wt_col)
  ests = list(twfe=twfe_est, naive=naive_est, or=or_est, ipw=ipw_est, dr=dr_est,
              dr_paramc=dr_est, dr_parami=dr_est)
  
  sand_var_info = get_sand_var_info_ctl(data=data, ps0_info=ps0_info, or0_info=or0_info, 
                                        wts=ps0_wts, est=dr_est)
  
  return(list(ests=ests, sand_var_info=sand_var_info))
}