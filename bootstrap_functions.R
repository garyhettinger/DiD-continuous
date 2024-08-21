get_boot_wts = function(seed, data) {
  #-----------------------------------------------------------------------------
  set.seed(seed)
  wts <- rexp(n=nrow(data)) # MCMCpack::rdirichlet(n=1, alpha=rep(1,nrow(data)))
  # normalize per treatment group (e.g., stratified bootstrap)
  wts[data$A == 1] <- wts[data$A == 1] * sum(data$A == 1) / sum(wts[data$A == 1])
  wts[data$A == 0] <- wts[data$A == 0] * sum(data$A == 0) / sum(wts[data$A == 0])
  #-----------------------------------------------------------------------------
  return(c(wts))
}

collect_sim_boots = function(data, nboots, D.vals, testD, or_good_dose, ps_good_dose, 
                             or_good_ctl, ps_good_ctl, wt_normD=T, wt_norm0=T) {
  boot_res = NULL
  for (i in 1:nboots) {
    data$boot_wts = get_boot_wts(i, data)
    est_info = get_est_info(data=data, D.vals=D.vals, testD=testD, trt_col="A", 
                            or_good_dose=or_good_dose, ps_good_dose=ps_good_dose, 
                            or_good_ctl=or_good_ctl, ps_good_ctl=ps_good_ctl, is_boot=T, 
                            wt_normD=wt_normD, wt_norm0=wt_norm0)
    est_dfs = list(psiD=data.frame(est_info$psiD), thetaD=data.frame(est_info$thetaD), theta0=data.frame(est_info$theta0))
    if (is.null(boot_res)) boot_res = est_dfs
    else for (est in names(est_dfs)) boot_res[[est]] = rbind(boot_res[[est]], est_dfs[[est]])
  }
  agg_res = lapply(boot_res, function(x) data.frame(cbind(t(apply(x, 2, function(y) mean(y, na.rm=T))), 
                                                          t(apply(x, 2, function(y) sd(y, na.rm=T))),
                                                          t(apply(x, 2, function(y) quantile(y, probs=0.025, na.rm=T))),
                                                          t(apply(x, 2, function(y) quantile(y, probs=0.975, na.rm=T))))))
  for (est in names(boot_res)) names(agg_res[[est]]) = c(paste0(names(boot_res[[est]]), ".mean"), paste0(names(boot_res[[est]]), ".sd"), 
                                                         paste0(names(boot_res[[est]]), ".lower"), paste0(names(boot_res[[est]]), ".upper"))
  return(agg_res)
}