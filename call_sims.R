library(dplyr)
source("sim_generation.R")
source("kernel_functions.R")
source("dose_component_functions.R")
source("ctl_component_functions.R")
source("bootstrap_functions.R")
source("variance_functions.R")

flatten_list_res = function(ests_list) {
  ests_df = NULL
  for (mdl in names(ests_list)) {
    mdl_ests = data.frame(t(ests_list[[mdl]]))
    names(mdl_ests) = paste0(mdl, 1:ncol(mdl_ests))
    ests_df = if (is.null(ests_df)) mdl_ests else cbind(ests_df, mdl_ests)
  }
  return(ests_df)
}

get_est_info = function(data, D.vals, testD, trt_col, or_good_dose, ps_good_dose, or_good_ctl, ps_good_ctl, 
                        is_boot=F, wt_col="boot_wts", wt_normD=T, wt_norm0=T) {
  est_dose_info = get_thetaD_ests(data=data[data[[trt_col]]==1,], D.vals=D.vals, 
                                  or_good=or_good_dose, ps_good=ps_good_dose, 
                                  is_boot=is_boot, wt_col=wt_col, wt_norm=wt_normD)
  est_ctl_info = get_theta0_ests(data=data, or_good=or_good_ctl, ps_good=ps_good_ctl,
                                 is_boot=is_boot, wt_col=wt_col, wt_norm=wt_norm0)
  thetaD = lapply(est_dose_info$ests, function(est) approx(x=D.vals, y=est, xout=testD)$y)
  theta0 = est_ctl_info$ests
  psiD = list()
  sand_var_info = merge_sand_var_info(svid=est_dose_info$sand_var_info, svic=est_ctl_info$sand_var_info, 
                                      testD=testD, thetaDests=thetaD$dr, D.vals=D.vals)
  for (model in names(thetaD)) psiD[[model]] = thetaD[[model]] - theta0[[model]]
  return(list(psiD=flatten_list_res(psiD), thetaD=flatten_list_res(thetaD), theta0=flatten_list_res(theta0), 
              sand_var_info=sand_var_info))
}

add_info = function(seed, res, info, ogd, pgd, ogc, pgc) {
  dfs = list(psiD=cbind(Seed=seed, OGD=ogd, PGD=pgd, OGC=ogc, PGC=pgc, data.frame(info$psiD)), 
             thetaD=cbind(Seed=seed, OGD=ogd, PGD=pgd, OGC=ogc, PGC=pgc, data.frame(info$thetaD)), 
             theta0=cbind(Seed=seed, OGD=ogd, PGD=pgd, OGC=ogc, PGC=pgc, data.frame(info$theta0)))
  if (is.null(res)) res = dfs
  else for (lbl in names(dfs)) res[[lbl]] = rbind(res[[lbl]], dfs[[lbl]])
  return(res)
}

collect_point_sims = function(start_sim, end_sim, n, testD, wt_normD=T, wt_norm0=T) {
  res = list()
  for (seed in start_sim:end_sim) {
    if (seed %% 1 == 0) print(seed)
    sim_info = gen_sim(seed=seed, n=n)
    for (ogd in c(T,F)) {
      for (pgd in c(T,F)) {
        for (ogc in c(T,F)) {
          for (pgc in c(T,F)) {
            est_info = get_est_info(data=sim_info$data, D.vals=sim_info$D.vals, testD=testD, trt_col="A", 
                                    or_good_dose=ogd, ps_good_dose=pgd, 
                                    or_good_ctl=ogc, ps_good_ctl=pgc, wt_normD=wt_normD, wt_norm0=wt_norm0)
            res = add_info(seed=seed, res=res, info=est_info, ogd=ogd, pgd=pgd, ogc=ogc, pgc=pgc)
          }
        }
      }
    }
  }
  return(res)
}

run_ci_sims = function(start_sim, end_sim, n, testD, nboots, get_vars=T, wt_normD=T, wt_norm0=T) {
  point_res = sand_var_res1 = sand_var_res2 = sand_var_res3 = sand_var_res4 = sand_var_res5 = sand_var_res6 = boot_res = NULL
  for (seed in start_sim:end_sim) {
    if (seed %% 1 == 0) print(seed)
    sim_info = gen_sim(seed=seed, n=n)
    for (og in c(T,F)) {
      for (pg in c(T,F)) {
        print("new run")
        est_info = get_est_info(data=sim_info$data, D.vals=sim_info$D.vals, testD=testD, trt_col="A", 
                                or_good_dose=og, ps_good_dose=pg, 
                                or_good_ctl=og, ps_good_ctl=pg, wt_normD=wt_normD, wt_norm0=wt_norm0)
        point_res = add_info(seed=seed, res=point_res, info=est_info, ogd=og, pgd=pg, ogc=og, pgc=pg)
        
        if (get_vars) {
          sand_var_res_info = get_sand_vars(est_info$sand_var_info)
          sand_var_res1 = add_info(seed=seed, res=sand_var_res1, info=sand_var_res_info$v1, ogd=og, pgd=pg, ogc=og, pgc=pg)
          sand_var_res2 = add_info(seed=seed, res=sand_var_res2, info=sand_var_res_info$v2, ogd=og, pgd=pg, ogc=og, pgc=pg)
          sand_var_res3 = add_info(seed=seed, res=sand_var_res3, info=sand_var_res_info$v3, ogd=og, pgd=pg, ogc=og, pgc=pg)
          sand_var_res4 = add_info(seed=seed, res=sand_var_res4, info=sand_var_res_info$v4, ogd=og, pgd=pg, ogc=og, pgc=pg)
          sand_var_res5 = add_info(seed=seed, res=sand_var_res5, info=sand_var_res_info$v5, ogd=og, pgd=pg, ogc=og, pgc=pg)
          sand_var_res6 = add_info(seed=seed, res=sand_var_res6, info=sand_var_res_info$v6, ogd=og, pgd=pg, ogc=og, pgc=pg)
        }
        
        if (nboots > 0) {
          boot_res_info = collect_sim_boots(data=sim_info$data, nboots=nboots, D.vals=sim_info$D.vals, 
                                            testD=testD, or_good_dose=og, ps_good_dose=pg, 
                                            or_good_ctl=og, ps_good_ctl=pg, wt_normD=wt_normD, wt_norm0=wt_norm0)
          boot_res = add_info(seed=seed, res=boot_res, info=boot_res_info, ogd=og, pgd=pg, ogc=og, pgc=pg)
        }
      }
    }
  }
  return(list(point=point_res, sand_var1=sand_var_res1, sand_var2=sand_var_res2, 
              sand_var3=sand_var_res3, sand_var4=sand_var_res4, 
              sand_var5=sand_var_res5, sand_var6=sand_var_res6, boot=boot_res))
}

