# install.packages(c("dplyr", "locpol", "WeightIt", "ggplot2"))
source("call_sims.R")
library(ggplot2)

# Get true parameters
ground_truth_info = get_true_params()

# Run simulation
sims_point_results = collect_point_sims(start_sim=6, end_sim=6, n=1000, 
                                        testD=ground_truth_info$D.vals, wt_normD=T, wt_norm0=T)
sims_point_results$psiD 
# muD=Outcome Regression Model for Dose is correct
# piD=Propensity score model for Dose is correct
# mu0=Outcome Regression Model for Control Treatment is correct
# piA=Binary Propensity score model for Treatment is correct

# Create Figure 2 with simulation results
plot_data = data.frame(Dose=rep(ground_truth_info$D.vals,4), 
                       outcome=c(unlist(ground_truth_info$trt_effect), 
                                 unlist(sims_point_results$psiD[1, paste0("dr",1:50)]),
                                 unlist(sims_point_results$psiD[1, paste0("naive",1:50)]),
                                 unlist(sims_point_results$psiD[1, paste0("twfe",1:50)])),
                       Estimator=rep(c("A", "B", "C", "D"), each=length(ground_truth_info$D.vals)))
p = plot_data %>% ggplot(aes(x=Dose,y=outcome, group=Estimator)) + 
  geom_line(aes(linetype=Estimator), size=1.25) +
  theme_minimal() + 
  scale_linetype_manual(values=c("solid","dashed","11","dotdash"), 
                        labels=c("Truth", "MR", "Naive", "TWFE")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),
        axis.text = element_text(size=12),
        legend.title = element_text(size=14), legend.text = element_text(size=12),
        legend.position = c(0.15, 0.15)) + guides(linetype=guide_legend(keywidth=3)) +
  xlab(expression(delta)) + ylab(expression(Psi(delta)))
print(p)


# Run simulation analyses with confidence intervals (time-consuming)
sims_ci_results = run_ci_sims(start_sim=6, end_sim=6, n=1000,
                              testD=ground_truth_info$D.vals, nboots=10, # best to use nboots>200
                              get_vars=T, wt_normD=T, wt_norm0=T)

# Variances using sandwich variance approach
sims_ci_results$sand_var$psiD

# SD's using Bootstrap variance approach
sims_ci_results$boot$psiD[,c("Seed", "muD", "piD", "mu0", "piA", paste0("dr", 1:length(ground_truth_info$D.vals), ".sd"))]

# Lower bounds using Bootstrap percentile approach
sims_ci_results$boot$psiD[,c("Seed", "muD", "piD", "mu0", "piA", paste0("dr", 1:length(ground_truth_info$D.vals), ".lower"))]

# Upper bounds using Bootstrap percentile approach
sims_ci_results$boot$psiD[,c("Seed", "muD", "piD", "mu0", "piA", paste0("dr", 1:length(ground_truth_info$D.vals), ".upper"))]
