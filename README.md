# DiD-continuous
This repository implements code from research on multiply robust difference-in-differences methodology for continuous exposures (https://doi.org/10.1093/biomtc/ujaf015).

## Example

An example simulated dataset is provided in `example_simulated_data.RData`. 

An example script calling a simulation and running an analysis for the Average Dose Effect on the Treated (ADT) is provided in `example_run.R`. An example figure is also produced in the script. To calculate pointwise confidence intervals, the `run_ci_sims()` function can be run with `nboots` parameter (bootstrap) set to greater than zero and/or the `get_vars` flag set to TRUE (sandwich). 

Results are given across model types mentioned in the manuscript and supplement (Two-Way Fixed Effects (twfe), Confounder-Naive (naive), Outcome Regression (or), Inverse Probability Weighted (IPW), and Multiply Robust (DR)) for all testing dose values (e.g., 1:50) and across each setting of nuisance model specification (muD, piD, mu0, and piA).

Additional results are represented for a parametric regression on dose-specific pseudo-outcomes when the regression function is correctly specified (dr_paramc) and incorrectly specified (dr_parami).

Results are given for the causal effect curve (psiD) and decomposed into the components related to trends under a given dose among the treated (thetaD) and under the control exposure (theta0).

## Functions

Relevant backend functions to fit outcome models, propensity score models, and nonparametric kernel regressions are found in `dose_component_functions.R`, `ctl_component_functions.R`, and `kernel_functions.R`. Backend functions for variance and confidence interval calculations are found in `sandwich_variance_functions.R` and `bootstrap_functions.R`. Functions to generate the simulations and to call multiple simulations are found in `sim_generation.R` and `call_sims.R`, respectively.
