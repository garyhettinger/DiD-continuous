# DiD-continuous
This repository implements code from research on multiply robust difference-in-differences methodology for continuous exposures (https://arxiv.org/abs/2401.14355).

## Example

An example simulated dataset is provided in `example_simulated_data.RData`. 

An example script calling a simulation and running an analysis for the Average Dose Effect on the Treated (ADT) is provided in `example_run.R`. An example figure is also produced in the script. To calculate pointwise confidence intervals, the `run_ci_sim()` function can be uncommented with `nboots` parameter (bootstrap) set to greater than zero and/or the `get_vars` flag set to TRUE (sandwich). 

## Functions

Relevant backend functions to fit outcome models, propensity score models, and nonparametric kernel regressions are found in `dose_component_functions.R`, `ctl_component_functions.R`, and `kernel_functions.R`. Backend functions for variance and confidence interval calculations are found in `variance_functions.R` and `bootstrap_functions.R`. Functions to generate the simulations and to call multiple simulations are found in `sim_generation.R` and `call_sims.R`, respectively.