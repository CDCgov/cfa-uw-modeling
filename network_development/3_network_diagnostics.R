library(epimodelcfa)
devtools::load_all("~/repos/cfa-epimodel")
folder_name <- "latest"
nets <- readRDS(here::here("input", "network_fits", folder_name, "nw.rds"))
yaml_params_loc <- here::here("input", "params", "nw_params.yaml")

# MCMC Diagnostics & GOF ---------------------------------------------------
## note: for some reason, the plots look worse when using SA than using MCMLE, but other metrics ok
## note: GOF metrics will be slightly different every time it is run, but should be similar to each other
## note: FYI GOF also takes a while to run

m1 <- mcmc.diagnostics(nets[[1]]$fit)
plot(m1)
g1 <- gof(nets[[1]]$fit, GOF = ~model)

m2 <- mcmc.diagnostics(nets[[2]]$fit)
plot(m2)
g2 <- gof(nets[[2]]$fit, GOF = ~model)

m3 <- mcmc.diagnostics(nets[[3]]$fit)
plot(m3)
g3 <- gof(nets[[3]]$fit, GOF = ~model)

# Network Diagnostics via Simulation --------------------------------------------------------
## Shared Parameters

ncores <- 5 #parallel::detectCores() - 2L
nsims <- ncores
nsteps <- 1000
main_casual_ergm_control <- control.simulate.formula(MCMC.prop = ~sparse, MCMC.burnin = 2e+05)
main_casual_tergm_control <- control.simulate.formula.tergm(MCMC.prop = ~ discord + sparse)

## Main network diagnostics
main_dynamic <- netdx(
  nets[[1]],
  dynamic = TRUE, nsims = nsims, nsteps = nsteps, ncores = ncores,
  set.control.ergm = main_casual_ergm_control,
  set.control.tergm = main_casual_tergm_control,
  keep.tedgelist = TRUE
)
main_dynamic
plot(main_dynamic)

## Casual network diagnostics
cas_dynamic <- netdx(
  nets[[2]],
  dynamic = TRUE, nsims = nsims, nsteps = nsteps, ncores = ncores,
  set.control.ergm = main_casual_ergm_control,
  set.control.tergm = main_casual_tergm_control,
  keep.tedgelist = TRUE
)
cas_dynamic
plot(cas_dynamic)

## Instantaneous network diagnostics
inst_static <- netdx(nets[[3]], dynamic = FALSE, nsims = 500)
inst_static

# Plot Degree Distributions --------------------------------------------
## Note: Because the ERGM terms used as target statistics in the model fit do not discretely represent each
## age and race degree combination, the degree distributions will not match the empirical distributions perfectly,
## but should be similar in shape.
##
## Note 2: The networks within the netest and netdx objects are independent of each other, and the netdx
## objects are simulated in a closed population. Thus, we will need to run a netsim simulation to see how
## the degree distributions change over time when simulated with vital dynamics and changing nodal degree terms
## that influence formation in other networks.

## Degree based on netest object fit
plot_final_degrees(nets[[1]], "main", yaml_params_loc)
plot_final_degrees(nets[[2]], "casual", yaml_params_loc)

## Degree based on netdx simulations
plot_final_degrees(main_dynamic, "main", yaml_params_loc)
plot_final_degrees(cas_dynamic, "casual", yaml_params_loc)
