# example simulation
# does not use any project-specific modules yet

library(epimodelcfa)
# load networks
folder_name <- "latest"
nets <- readRDS(here::here("input", "network_fits", folder_name, "nw.rds"))

# Ensure the 'local_tests' directory exists
local_tests_dir <- here::here("local_tests")
if (!dir.exists(local_tests_dir)) {
  dir.create(local_tests_dir, recursive = TRUE)
}


ncores <- max(1, parallel::detectCores() - 2L) # Reserve two cores by default
nsims <- ncores
years <- 1
nsteps <- 365 * years

# specify params & simulation controls
params <- EpiModel::param.net(
  units_per_year = 365,
  inf.prob = 0,
  rec.rate = 1 / 60,
  act.rate = 1,
  # will add infProbMTF = 1,
  # will add infProbFTM = 1,
  exitAge = 50,
  arrivalType = "departures",
  entryAge = 15,
  entryRaceNames = c("B", "H", "O", "W"),
  entryRaceProbs = c(0.12, 0.19, 0.11, 0.58),
  entryFemaleProb = 0.5,
)

inits <- EpiModel::init.net(i.num = 0)

controls <- EpiModel::control.net(
  nsims = nsims, nsteps = nsteps, ncores = nsims,
  verbose = FALSE,
  save.nwstats = FALSE,
  tergmLite = TRUE,
  resimulate.network = TRUE,
  dat.updates = epimodelcfa::resimnet_updates_sti,
  initialize.FUN = EpiModel::initialize.net,
  resim_nets.FUN = EpiModel::resim_nets,
  summary_nets.FUN = EpiModel::summary_nets,
  departures.FUN = epimodelcfa::mod_departures,
  arrivals.FUN = epimodelcfa::mod_arrivals,
  nwupdate.FUN = EpiModel::nwupdate.net,
  prevalence.FUN = EpiModel::prevalence.net,
  aging.FUN = epimodelcfa::mod_aging,
  verbose.FUN = EpiModel::verbose.net,
  epi.by = c("female"),
  tergmLite.track.duration = TRUE,
  save.run = TRUE,
  set.control.ergm = control.simulate.formula(
    MCMC.prop = ~sparse,
    MCMC.burnin = 2e+05
  )
)

# simulate
t1 <- Sys.time()
sim <- netsim(nets, params, inits, controls)
dur1 <- Sys.time() - t1
sim$simdur <- dur1



# Save the simulation object to a file
saveRDS(sim, file = file.path(local_tests_dir, paste0("sim_", years, "_years.rds")))

# Plotting and summarizing the simulation results (relationship stats)
plot_edges_history(sim, "main", "percent")
plot_edges_history(sim, "casual", "percent")

yaml_params_loc <- here::here("input", "params", "nw_params.yaml")
plot_final_degrees(sim, "main", yaml_params_loc)
plot_final_degrees(sim, "casual", yaml_params_loc)

get_mean_durations(sim, yaml_params_loc = yaml_params_loc)
