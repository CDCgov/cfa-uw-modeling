# This script simulates a network with
# vital dynamics and interdependent networks.
# No infection dynamics

library(epimodelcfa)
# load networks

nets <- readRDS(here::here(
  "input",
  "network_fits",
  "nw.rds"
))

yaml_params_loc <- here::here(
  "input",
  "params",
  "nw_params.yaml"
)

ncores <- parallel::detectCores()
nsims <- parallel::detectCores()
years <- 5
nsteps <- 365 * years

# specify params & simulation controls
params <- EpiModel::param.net(
  units_per_year = 365,
  inf.prob = 0,
  rec.rate = 1 / 60,
  act.rate = 1,
  vital_dynamics = TRUE,
  exit_age = 50,
  entry_age = 15,
  arrivalType = "departures",
  entry_female_prob = 0.5,
  entry_race_names = c("B", "H", "O", "W"),
  entry_race_probs = c(0.12, 0.19, 0.11, 0.58)
)

inits <- EpiModel::init.net(i.num = 0)

controls <- EpiModel::control.net(
  nsims = nsims,
  nsteps = nsteps,
  ncores = nsims,
  verbose = FALSE,
  save.nwstats = TRUE,
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
  edge_tracker.FUN = epimodelcfa::track_edges,
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
sim_name <- paste0("sim_vitals_", nsteps, ".rds")
saveRDS(
  sim,
  file = file.path(
    "input",
    "network_fits",
    sim_name
  )
)

# Plotting and summarizing the simulation results (relationship stats)
# Load the simulation object (unless continuing from above)
sim <- readRDS(file.path("input", "network_fits", sim_name))

# Make plots
p1 <- plot_edges_history(sim, "main", "percent")
p2 <- plot_edges_history(sim, "casual", "percent")
p3 <- plot_edges_history(sim, "main", "absolute")
p4 <- plot_edges_history(sim, "casual", "absolute")
p5 <- plot_final_degrees(sim, "main", yaml_params_loc)
p6 <- plot_final_degrees(sim, "casual", yaml_params_loc)

# Save the plots to a PDF file, one per drate adjustment
p <- list(p1, p2, p3, p4, p5, p6)
pdf(
  file.path(
    "input",
    "network_fits",
    paste0("sim_vitals_", nsteps, ".pdf")
  ),
  width = 10,
  height = 6
)
print(p)
dev.off()
