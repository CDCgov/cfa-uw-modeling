# This script simulates a network with vital dynamics and interdependent networks.
# No infection dynamics
# Goal: evaluate where networks find equilibrium

library(epimodelcfa)
# load networks
# adj b/c people aging out more likely to be in a main rel
main_drate_adjustment <- 1.4
drate_labels <- (main_drate_adjustment - 1) * 100 # labels for main_drate_adjustment for file saves

# Ensure the 'local_tests' directory exists
local_tests_dir <- here::here("local_tests")
if (!dir.exists(local_tests_dir)) {
  dir.create(local_tests_dir, recursive = TRUE)
}

nets <- readRDS(here::here("input", "network_fits", "latest", "nw.rds"))

ncores <- 20
nsims <- ncores
years <- 25
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
sim_name <- paste0("sim_burnin_", drate_labels[i], ".rds")
saveRDS(sim, file = file.path(local_tests_dir, sim_name))


# Plotting and summarizing the simulation results (relationship stats)
# Load the simulation object (unless continuing from above)
sim <- readRDS(file.path(local_tests_dir, sim_name))

# Make plots
p1 <- plot_edges_history(sim, "main", "percent")
p2 <- plot_edges_history(sim, "casual", "percent")
yaml_params_loc <- here::here("input", "params", "nw_params.yaml")
p3 <- plot_final_degrees(sim, "main", yaml_params_loc)
p4 <- plot_final_degrees(sim, "casual", yaml_params_loc)

# Save the plots to a PDF file, one per drate adjustment
p <- list(p1, p2, p3, p4)
pdf(file.path(local_tests_dir, paste0("sim_burnin_", drate_labels[i], ".pdf")), width = 10, height = 6)
p
dev.off()
rm(sim, p)


get_mean_durations(sim, yaml_params_loc = yaml_params_loc)
