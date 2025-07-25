# example simulation
# does not use any project-specific modules yet

library(EpiModel)
devtools::load_all("epimodel-sti")
# load networks
folder_name <- "latest"
nets <- readRDS(here::here("networks", "fits", folder_name, "nw.rds"))

ncores <- max(1, parallel::detectCores() - 2L) # Reserve two cores by default
nsims <- ncores
nsteps <- 365 * 30

# specify params & simulation controls
params <- param.net(
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

inits <- init.net(i.num = 0)

controls <- control.net(
  nsims = nsims, nsteps = nsteps, ncores = nsims,
  verbose = FALSE,
  save.nwstats = FALSE,
  tergmLite = TRUE,
  resimulate.network = TRUE,
  dat.updates = resimnet_updates_sti,
  initialize.FUN = initialize.net,
  resim_nets.FUN = resim_nets,
  summary_nets.FUN = summary_nets,
  departures.FUN = mod_departures,
  arrivals.FUN = mod_arrivals,
  nwupdate.FUN = nwupdate.net,
  prevalence.FUN = prevalence.net,
  aging.FUN = mod_aging,
  verbose.FUN = verbose.net,
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

# Ensure the 'localtests' directory exists
localtests_dir <- here::here("localtests")
if (!dir.exists(localtests_dir)) {
  dir.create(localtests_dir, recursive = TRUE)
}

# Save the simulation object to a file
saveRDS(sim, file = here::here("localtests", "sim_30yrs.rds"))

plot(sim, y = c("edges_main"), sim.lines = TRUE)
plot(sim, y = c("edges_casual"), sim.lines = TRUE)


e <- get_edges_history(sim)

e |>
  ggplot2::ggplot(ggplot2::aes(x = time, y = main_diff_perc, color = sim)) +
  ggplot2::geom_line() +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
  ggplot2::geom_line(ggplot2::aes(y = mean_main_diff_perc), color = "black", linewidth = 1) +
  ggplot2::labs(title = "Main Network Edges Over Time")

e |>
  ggplot2::ggplot(ggplot2::aes(x = time, y = casual_diff_perc, color = sim)) +
  ggplot2::geom_line() +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
  ggplot2::geom_line(ggplot2::aes(y = mean_casual_diff_perc), color = "black", linewidth = 1) +
  ggplot2::labs(title = "Casual Network Edges Over Time")

s <- summarize_final_degrees(sim)
t <- get_target_degrees_age_race(here::here("networks", "params", "nw_params.yaml"))


s |>
  dplyr::filter(type == "main") |>
  ggplot2::ggplot(ggplot2::aes(x = age, y = degree, color = data)) +
  ggplot2::geom_point() +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = IQR1, ymax = IQR3), width = 0.2) +
  ggplot2::geom_point(
    data = t |> dplyr::filter(type == "main"),
    ggplot2::aes(x = age, y = degree, color = data), linewidth = 3
  ) +
  ggplot2::facet_wrap(~race)

s |>
  dplyr::filter(type == "casual") |>
  ggplot2::ggplot(ggplot2::aes(x = age, y = degree, color = data)) +
  ggplot2::geom_point() +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = IQR1, ymax = IQR3), width = 0.2) +
  ggplot2::geom_point(
    data = t |> dplyr::filter(type == "casual"),
    ggplot2::aes(x = age, y = degree, color = data), linewidth = 3
  ) +
  ggplot2::facet_wrap(~race)



get_mean_durations(sim)
