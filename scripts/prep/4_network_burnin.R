# example simulation
# does not use any project-specific modules yet

library(EpiModel)
devtools::load_all("epimodel-sti")
# load networks
folder_name <- "latest"
nets <- readRDS(here::here("networks", "fits", folder_name, "nw.rds"))

ncores <- parallel::detectCores() - 1L
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
  entryFemaleAgeAdj = 1.8331817
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

get_edge_stats <- function(sim, attrs = c("edges_main", "edges_casual")) {
  # extract edges target for given network
  target_main <- sim$nwparam[[1]]$target.stats[1]
  target_casual <- sim$nwparam[[2]]$target.stats[1]

  sim |>
    as.data.frame() |>
    dplyr::select(time, sim, dplyr::all_of(attrs)) |>
    dplyr::mutate(
      main_diff = edges_main - target_main,
      casual_diff = edges_casual - target_casual,
      main_diff_perc = (main_diff) / target_main * 100,
      casual_diff_perc = (casual_diff) / target_casual * 100
    ) |>
    dplyr::group_by(time) |>
    dplyr::mutate(
      mean_main_diff_perc = mean(main_diff_perc, na.rm = TRUE),
      mean_casual_diff_perc = mean(casual_diff_perc, na.rm = TRUE)
    )
}

e <- get_edge_stats(sim)

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

# frequency of rels by age in networks at end of simulation
summarize_final_degrees <- function(sim) {
  simdat <- NULL

  for (i in seq_len(nsims)) {
    this_sim <- paste0("sim", i)
    d <- data.frame(
      age = floor(sim[["network"]][[this_sim]][[1]] %v% "age"),
      race = sim[["network"]][[this_sim]][[1]] %v% "race",
      deg_main = get_degree(sim[["network"]][[this_sim]][[1]]),
      deg_cas = get_degree(sim[["network"]][[this_sim]][[2]]),
      sim = this_sim
    )
    simdat <- rbind(simdat, d)
  }

  sim_summary <- simdat |>
    dplyr::group_by(sim, age, race) |>
    dplyr::summarize( # calc mean degree by sim, age, race
      main = mean(deg_main),
      casual = mean(deg_cas),
      .groups = "drop"
    ) |>
    tidyr::pivot_longer(
      cols = c("main", "casual"),
      names_to = "type",
      values_to = "deg"
    ) |>
    dplyr::group_by(age, race, type) |>
    dplyr::summarize( # mean and sd of mean degree across sims
      degree = mean(deg),
      IQR1 = quantile(deg, 1 / 4),
      IQR3 = quantile(deg, 3 / 4),
      .groups = "drop"
    ) |>
    dplyr::mutate(data = "simulated")
}

get_target_degrees <- function(target_yaml_file) {
  x <- yaml::read_yaml(target_yaml_file)
  dat <- data.frame(
    main = x$main$nodefactor$age_race,
    casual = x$casual$nodefactor$age_race,
    age = rep(15:49, 4),
    race = rep(c("B", "H", "O", "W"), each = 35)
  )

  dat |>
    dplyr::group_by(age, race) |>
    dplyr::summarize(
      main = mean(main),
      casual = mean(casual),
      .groups = "drop"
    ) |>
    dplyr::mutate(data = "targets") |>
    tidyr::pivot_longer(
      cols = c("main", "casual"),
      names_to = "type",
      values_to = "degree"
    )
}
s <- summarize_final_degrees(sim)
t <- get_target_degrees(here::here("networks", "params", "nw_params.yaml"))

# targets from yaml file

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

# mean rel durs at end (may not match targets if simulation is not long enough)
m <- sim$network$sim1[[1]] %n% "lasttoggle"
mdurs <- abs(m[, 3] - nsteps)
summary(mdurs)

c <- sim$network$sim4[[2]] %n% "lasttoggle"
cdurs <- abs(c[, 3] - nsteps)
summary(cdurs)
