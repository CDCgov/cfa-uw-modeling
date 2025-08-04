library(epimodelcfa)

folder_name <- "latest"
nets <- readRDS(here::here("networks", "fits", folder_name, "nw.rds"))

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

ncores <- parallel::detectCores() - 1L
nsims <- ncores
nsteps <- 1000
main_casual_ergm_control <- control.simulate.formula(MCMC.prop = ~sparse, MCMC.burnin = 2e+05)
main_casual_tergm_control <- control.simulate.formula.tergm(MCMC.prop = ~ discord + sparse)

## Main network diagnostics
main_dynamic <- netdx(
  nets[[1]],
  dynamic = TRUE, nsims = nsims, nsteps = nsteps, ncores = ncores,
  set.control.ergm = main_casual_ergm_control,
  set.control.tergm = main_casual_tergm_control
)
main_dynamic
plot(main_dynamic)

## Casual network diagnostics
cas_dynamic <- netdx(
  nets[[2]],
  dynamic = TRUE, nsims = nsims, nsteps = nsteps, ncores = ncores,
  set.control.ergm = main_casual_ergm_control,
  set.control.tergm = main_casual_tergm_control
)
cas_dynamic
plot(cas_dynamic)

## Instantaneous network diagnostics
inst_static <- netdx(nets[[3]], dynamic = FALSE, nsims = 500)
inst_static

# Plot compare to empirical dists --------------------------------------------

## Target Parameters
x <- yaml::read_yaml(here::here("networks", "params", "nw_params.yaml"))
params <- data.frame(
  main_targets = x$main$nodefactor$age_race,
  casual_targets = x$casual$nodefactor$age_race,
  age = rep(15:49, 4),
  race = rep(c("B", "H", "O", "W"), each = length(15:49)),
  names = paste0(rep(15:49, 4), rep(c("B", "H", "O", "W"), each = length(15:49)))
)

## Simulated Degrees extracted from network fits
### get info from inst network ()
dat <- data.frame(
  age = floor(nets[[1]]$newnetwork %v% "age"),
  age_group = nets[[1]]$newnetwork %v% "age_group",
  race = nets[[1]]$newnetwork %v% "race",
  d = get_degree(nets[[1]]$newnetwork),
  c = get_degree(nets[[2]]$newnetwork)
)

## Combine and summarize
comp <- dat |>
  dplyr::group_by(age, race) |>
  dplyr::summarize(
    main_fit = mean(d),
    casual_fit = mean(c)
  ) |>
  dplyr::left_join(params, by = c("age", "race")) |>
  tidyr::pivot_longer(
    cols = c("main_targets", "main_fit", "casual_targets", "casual_fit"),
    names_to = c("type"), values_to = "vals"
  )

## Plot main
comp |>
  dplyr::filter(type %in% c("main_targets", "main_fit")) |>
  ggplot2::ggplot(ggplot2::aes(x = age, y = vals, col = type)) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(span = 0.75) +
  ggplot2::facet_wrap(~race)

## Plot casual
comp |>
  dplyr::filter(type %in% c("casual_targets", "casual_fit")) |>
  ggplot2::ggplot(ggplot2::aes(x = age, y = vals, col = type)) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(span = 0.75) +
  ggplot2::facet_wrap(~race)
