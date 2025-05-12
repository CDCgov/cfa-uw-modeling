# example simulation
# does not use any project-specific modules yet

library(EpiModel)
devtools::load_all("epimodel-sti")
# load networks
data(nw_50000)
data(mixmats)

main_mixmat_race <- mixmats[[1]]
cas_mixmat_race <- mixmats[[3]]
cas_mixmat_ag <- mixmats[[4]]

nsims <- 10
nsteps <- 365 * 10
# specify params & simulation controls
params <- param.net(
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

ergm_control_main <- control.simulate.formula.ergm(
  MCMC.prop = ~ sparse + strat(attr = ~race, pmat = main_mixmat_race),
  MCMC.interval = 2048 * 10
)

tergm_control_main <- control.simulate.formula.tergm(
  MCMC.prop = ~ sparse + strat(attr = ~race, pmat = main_mixmat_race),
  MCMC.maxchanges = Inf
)

ergm_control_casual <- control.simulate.formula.ergm(
  MCMC.prop = ~ sparse +
    strat(attr = ~race, pmat = cas_mixmat_race) +
    strat(attr = ~age_group, pmat = cas_mixmat_ag),
  MCMC.interval = 2048 * 10
)

tergm_control_casual <- control.simulate.formula.tergm(
  MCMC.prop = ~ sparse +
    strat(attr = ~race, pmat = cas_mixmat_race) +
    strat(attr = ~age_group, pmat = cas_mixmat_ag),
  MCMC.maxchanges = Inf
)

ergm_control <- list(ergm_control_main, ergm_control_casual, ergm_control_casual)
tergm_control <- list(tergm_control_main, tergm_control_casual, tergm_control_casual)
class(ergm_control) <- "multilayer"
class(tergm_control) <- "multilayer"

controls <- control.net(
  nsims = nsims, nsteps = nsteps, ncores = nsims,
  verbose = FALSE,
  save.nwstats = TRUE, nwstats.formula = "formation",
  tergmLite = TRUE, resimulate.network = TRUE,
  dat.updates = resimnet_updates_sti,
  initialize.FUN = initialize.net,
  resim_nets.FUN = resim_nets,
  summary_nets.FUN = summary_nets,
  infection.FUN = infection.net,
  recovery.FUN = recovery.net,
  departures.FUN = mod_departures,
  arrivals.FUN = mod_arrivals,
  nwupdate.FUN = nwupdate.net,
  prevalence.FUN = prevalence.net,
  aging.FUN = mod_aging,
  verbose.FUN = verbose.net,
  epi.by = c("female"),
  set.control.ergm = ergm_control,
  set.control.tergm = tergm_control,
  tergmLite.track.duration = TRUE,
  save.run = TRUE,
  cumulative.edgelist = TRUE
)

# simulate
t1 <- Sys.time()
sim <- netsim(nw_50000, params, inits, controls)
dur1 <- Sys.time() - t1

saveRDS(sim, file = here::here("localtests", "sim_10yrs_mcmcupdates.rds"))
# plot prevalence and daily new infecteds

plot(sim, y = c("a.flow", "d.flow", "meanAge"))


nws <- sim$stats$nwstats
mfull <- NULL
for (i in seq_len(nsims)) {
  m <- nws[[i]][[1]]
  m$day <- seq_len(nsteps)
  m$sim <- i
  mfull <- rbind(mfull, m)
}

cfull <- NULL
for (i in seq_len(nsims)) {
  c <- nws[[i]][[2]]
  c$day <- seq_len(nsteps)
  c$sim <- i
  cfull <- rbind(cfull, c)
}

ifull <- NULL
for (i in seq_len(nsims)) {
  inst <- nws[[i]][[3]]
  inst$day <- seq_len(nsteps)
  inst$sim <- inst
  ifull <- rbind(ifull, inst)
}

library(ggplot2)
mfull |>
  tidyr::pivot_longer(-c("day", "sim"), names_to = "ergm_term", values_to = "val") |>
  ggplot(aes(x = day, y = val, group = sim)) +
  geom_smooth() +
  facet_wrap(~ergm_term, scales = "free")

cfull |>
  tidyr::pivot_longer(-c("day", "sim"), names_to = "ergm_term", values_to = "val") |>
  ggplot(aes(x = day, y = val, group = sim)) +
  geom_smooth() +
  facet_wrap(~ergm_term, scales = "free")

ifull |>
  tidyr::pivot_longer(-c("day", "sim"), names_to = "ergm_term", values_to = "val") |>
  ggplot(aes(x = day, y = val)) +
  geom_smooth() +
  facet_wrap(~ergm_term, scales = "free")

# checking age dist
age <- floor(sim$run$sim1$attr$age)
deg_main <- sim$run$sim1$attr$deg_main
main_el <- as.vector(sim$run$sim1$el[[1]])
main_ages <- age[main_el]
main_activity <- table(main_ages) / table(age)

cas_el <- as.vector(sim$run$sim1$el[[2]])
cas_ages <- age[cas_el]
cas_activity <- table(cas_ages) / table(age)

# data expectations
x <- yaml::read_yaml(here::here("params", paste0("nw_params_", "predicted", ".yaml")))
dat <- data.frame(
  main = x$main$nodefactor$age_race,
  casual = x$casual$nodefactor$age_race,
  age = rep(15:49, each = 4),
  race = rep(c("B", "H", "O", "W"), times = 35)
)

datsum <- dat |>
  dplyr::group_by(age) |>
  dplyr::summarize(
    main = mean(main, vartype = NULL),
    casual = mean(casual, vartype = NULL)
  )

datmain <- cbind(datsum[c("age", "main")], main_activity)
datcas <- cbind(datsum[c("age", "casual")], cas_activity)

plot(x = datmain$age, y = datmain$main, col = "blue")
points(x = datmain$age, y = datmain$Freq, col = "red")

plot(x = datcas$age, y = datcas$casual, col = "blue")
points(x = datcas$age, y = datcas$Freq, col = "red")

# mean rel durs at end
cuml.el <- sim$run$sim1$el_cuml_cur
m.el <- cuml.el[[1]]
m.el$dur <- nsteps - m.el$start
summary(m.el$dur)
x$main$duration$duration

c.el <- cuml.el[[2]]
c.el$dur <- nsteps - c.el$start
summary(c.el$dur)
x$casual$duration$duration
