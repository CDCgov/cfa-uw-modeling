# example simulation
# does not use any project-specific modules yet

library(EpiModel)
devtools::load_all("epimodel-sti")
# load networks
data(nw_50000)
data(mixmats)

main_mixmat_race <- mixmats[[1]]
main_mixmat_ag <- mixmats[[2]]
cas_mixmat_race <- mixmats[[3]]
cas_mixmat_ag <- mixmats[[4]]

nsims <- 1
nsteps <- 365
# specify params & simulation controls
params <- param.net(
  inf.prob = 0,
  rec.rate = 1 / 60,
  act.rate = 1,
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
  verbose.FUN = verbose.net
)

# simulate
t1 <- Sys.time()
sim <- netsim(nw_50000, params, inits, controls)
dur1 <- Sys.time() - t1

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
  dplyr::select(-c("offset(nodefactor.floor(age).15)", "offset(nodefactor.floor(age).16)")) |>
  tidyr::pivot_longer(-c("day", "sim"), names_to = "ergm_term", values_to = "val") |>
  ggplot(aes(x = day, y = val)) +
  geom_smooth() +
  facet_wrap(~ergm_term, scales = "free")

cfull |>
  tidyr::pivot_longer(-c("day", "sim"), names_to = "ergm_term", values_to = "val") |>
  ggplot(aes(x = day, y = val)) +
  geom_smooth() +
  facet_wrap(~ergm_term, scales = "free")

ifull |>
  tidyr::pivot_longer(-c("day", "sim"), names_to = "ergm_term", values_to = "val") |>
  ggplot(aes(x = day, y = val)) +
  geom_smooth() +
  facet_wrap(~ergm_term, scales = "free")
