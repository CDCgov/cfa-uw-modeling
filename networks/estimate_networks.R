library(EpiModel)
devtools::load_all("epimodel-sti")
library(ggplot2)

# Required objects
this_seed <- 12345
ncores <- parallel::detectCores() - 1L
estimate_type <- "predicted" # alt: empirical
x <- yaml::read_yaml(here::here("params", paste0("nw_params_", estimate_type, ".yaml")))
main_mixmat_race <- readRDS(here::here("params", paste0("main_mixmat_", estimate_type, ".rds")))
cas_mixmat_race <- readRDS(here::here("params", paste0("casual_mixmat_", estimate_type, ".rds")))
main_mixmat_ag <- readRDS(here::here("params", paste0("main_mixmat_ag_", estimate_type, ".rds")))
cas_mixmat_ag <- readRDS(here::here("params", paste0("casual_mixmat_ag_", estimate_type, ".rds")))
nw <- generate_init_network(x, seed = this_seed, assign_deg_casual = TRUE)

# Estimate Networks ----------------------------------------------
## Main ---------------------------------------------------------
main_form <- ~ edges +
  nodematch(~race, diff = TRUE, levels = c(1:2, 4)) +
  nodefactor(~race, levels = c(1:2, 4)) +
  nodecov(~age) +
  nodecov(~agesq) +
  absdiff(~ sqrt(age_adj)) +
  nodefactor(~ deg_casual > 0) +
  offset(nodefactor(~ floor(age), levels = 1:2))

## --- main targets ------------------------
main_edges <- calc_targets(nw, x, "main", "edges")
main_nm_race <- calc_targets(nw, x, "main", "nodematch", "race")
main_nf_race <- calc_targets(nw, x, "main", "nodefactor", "race")
main_nc_age <- calc_targets(nw, x, "main", "nodecov", "age")
main_nc_agesq <- calc_targets(nw, x, "main", "nodecov", "age", attr_squared = TRUE)
main_agemix <- calc_targets(nw, x, "main", "absdiff_sqrt_age")
main_cross <- calc_targets(nw, x, "main", "cross_network")

main_targets <- unname(c(
  main_edges,
  main_nm_race[c(1:2, 4)],
  main_nf_race[c(1:2, 4)],
  main_nc_age,
  main_nc_agesq,
  main_agemix,
  main_cross
))

main_offset <- c(rep(-Inf, 2))

main_diss <- dissolution_coefs(
  dissolution = ~ offset(edges),
  duration = x$main$duration$duration$overall,
  d.rate = 0
)

# no same-sex ties, no concurrency in main network
main_constraints <- ~
  bd(maxout = 1) +
    sparse +
    blocks(attr = ~female, levels2 = diag(TRUE, 2)) +
    strat(attr = ~race, pmat = main_mixmat_race)

main_netest <- EpiModel::netest(
  nw = nw,
  formation = main_form,
  target.stats = main_targets,
  coef.form = main_offset,
  coef.diss = main_diss,
  constraints = main_constraints,
  set.control.ergm = control.ergm(
    MCMC.interval = 2048 * 10,
    MCMC.samplesize = 2048 * 10,
    SA.interval = 2048 * 10,
    SA.samplesize = 2048 * 10,
    main.method = "Stochastic-Approximation"
  )
)

main_dynamic <- netdx(
  main_netest,
  dynamic = TRUE, nsims = 10, nsteps = 1000,
  ncores = ncores
)
main_dynamic

if (!dir.exists(here::here("networks", "fits", Sys.Date()))) {
  dir.create(here::here("networks", "fits", Sys.Date()))
}
saveRDS(main_netest, here::here("networks", "fits", Sys.Date(), "main.rds"))
saveRDS(nw, here::here("networks", "fits", Sys.Date(), "nw.rds"))

nw <- set.vertex.attribute(nw, "deg_main", get_degree(main_netest$newnetwork))

## Casual ----------------------------------------------------

cas_form <- ~ edges +
  nodefactor(~age_group, levels = 1:6) +
  nodefactor(~race, levels = c(1:2, 4)) +
  nodematch(~race, diff = TRUE, levels = c(1, 3:4)) +
  concurrent() +
  nodefactor(~deg_main)

## -- casual targets -------------------------------------
cas_edges <- calc_targets(nw, x, "casual", "edges")
cas_nf_agegp <- calc_targets(nw, x, "casual", "nodefactor", "age_group")
cas_nf_race <- calc_targets(nw, x, "casual", "nodefactor", "race")
cas_nm_race <- calc_targets(nw, x, "casual", "nodematch", "race")
cas_conc <- calc_targets(nw, x, "casual", "concurrent")
cas_cross <- calc_targets(nw, x, "casual", "cross_network")

cas_targets <- unname(c(
  cas_edges,
  cas_nf_agegp[1:6],
  cas_nf_race[c(1:2, 4)],
  cas_nm_race[c(1, 3:4)],
  cas_conc,
  cas_cross
))

cas_offset <- c()

cdurs <- c(x$casual$duration$duration$agegrp_match[8], x$casual$duration$duration$agegrp_match[c(1:3, 5:7)])
cas_diss <- dissolution_coefs(
  dissolution = ~ offset(edges) + offset(nodematch(~age_group, diff = TRUE, levels = c(1:3, 5:7))),
  duration = cdurs,
  d.rate = 0
)

cas_diss <- dissolution_coefs(
  dissolution = ~ offset(edges),
  duration = x$casual$duration$duration$overall,
  d.rate = 0
)

cas_constraints <- ~
  sparse +
    blocks(attr = ~female, levels2 = diag(TRUE, 2)) +
    strat(attr = ~race, pmat = cas_mixmat_race) +
    strat(attr = ~age_group, pmat = cas_mixmat_ag)

cas_netest <- EpiModel::netest(
  nw = nw,
  formation = cas_form,
  target.stats = cas_targets,
  coef.form = cas_offset,
  coef.diss = cas_diss,
  constraints = cas_constraints,
  set.control.ergm = control.ergm(
    MCMC.interval = 2048 * 10,
    MCMC.samplesize = 2048 * 10,
    SA.interval = 2048 * 10,
    SA.samplesize = 2048 * 10,
    main.method = "Stochastic-Approximation"
  )
)

cas_dynamic <- netdx(
  cas_netest,
  dynamic = TRUE, nsims = 10, nsteps = 1000, ncores = ncores
)
cas_dynamic

saveRDS(cas_netest, here::here("networks", "fits", Sys.Date(), "casual.rds"))

nw <- set.vertex.attribute(nw, "deg_casual", get_degree(cas_netest$newnetwork))


# Inst network ----------------------------------------------------
inst_form <- ~ edges +
  nodefactor(~race, levels = c(1:2, 4)) +
  nodefactor(~age_group, levels = 1:3) +
  absdiff(~ sqrt(age_adj))

inst_edges <- calc_targets(nw, x, "inst", "edges", inst_correct = TRUE, time_unit = "days")
inst_nf_race <- calc_targets(nw, x, "inst", "nodefactor", "race", inst_correct = TRUE, time_unit = "days")
inst_nf_agegrp <- calc_targets(nw, x, "inst", "nodefactor", "age_group",
  inst_correct = TRUE, time_unit = "days"
)
inst_agemix <- calc_targets(nw, x, "inst", "absdiff_sqrt_age", inst_correct = TRUE, time_unit = "days")

inst_targets <- unname(c(
  inst_edges,
  inst_nf_race[c(1:2, 4)],
  inst_nf_agegrp[1:3],
  inst_agemix
))

inst_offset <- c()

inst_diss <- dissolution_coefs(
  dissolution = ~ offset(edges),
  duration = 1,
  d.rate = 0
)

inst_constraints <- ~
  sparse +
    blocks(attr = ~female, levels2 = diag(TRUE, 2))


inst_netest <- EpiModel::netest(
  nw = nw,
  formation = inst_form,
  target.stats = inst_targets,
  coef.form = inst_offset,
  coef.diss = inst_diss,
  constraints = inst_constraints,
  set.control.ergm = control.ergm(
    MCMC.interval = 2048 * 10,
    MCMC.samplesize = 2048 * 10,
    SA.interval = 2048 * 10,
    SA.samplesize = 2048 * 10,
    main.method = "Stochastic-Approximation"
  )
)

inst_static <- netdx(inst_netest, dynamic = FALSE, nsims = 500)
inst_static

saveRDS(inst_netest, here::here("networks", "fits", Sys.Date(), "inst.rds"))

##### Save out
est <- list(main_netest, cas_netest, inst_netest)
save(est, file = here::here("networks", "fits", Sys.Date(), "nw.rda"))

#### Plot compare to empirical dists --------------------------------------------

## empirical dist
emp <- data.frame(
  main = x$main$nodefactor$age_race,
  cas = x$casual$nodefactor$age_race,
  age = rep(15:49, each = 4),
  race = rep(c("B", "H", "O", "W"), 7),
  names = paste0(rep(15:49, each = 4), c("B", "H", "O", "W"))
)

# modeled dist
dat <- data.frame(
  age = floor(nw %v% "age"),
  age_group = nw %v% "age_group",
  race = nw %v% "race",
  d = nw %v% "deg_main",
  c = nw %v% "deg_casual"
)


comp <- dat |>
  dplyr::group_by(age, race) |>
  dplyr::summarize(
    main_mod = mean(d),
    cas_mod = mean(c)
  ) |>
  dplyr::left_join(emp, by = c("age", "race")) |>
  tidyr::pivot_longer(
    cols = c("main", "main_mod", "cas", "cas_mod"),
    names_to = c("type"), values_to = "vals"
  )

comp |>
  dplyr::filter(type %in% c("main", "main_mod")) |>
  ggplot2::ggplot(ggplot2::aes(x = age, y = vals, col = type)) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(span = 0.75) +
  ggplot2::facet_wrap(~race)

comp |>
  dplyr::filter(type %in% c("cas", "cas_mod")) |>
  ggplot2::ggplot(ggplot2::aes(x = age, y = vals, col = type)) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(span = 0.75) +
  ggplot2::facet_wrap(~race)
