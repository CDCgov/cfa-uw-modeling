library(EpiModel)
setwd("~/cfa-uw-modeling")
devtools::load_all("epimodel-sti")

# Required objects
this_seed <- 12345
ncores <- parallel::detectCores() - 1L
estimate_type <- "predicted"
x <- yaml::read_yaml(here::here("params", paste0("nw_params_", estimate_type, ".yaml")))
x$pop$size <- 50000
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

##### Save out
nw_50000 <- list(main_netest, cas_netest, inst_netest)
mixmats <- list(main_mixmat_race, main_mixmat_ag, cas_mixmat_race, cas_mixmat_ag)
setwd("epimodel-sti")
usethis::use_data(nw_50000, mixmats, overwrite = TRUE)
