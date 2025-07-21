library(EpiModel)
devtools::load_all("epimodel-sti")

# Required objects --------------------------------------------
this_seed <- 11111
dept_rate <- (1 / (50 - 15)) * (1 / 365)
x <- yaml::read_yaml(here::here("networks", "params", "nw_params.yaml"))

# Generate initial network ----------------------------------------
## generates empty network with nodal attributes reflecting pop specs
## and initial guess at casual degree distribution to fit main network
nw <- generate_init_network(
  x,
  seed = this_seed,
  assign_deg_casual = FALSE,
  assign_deg_main = TRUE,
  olderpartner = TRUE
)

# Estimate Networks ----------------------------------------------
## Shared control settings ---------------------------
ergm_controls <- control.ergm(
  seed = this_seed,
  MCMC.prop = ~sparse,
  MCMC.interval = 2048 * 10,
  MCMC.samplesize = 2048 * 10,
  SA.interval = 2048 * 10,
  SA.samplesize = 2048 * 10,
  main.method = "Stochastic-Approximation"
)
## Casual ----------------------------------------------------
### Formation model
cas_form <- ~ edges +
  nodematch(~race, diff = TRUE, levels = c(1:2, 4)) +
  nodefactor(~race, levels = c(1:2, 4)) +
  nodefactor(~age_group, levels = 1:5) +
  concurrent() +
  offset(nodefactor(~deg_main))

### Calc and Compile Targets
cas_edges <- calc_targets(nw, x, "casual", "edges")
cas_nf_agegp <- calc_targets(nw, x, "casual", "nodefactor", "age_group")
cas_nf_race <- calc_targets(nw, x, "casual", "nodefactor", "race")
cas_nm_race <- calc_targets(nw, x, "casual", "nodematch", "race", diff = TRUE)
cas_nm_ag <- calc_targets(nw, x, "casual", "nodematch", "age_group", diff = TRUE)
cas_nc_age <- calc_targets(nw, x, "casual", "nodecov", "age")
cas_nc_agesq <- calc_targets(nw, x, "casual", "nodecov", "age", attr_squared = TRUE)
cas_conc <- calc_targets(nw, x, "casual", "concurrent")

cas_targets <- unname(c(
  cas_edges,
  cas_nm_race[c(1:2, 4)],
  cas_nf_race[c(1:2, 4)],
  cas_nf_agegp[1:5],
  cas_conc
))

### Additional Arguments
cas_offset <- c(-Inf)

cas_diss <- dissolution_coefs(
  dissolution = ~ offset(edges),
  duration = x$casual$duration$overall,
  d.rate = 0 # few people aging out have casual rels, don't need to correct duration
)

cas_constraints <- ~
  sparse +
    blocks(attr = ~female, levels2 = diag(TRUE, 2)) +
    strat(attr = ~race, pmat = matrix( # should be same as x$casual$mixmat$race
      c(
        0.1354340, 0.0435742, 0.0055867, 0.0050679,
        0.0435742, 0.0970881, 0.0294654, 0.0386104,
        0.0055867, 0.0294654, 0.0187281, 0.0441261,
        0.0050679, 0.0386104, 0.0441261, 0.4158885
      ),
      byrow = TRUE,
      nrow = 4, ncol = 4
    ))

### Fit
cas_netest <- EpiModel::netest(
  nw = nw,
  formation = cas_form,
  target.stats = cas_targets,
  coef.form = cas_offset,
  coef.diss = cas_diss,
  constraints = cas_constraints,
  set.control.ergm = ergm_controls
)

### Set degree attributes -----------------------------------
nw <- set.vertex.attribute(nw, "deg_casual", get_degree(cas_netest$newnetwork))

## Main ---------------------------------------------------------
### Formation model
main_form <- ~ edges +
  nodematch(~race, diff = TRUE, levels = c(1:2, 4)) +
  nodefactor(~race, levels = c(1:2, 4)) +
  nodecov(~age) +
  nodecov(~agesq) +
  absdiff(~ sqrt(age)) +
  offset(nodefactor(~ deg_casual > 0)) +
  offset(nodefactor(~ olderpartner > 0))

### Calc and Compile Targets
main_edges <- calc_targets(nw, x, "main", "edges")
main_nm_race <- calc_targets(nw, x, "main", "nodematch", "race", diff = TRUE)
main_nf_race <- calc_targets(nw, x, "main", "nodefactor", "race")
main_nc_age <- calc_targets(nw, x, "main", "nodecov", "age")
main_nc_agesq <- calc_targets(nw, x, "main", "nodecov", "age", attr_squared = TRUE)
main_agemix <- calc_targets(nw, x, "main", "absdiff_sqrt_age")

main_targets <- unname(c(
  main_edges,
  main_nm_race[c(1:2, 4)],
  main_nf_race[c(1:2, 4)],
  main_nc_age,
  main_nc_agesq,
  main_agemix
))

### Additional Arguments
main_offset <- rep(-Inf, 2)

main_diss <- dissolution_coefs(
  dissolution = ~ offset(edges),
  duration = x$main$duration$overall,
  d.rate = dept_rate * 1.25 # to account for more relationships among older people in main netwrok
)

main_constraints <- ~
  bd(maxout = 1) +
    sparse +
    blocks(attr = ~female, levels2 = diag(TRUE, 2)) +
    strat(attr = ~race, pmat = matrix( # should be same as x$main$mixmat$race
      c(
        0.057178141, 0.02842713, 0.00417141, 0.003909043,
        0.028427131, 0.08849135, 0.02961956, 0.039865982,
        0.004171410, 0.02961956, 0.02019247, 0.049425911,
        0.003909043, 0.03986598, 0.04942591, 0.523299963
      ),
      byrow = TRUE,
      nrow = 4, ncol = 4
    ))

### Fit
main_netest <- EpiModel::netest(
  nw = nw,
  formation = main_form,
  target.stats = main_targets,
  coef.form = main_offset,
  coef.diss = main_diss,
  constraints = main_constraints,
  set.control.ergm = ergm_controls
)

### Set degree attributes -----------------------------------
nw <- set.vertex.attribute(nw, "deg_main", get_degree(main_netest$newnetwork))

## Inst network ----------------------------------------------------
### Formation model
inst_form <- ~ edges +
  nodefactor(~race, levels = c(1:2, 4)) +
  nodefactor(~age_group, levels = 1:3) +
  absdiff(~ sqrt(age))

### Calc and Compile Targets
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

### Additional Arguments
inst_offset <- c() # no offsets for inst network in formation model
inst_diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 1, d.rate = 0)
inst_constraints <- ~ sparse + blocks(attr = ~female, levels2 = diag(TRUE, 2))

### Fit
inst_netest <- EpiModel::netest(
  nw = nw,
  formation = inst_form,
  target.stats = inst_targets,
  coef.form = inst_offset,
  coef.diss = inst_diss,
  constraints = inst_constraints,
  set.control.ergm = ergm_controls
)

# Save out ----------------------------------------------------------
est <- list(main_netest, cas_netest, inst_netest)
if (!dir.exists(here::here("networks", "fits", Sys.Date()))) {
  dir.create(here::here("networks", "fits", Sys.Date()), recursive = TRUE)
}
saveRDS(est, file = here::here("networks", "fits", Sys.Date(), "nw.rds"))
