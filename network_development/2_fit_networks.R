library(epimodelcfa)

# -----------------------------------------------------------------
# Required inputs & setup -----------------------------------------
# -----------------------------------------------------------------

## Set out folder
out_folder <- here::here("input", "network_fits")
if (!dir.exists(out_folder)) {
  dir.create(out_folder, recursive = TRUE)
}

## Set seed
this_seed <- 11111

## Caclulate departure rate based on rate of aging out of pop
age_min <- 15 # Age at entry
age_max <- 50 # Age at exit
units_per_year <- 365 # Days

drate <- (1 / (age_max - age_min)) * (1 / units_per_year)

## Load behavioral params estimated from data
x <- yaml::read_yaml(here::here(
  "input",
  "params",
  "nw_params.yaml"
))


# -----------------------------------------------------------------
# Generate initial network ----------------------------------------
# -----------------------------------------------------------------

## Generates empty network with nodal attributes reflecting pop specs
nw <- generate_init_network(x, seed = this_seed)
## Set deg_main to 0 for preliminary casual network fitting
set.vertex.attribute(nw, "deg_main", rep(0, x$pop$size))

# -----------------------------------------------------------------
# Estimate Networks -----------------------------------------------
# -----------------------------------------------------------------

## Workflow:
## 1) Fit preliminary casual network to set plausible degree distribution
## 2) Add deg_casual attribute to nw object based on this fit
## 3) Fit main network (includes cross-network degree term in formation)
## 4) Update deg_main attribute based on this fit
## 5) Fit final casual network (incl cross-network degree term)
## 6) Fit inst network
## 7) Combine main, casual, and inst fits into list object, save out

## Prelim Casual Network -------------------------------------------
### Formation component
### (no cross-network term in prelim model)
cas_prelim_form <- ~ edges +
  nodematch(~age_group, diff = TRUE, levels = 1:6) +
  nodefactor(~ floor(age), levels = 1) +
  nodefactor(~age_group, levels = 1:6) +
  nodefactor(~race, levels = 1) +
  nodematch(~race) +
  concurrent()

### Calc and Compile Targets
cas_edges <- calc_targets(nw, x, "casual", "edges")
cas_nf_age <- calc_targets(nw, x, "casual", "nodefactor", "age")
cas_nf_age_group <- calc_targets(nw, x, "casual", "nodefactor", "age_group")
cas_nf_race <- calc_targets(nw, x, "casual", "nodefactor", "race")
cas_nm_race <- calc_targets(nw, x, "casual", "nodematch", "race")
cas_nm_age_group <- calc_targets(
  nw,
  x,
  "casual",
  "nodematch",
  "age_group",
  diff = TRUE
)
cas_conc <- calc_targets(nw, x, "casual", "concurrent")

cas_prelim_targets <- unname(c(
  cas_edges,
  cas_nm_age_group[1:6],
  cas_nf_age[1],
  cas_nf_age_group[1:6],
  cas_nf_race[1],
  cas_nm_race,
  cas_conc
))

## Dissolution component
### Re-arrange targets so that non-age-group-matched
### duration estimate is first in vector instead of last
durs <- c(
  x$casual$duration$by_age_group[8],
  x$casual$duration$by_age_group[-c(7:8)]
)

cas_diss <- dissolution_coefs(
  dissolution = ~ offset(edges) +
    offset(nodematch(~age_group, diff = TRUE, levels = 1:6)),
  duration = durs,
  d.rate = drate
)

## Network constraints
### no same-sex ties, max degree = 3, sparse network hint
cas_constraints <- ~ bd(maxout = 3) +
  sparse +
  blocks(attr = ~female, levels2 = diag(TRUE, 2))

### Fit Model
cas_prelim_netest <- EpiModel::netest(
  nw = nw,
  formation = cas_prelim_form,
  target.stats = cas_prelim_targets,
  coef.diss = cas_diss,
  constraints = cas_constraints,
  set.control.ergm = control.ergm(
    seed = this_seed,
    MCMC.prop = ~sparse,
    MCMC.interval = 2048 * 5,
    main.method = "Stochastic-Approximation"
  )
)

### Set deg_casual degree attributes to base nw --------------------
nw <- set.vertex.attribute(
  nw,
  "deg_casual",
  get_degree(cas_prelim_netest$newnetwork)
)

## Main Model -----------------------------------------------------
### Formation component
main_form <- ~ edges +
  nodematch(~age_group, diff = TRUE, levels = 1:6) +
  nodefactor(~age_group, levels = 1:6) +
  nodefactor(~race, levels = 1) +
  nodematch(~race) +
  absdiff(~ sqrt(age)) +
  nodefactor(~ deg_casual > 0)

### Calc and Compile Targets
main_edges <- calc_targets(nw, x, "main", "edges")
main_nf_age <- calc_targets(nw, x, "main", "nodefactor", "age")
main_nf_age_group <- calc_targets(nw, x, "main", "nodefactor", "age_group")
main_nf_race <- calc_targets(nw, x, "main", "nodefactor", "race")
main_nm_race <- calc_targets(nw, x, "main", "nodematch", "race")
main_nm_age_group <- calc_targets(
  nw,
  x,
  "main",
  "nodematch",
  "age_group",
  diff = TRUE
)
main_agemix <- calc_targets(nw, x, "main", "absdiff_sqrt_age")

main_targets <- unname(c(
  main_edges,
  main_nm_age_group[1:6],
  main_nf_age_group[1:6],
  main_nf_race[1],
  main_nm_race,
  main_agemix,
  x$main$cross * main_edges * 2
))

### Dissolution component
mdurs <- c(
  x$main$duration$by_age_group[8],
  x$main$duration$by_age_group[-c(7:8)]
)
main_diss <- dissolution_coefs(
  dissolution = ~ offset(edges) +
    offset(nodematch(~age_group, diff = TRUE, levels = 1:6)),
  duration = mdurs,
  d.rate = drate
)

### Main constraints
### no same-sex ties, max degree = 1, sparse network hint
main_constraints <- ~ bd(maxout = 1) +
  sparse +
  blocks(attr = ~female, levels2 = diag(TRUE, 2))

### Fit Model
main_netest <- EpiModel::netest(
  nw = nw,
  formation = main_form,
  target.stats = main_targets,
  coef.diss = main_diss,
  constraints = main_constraints,
  set.control.ergm = control.ergm(
    seed = this_seed,
    MCMC.prop = ~sparse,
    MCMC.interval = 2048 * 10,
    main.method = "Stochastic-Approximation"
  )
)

### Set deg_main degree attributes -----------------------------------
nw <- set.vertex.attribute(nw, "deg_main", get_degree(main_netest$newnetwork))

## Fit final casual network
### now includes cross-network degree term

### Updated formation component
cas_form <- ~ edges +
  nodematch(~age_group, diff = TRUE, levels = 1:6) +
  nodefactor(~ floor(age), levels = 1) +
  nodefactor(~age_group, levels = 1:6) +
  nodefactor(~race, levels = 1) +
  nodematch(~race) +
  concurrent() +
  nodefactor(~deg_main)

cas_targets <- unname(c(
  cas_edges,
  cas_nm_age_group[1:6],
  cas_nf_age[1],
  cas_nf_age_group[1:6],
  cas_nf_race[1],
  cas_nm_race,
  cas_conc,
  x$casual$cross * cas_edges * 2
))

## Fit Model
cas_netest <- EpiModel::netest(
  nw = nw,
  formation = cas_form,
  target.stats = cas_targets,
  coef.diss = cas_diss,
  constraints = cas_constraints,
  set.control.ergm = control.ergm(
    seed = this_seed,
    MCMC.prop = ~sparse,
    MCMC.interval = 2048 * 5,
    main.method = "Stochastic-Approximation"
  )
)

## Inst network ----------------------------------------------------
### Formation component
inst_form <- ~ edges +
  nodefactor(~race, levels = 1) +
  nodefactor(~age_group, levels = 1:3) +
  absdiff(~ sqrt(age))

### Calc and Compile Targets
inst_edges <- calc_targets(
  nw,
  x,
  "inst",
  "edges",
  inst_correct = TRUE,
  time_unit = "days"
)
inst_nf_race <- calc_targets(
  nw,
  x,
  "inst",
  "nodefactor",
  "race",
  inst_correct = TRUE,
  time_unit = "days"
)
inst_nf_agegrp <- calc_targets(
  nw,
  x,
  "inst",
  "nodefactor",
  "age_group",
  inst_correct = TRUE,
  time_unit = "days"
)
inst_agemix <- calc_targets(
  nw,
  x,
  "inst",
  "absdiff_sqrt_age",
  inst_correct = TRUE,
  time_unit = "days"
)

inst_targets <- unname(c(
  inst_edges,
  inst_nf_race[1],
  inst_nf_agegrp[1:3],
  inst_agemix
))

### Dissolution component
inst_diss <- dissolution_coefs(
  dissolution = ~ offset(edges),
  duration = 1,
  d.rate = 0
)

### Constraints
inst_constraints <- ~ sparse + blocks(attr = ~female, levels2 = diag(TRUE, 2))

### Fit
inst_netest <- EpiModel::netest(
  nw = nw,
  formation = inst_form,
  target.stats = inst_targets,
  coef.diss = inst_diss,
  constraints = inst_constraints,
  set.control.ergm = control.ergm(
    seed = this_seed,
    MCMC.prop = ~sparse,
    MCMC.interval = 2048,
    main.method = "Stochastic-Approximation"
  )
)

# ------------------------------------------------------------------
# Save out ----------------------------------------------------------
# ------------------------------------------------------------------

est <- list(main_netest, cas_netest, inst_netest)
saveRDS(
  est,
  file = file.path(out_folder, "nw.rds")
)
