library(EpiModel)
devtools::load_all("epimodel-sti")

# Required objects
this_seed <- 12345
ncores <- parallel::detectCores() - 1L
x <- yaml::read_yaml(here::here("params", "nw_params.yaml"))
main_mixmat <- readRDS(here::here("params", "main_mixmat.rds"))
cas_mixmat <- readRDS(here::here("params", "casual_mixmat.rds"))
nw <- generate_init_network(x, seed = this_seed, assign_deg_casual = TRUE)
# Estimate Networks ----------------------------------------------
## Main ---------------------------------------------------------
main_form <- ~ edges +
  nodefactor(~race, levels = c(1:2, 4)) +
  nodefactor(~ floor(age), levels = 3:34) +
  nodematch(~race, diff = TRUE, levels = c(1:2, 4)) +
  absdiff(~ sqrt(age_adj)) +
  offset(nodefactor(~ deg_casual > 0)) +
  offset(nodefactor(~ floor(age), levels = 1:2))

## --- main targets ------------------------
main_edges <- calc_targets(nw, x, "main", "edges")
main_nf_race <- calc_targets(nw, x, "main", "nodefactor", "race")
main_nf_age <- calc_targets(nw, x, "main", "nodefactor", "age")
main_nm_race <- calc_targets(nw, x, "main", "nodematch", "race")
main_agemix <- calc_targets(nw, x, "main", "absdiff_sqrt_age")

main_targets <- unname(c(
  main_edges,
  main_nf_race[c(1:2, 4)],
  main_nf_age[-c(1:2, 35)],
  main_nm_race[c(1:2, 4)],
  main_agemix
))

main_offset <- c(rep(-Inf, 3))

main_diss <- dissolution_coefs(
  dissolution = ~ offset(edges),
  duration = x$main$duration$duration,
  d.rate = 0
)

# no same-sex ties, no concurrency in main network
main_constraints <- ~
  bd(maxout = 1) +
    sparse +
    blocks(attr = ~female, levels2 = diag(TRUE, 2)) +
    strat(attr = ~race, pmat = main_mixmat)

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

nw %v% "deg_main" <- get_degree(main_netest$newnetwork)

## Casual ----------------------------------------------------

cas_form <- ~ edges +
  nodefactor(~race, levels = c(1:2, 4)) +
  nodefactor(~ floor(age), levels = c(1:8)) +
  nodematch(~race, diff = TRUE, levels = c(1:2, 4)) +
  absdiff(~ sqrt(age_adj)) +
  concurrent() +
  offset(nodefactor(~deg_main))

## -- casual targets -------------------------------------
cas_edges <- calc_targets(nw, x, "casual", "edges")
cas_nf_race <- calc_targets(nw, x, "casual", "nodefactor", "race")
cas_nf_age <- calc_targets(nw, x, "casual", "nodefactor", "age")
cas_nm_race <- calc_targets(nw, x, "casual", "nodematch", "race")
cas_agemix <- calc_targets(nw, x, "casual", "absdiff_sqrt_age")
cas_conc <- calc_targets(nw, x, "casual", "concurrent")

cas_targets <- unname(c(
  cas_edges,
  cas_nf_race[c(1:2, 4)],
  cas_nf_age[c(1:8)],
  cas_nm_race[c(1:2, 4)],
  cas_agemix,
  cas_conc
))

cas_offset <- c(-Inf)

cas_diss <- dissolution_coefs(
  dissolution = ~ offset(edges),
  duration = x$casual$duration$duration,
  d.rate = 0
)

cas_constraints <- ~
  bd(maxout = 3) +
    sparse +
    blocks(attr = ~female, levels2 = diag(TRUE, 2)) +
    strat(attr = ~race, pmat = cas_mixmat)


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

nw %v% "deg_casual" <- get_degree(cas_netest$newnetwork)

# Inst network ----------------------------------------------------
inst_form <- ~ edges +
  nodefactor(~race, levels = c(1:2, 4)) +
  nodefactor(~age_group, levels = 1:3) +
  absdiff(~ sqrt(age_adj))

inst_edges <- calc_targets(nw, x, "inst", "edges", inst_correct = TRUE)
inst_nf_race <- calc_targets(nw, x, "inst", "nodefactor", "race", inst_correct = TRUE)
inst_nf_agegrp <- calc_targets(nw, x, "inst", "nodefactor", "age_group",
  inst_correct = TRUE
)
inst_agemix <- calc_targets(nw, x, "inst", "absdiff_sqrt_age", inst_correct = TRUE)

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


casdx <- cas_dynamic$nw
maindx <- main_dynamic$nw
datdx <- data.frame(
  age = floor(casdx %v% "age"),
  age_group = casdx %v% "age_group",
  race = casdx %v% "race",
  d = get_degree(maindx),
  c = get_degree(casdx)
)

compdx <- datdx |>
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

compdx |>
  dplyr::filter(type %in% c("main", "main_mod")) |>
  ggplot2::ggplot(ggplot2::aes(x = age, y = vals, col = type)) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(span = 0.75) +
  ggplot2::facet_wrap(~race)

compdx |>
  dplyr::filter(type %in% c("cas", "cas_mod")) |>
  ggplot2::ggplot(ggplot2::aes(x = age, y = vals, col = type)) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(span = 0.75) +
  ggplot2::facet_wrap(~race)

compdx |>
  ggplot2::ggplot(ggplot2::aes(x = age, y = vals, col = type)) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(span = 0.75) +
  ggplot2::facet_wrap(~race) +
  ggplot2::ggtitle("Comparing Empirical Relationship Activity to Fit Networks") +
  ggplot2::ylab("Proportion of persons in a given relationship") +
  ggplot2::xlab("Age")

ggplot2::ggsave(
  filename = here::here("networks", "fits", Sys.Date(), "post_diagnostics_comps.PNG"),
  width = 10, height = 8
)
