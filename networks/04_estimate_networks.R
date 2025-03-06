library(EpiModel)
source(here::here("networks", "helper_functions.R"))

# Required objects in environemnt: params, starting network, seed
this_seed <- 12345
x <- yaml::read_yaml(here::here("params", "nw_params.yaml"))
nw <- generate_init_network(x, seed = this_seed, deg_casual = TRUE)

# Estimate Networks ----------------------------------------------
## Main ---------------------------------------------------------
main_form <- ~ edges +
  nodematch(~race, diff = TRUE, levels = c(1:2, 4)) +
  nodefactor(~race, levels = c(1:2, 4)) +
  nodefactor(~ floor(age), levels = 3:34) +
  absdiff(~ sqrt(age_adj)) +
  offset(nodefactor(~ deg_casual > 0)) +
  offset(nodefactor(~ floor(age), levels = 1:2))

## --- main targets ------------------------
main_edges <- calc_targets(rel = "main", count_type = "edges")
main_nodematch_race_diff <- calc_targets(rel = "main", count_type = "nodematch", attr_name = "race")
main_nodefactor_race <- calc_targets(rel = "main", count_type = "nodefactor", attr_name = "race")
main_nodefactor_age <- calc_targets(rel = "main", count_type = "nodefactor", attr_name = "age")
main_absdiff_sqrt_age_adj <- calc_targets(rel = "main", count_type = "absdiff_sqrt_age")


l <- readRDS(here::here("data", "nsfg_long.rds"))
lsvy <- srvyr::as_survey_design(l, weights = weight)
m_mixmat_unadj <- lsvy |>
  dplyr::filter(rel2 == "Marriage/Cohab", !is.na(alter_race)) |>
  dplyr::mutate(race_combo = paste0(race, alter_race)) |>
  dplyr::group_by(race, alter_race) |>
  dplyr::summarize(num = srvyr::survey_total(
    na.rm = TRUE,
    vartype = NULL
  )) |>
  dplyr::ungroup() |>
  dplyr::mutate(prop = num / sum(num)) |>
  dplyr::select(-num) |>
  tidyr::pivot_wider(names_from = alter_race, values_from = prop) |>
  dplyr::select(-race) |>
  as.matrix()

main_race_mixmat <- matrix_symmetrical(m_mixmat_unadj, nrow(m_mixmat_unadj))

main_targets <- unname(c(
  main_edges,
  main_nodematch_race_diff[c(1:2, 4)],
  main_nodefactor_race[c(1:2, 4)],
  main_nodefactor_age[-c(1:2, 35)],
  main_absdiff_sqrt_age_adj
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
    strat(attr = ~race, pmat = main_race_mixmat)

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
  dynamic = TRUE, nsims = 5, nsteps = 500, ncores = 5
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
  nodematch(~race, diff = TRUE, levels = c(1:2, 4)) +
  nodefactor(~race, levels = c(1:2, 4)) +
  nodefactor(~ floor(age), levels = c(1:8)) +
  absdiff(~ sqrt(age_adj)) +
  concurrent() +
  offset(nodefactor(~deg_main))

## -- casual targets -------------------------------------
cas_edges <- calc_targets(rel = "casual", count_type = "edges")
cas_nodefactor_race <- calc_targets(rel = "casual", count_type = "nodefactor", attr_name = "race")
cas_nodematch_race_diff <- calc_targets(rel = "casual", count_type = "nodematch", attr_name = "race")
cas_nodefactor_age <- calc_targets(rel = "casual", count_type = "nodefactor", attr_name = "age")
cas_absdiff_sqrt_age_adj <- calc_targets(rel = "casual", count_type = "absdiff_sqrt_age")
cas_conc <- calc_targets(rel = "casual", count_type = "concurrent")


c_mixmat_unadj <- lsvy |>
  dplyr::filter(rel2 == "Casual/Other", !is.na(alter_race) & curr == 1) |>
  dplyr::mutate(race_combo = paste0(race, alter_race)) |>
  dplyr::group_by(race, alter_race) |>
  dplyr::summarize(num = srvyr::survey_total(
    na.rm = TRUE,
    vartype = NULL
  )) |>
  dplyr::ungroup() |>
  dplyr::mutate(prop = num / sum(num)) |>
  dplyr::select(-num) |>
  tidyr::pivot_wider(names_from = alter_race, values_from = prop) |>
  dplyr::select(-race) |>
  as.matrix()

c_race_mixmat <- matrix_symmetrical(c_mixmat_unadj, nrow(c_mixmat_unadj))


cas_targets <- unname(c(
  cas_edges,
  cas_nodematch_race_diff[c(1:2, 4)],
  cas_nodefactor_race[c(1:2, 4)],
  cas_nodefactor_age[c(1:8)],
  cas_absdiff_sqrt_age_adj,
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
    strat(attr = ~race, pmat = c_race_mixmat)


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
  dynamic = TRUE, nsims = 5, nsteps = 500, ncores = 5
)
cas_dynamic

saveRDS(cas_netest, here::here("networks", "fits", Sys.Date(), "casual.rds"))

nw %v% "deg_casual" <- get_degree(cas_netest$newnetwork)

# Inst network ----------------------------------------------------
inst_form <- ~ edges +
  nodefactor(~race, levels = c(1:2, 4)) +
  nodefactor(~age_group, levels = 1:3) +
  absdiff(~ sqrt(age_adj))

inst_edges <- calc_targets(rel = "inst", count_type = "edges", inst_correct = TRUE)
inst_nodefactor_race <- calc_targets(rel = "inst", count_type = "nodefactor", attr_name = "race", inst_correct = TRUE)
inst_nodefactor_agegroup <- calc_targets(
  rel = "inst", count_type = "nodefactor",
  attr_name = "age_group", inst_correct = TRUE
)
inst_absdiff_sqrt_age_adj <- calc_targets(rel = "inst", count_type = "absdiff_sqrt_age", inst_correct = TRUE)

inst_targets <- unname(c(
  inst_edges,
  inst_nodefactor_race[c(1:2, 4)],
  inst_nodefactor_agegroup[1:3],
  inst_absdiff_sqrt_age_adj
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
