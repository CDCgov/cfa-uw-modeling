# Here we set up list with information about desired population characteristics
# Will be used to generate starting network and ergm target statistics
# for network fitting process
# output: yaml config file

# -----------------------------------------------------------------
# Set up list and declare estimate type, load data
# -----------------------------------------------------------------
estimate_type <- "predicted" # alt: empirical"
w <- readRDS(here::here("data", "nsfg_wide.rds"))
l <- readRDS(here::here("data", "nsfg_long.rds"))

out <- list()
out$pop$size <- 100000

out$pop$female$levels <- 0:1
out$pop$female$dist <- c(0.5, 0.5)

out$pop$age$min <- 15
out$pop$age$max <- 49
out$pop$age_group$levels <- 1:7
out$pop$age_group$dist <- c(rep(0.143, 6), 0.142) # even ages with slightly fewer 45-49

out$pop$race$levels <- c("B", "H", "O", "W")
out$pop$race$dist <- c(0.12, 0.19, 0.11, 0.58) # based on 2020 census

# rake surveys to population margins
pop_agegroup <- data.frame(age_group = names(table(w$age_group)), Freq = out$pop$age_group$dist * 144764299)
pop_race <- data.frame(race = out$pop$race$levels, Freq = out$pop$race$dist * 144764299)

wsvy <- srvyr::as_survey_design(w, weights = weight)
lsvy <- srvyr::as_survey_design(l, weights = weight)

wsvy <- survey::rake(wsvy, list(~age_group, ~race), list(pop_agegroup, pop_race))
lsvy <- survey::rake(lsvy, list(~age_group, ~race), list(pop_agegroup, pop_race))

## If estimate_type = "empirical", use weighted survey data to calculate estimates
## If estimate_type = "predicted", use glms based on weighted survey data to calculate smoothed estimates

if (estimate_type == "empirical") {
  # -----------------------------------------------------------------
  # Estimate Partnership Frequencies
  # -----------------------------------------------------------------

  ## 1. Proportion of respondents in a relationship, by type, age, and race
  ## 2. Proportion of respondents in a relationship, by main / casual (matrix)
  ## (used for within- and cross-network concurrency)
  ## 3. Instantaneous activity in last year

  ## 1. Proportion of respondents in a relationship, by type, age, and race
  ### fx accounts for casual network having >1 active rel, and >1 inst in last year

  target_props <- function(x, attr, degtype, attr2) {
    out <- x |>
      dplyr::mutate(age = floor(age)) |>
      dplyr::group_by({{ attr }}, {{ attr2 }}) |>
      dplyr::summarize(
        sum = srvyr::survey_total({{ degtype }},
          vartype = NULL
        ),
        pop = srvyr::survey_total(vartype = NULL)
      ) |>
      dplyr::mutate(prop = sum / pop) |> # nolint
      dplyr::select(-c("sum", "pop"))

    return(out)
  }

  out$main$nodefactor$age_race <- unname(unlist(
    target_props(wsvy, attr = age, attr2 = race, degtype = deg_main)$prop
  ))
  out$casual$nodefactor$age_race <- unname(unlist(
    target_props(wsvy, attr = age, attr2 = race, degtype = deg_casual)$prop
  ))
  out$inst$nodefactor$age_race <- unname(unlist(
    target_props(wsvy, attr = age, attr2 = race, degtype = deg_inst_high)$prop
  ))

  ## 2. Proportion of respondents in a relationship, by main / casual (matrix)
  degmat <- wsvy |>
    dplyr::group_by(deg_main, deg_casual) |>
    dplyr::summarize(sum = srvyr::survey_total(vartype = NULL)) |>
    dplyr::ungroup() |>
    dplyr::mutate(prop = sum / sum(sum)) |>
    dplyr::select(-sum) |>
    tidyr::pivot_wider(names_from = deg_casual, values_from = prop) |>
    dplyr::select(-deg_main) |>
    as.matrix()

  ## add overall concurrency types to out list
  out$main$has_casual <- degmat[2, 2] + degmat[2, 3]
  out$casual$has_main <- degmat[2, 2] + degmat[2, 3]
  out$casual$concurrent <- sum(colSums(degmat, na.rm = TRUE)[3:4])

  ## 3. Proportion of respondents with an inst rels
  inst <- wsvy |>
    dplyr::group_by(deg_inst_high) |>
    dplyr::summarize(prop = srvyr::survey_prop(vartype = NULL))

  inst_activity_year <-
    inst$prop[2] + 2 * inst$prop[3] + 3 * inst$prop[4] +
    4 * inst$prop[5] + 5 * inst$prop[6] + 7 * inst$prop[8]

  ### add to out list
  out$inst$density <- inst_activity_year
  out$inst$summary_time <- "year"

  # -----------------------------------------------------------------
  # Estimate Within-Partnership Characteristics
  # -----------------------------------------------------------------
  ## If estimate_type = "empirical", use weighted survey data to calculate estimates
  ## If estimate_type = "predicted", use glms based on weighted survey data to calculate smoothed estimates
  ## NOTE: for within partnership characteristcs, both methods produce similar results relative to above

  # Mean absolute value of difference in sqrt of age between partners (for absdiff ERGM term)
  # with female age adjusted to account for asymmmetry (they are usually younger than male partners)
  adjs <- lsvy |>
    dplyr::summarize(mean = srvyr::survey_mean(asym_agediff, na.rm = TRUE, vartype = NULL))

  lsvy <- lsvy |>
    dplyr::mutate(
      age_adj1 = ifelse(female == 1, age + as.numeric(adjs), age),
      age_adj2 = ifelse(female == 0, alter_age + as.numeric(adjs), alter_age)
    )

  agemix_vars <- lsvy |>
    dplyr::mutate(diff_sqrt_age = abs(sqrt(age_adj1) - sqrt(age_adj2))) |>
    dplyr::group_by(rel2) |>
    dplyr::summarize(abs_diff_sqrt_age = srvyr::survey_mean(diff_sqrt_age,
      na.rm = TRUE,
      vartype = NULL
    ))

  out$pop$age$female_age_adj <- as.numeric(adjs) # for pop setup
  out$main$absdiff_sqrt_age <- as.numeric(agemix_vars[1, 2]) # diff in sqrt age, main
  out$casual$absdiff_sqrt_age <- as.numeric(agemix_vars[2, 2]) # diff in sqrt age, casual
  # no information about age mixing for non-current (one-time) partners
  # use current casuals as proxy
  out$inst$absdiff_sqrt_age <- as.numeric(agemix_vars[2, 2])

  # Race Matching (for nodematch ERGM term)
  main_rm <- lsvy |>
    dplyr::filter(rel2 == "Marriage/Cohab") |>
    dplyr::group_by(race) |>
    dplyr::summarize(prop = srvyr::survey_mean(race_match,
      na.rm = TRUE,
      vartype = NULL
    ))

  casual_rm <- lsvy |>
    dplyr::filter(rel2 == "Casual/Other", curr == 0) |>
    dplyr::group_by(race) |>
    dplyr::summarize(prop = srvyr::survey_mean(race_match,
      na.rm = TRUE,
      vartype = NULL
    ))

  out$main$nodematch$race <- as.numeric(main_rm$prop)
  out$casual$nodematch$race <- as.numeric(casual_rm$prop)

  # --------------------------------------------------
  # Estimate Relationship Duration
  # --------------------------------------------------
  out$main$duration$duration <- lsvy |>
    dplyr::filter(rel2 == "Marriage/Cohab") |>
    dplyr::summarize(dur = srvyr::survey_median(partdur,
      na.rm = TRUE,
      vartype = NULL
    ) * 4) |>
    as.numeric()

  out$casual$duration$duration <- lsvy |>
    dplyr::filter(rel2 == "Casual/Other", curr == 1) |>
    dplyr::summarize(dur = srvyr::survey_median(partdur,
      na.rm = TRUE,
      vartype = NULL
    ) * 4) |>
    as.numeric()
  out$main$duration$metric <- "weeks"
  out$casual$duration$metric <- "weeks"

  ## race mixing to help ergm fit
  devtools::load_all("epimodel-sti")

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

  # not sure why there's still an NA column, so remove
  m_mixmat_unadj <- m_mixmat_unadj[, -5]

  main_race_mixmat <- matrix_symmetrical(m_mixmat_unadj)

  c_mixmat_unadj <- lsvy |>
    dplyr::filter(rel2 == "Casual/Other", !is.na(alter_race), curr == 1) |>
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

  c_mixmat_unadj <- c_mixmat_unadj[, -5]
  casual_race_mixmat <- matrix_symmetrical(c_mixmat_unadj)

  saveRDS(main_race_mixmat, here::here("params", "main_mixmat_empirical.rds"))
  saveRDS(casual_race_mixmat, here::here("params", "casual_mixmat_empirical.rds"))
}

if (estimate_type == "predicted") {
  # -----------------------------------------------------------------
  # Estimate Partnership Frequencies
  # -----------------------------------------------------------------
  # First: cross- and within- network concurrenct
  wsvy <- wsvy |> dplyr::mutate(
    agesq = age^2,
    cross_network = ifelse(deg_main > 0 & deg_casual > 0, 1, 0),
    conc = ifelse(deg_casual > 1, 1, 0)
  )
  has_cn_glm <- survey::svyglm(cross_network ~ age + race, design = wsvy, family = quasibinomial())
  wsvy$variables$has_cn <- predict(has_cn_glm, newdata = wsvy$variables, type = "response")

  cn_prop <- wsvy |>
    dplyr::summarize(has_cn = srvyr::survey_mean(has_cn, vartype = NULL)) |>
    dplyr::pull()

  has_conc_glm <- survey::svyglm(conc ~ age + race, design = wsvy, family = quasibinomial())
  wsvy$variables$has_conc <- predict(has_conc_glm, newdata = wsvy$variables, type = "response")

  conc_prop <- wsvy |>
    dplyr::summarize(has_conc = srvyr::survey_mean(has_conc, vartype = NULL)) |>
    dplyr::pull()

  ## Update Params
  out$main$cross_network <- cn_prop
  out$casual$cross_network <- cn_prop
  out$casual$concurrent <- conc_prop

  # Second: smoothed predicted activity by age/race vals
  cas_glm <- survey::svyglm(deg_casual ~ age + agesq + deg_main,
    design = wsvy, family = quasipoisson()
  )

  wsvy$variables$deg_casual_pred <- predict(cas_glm, newdata = wsvy$variables, type = "response")
  wsvy$variables$deg_casual2 <- rpois(nrow(wsvy$variables), wsvy$variables$deg_casual_pred)

  main_glm <- survey::svyglm(deg_main ~ age + race + agesq + female + deg_casual,
    design = wsvy, family = quasibinomial()
  )
  wsvy$variables$deg_main_prob <- predict(main_glm, newdata = wsvy$variables, type = "response")

  inst_glm <- survey::svyglm(deg_inst_high ~ age + agesq + race + female, design = wsvy, family = quasipoisson())
  wsvy$variables$deg_inst_prob <- predict(inst_glm, newdata = wsvy$variables, type = "response")

  wsvy$variables <- wsvy$variables |>
    dplyr::mutate(
      deg_main2 = rbinom(nrow(wsvy$variables), 1, deg_main_prob),
      deg_inst2 = rpois(nrow(wsvy$variables), deg_inst_prob)
    )

  emp <- wsvy |>
    dplyr::mutate(age_floor = floor(age)) |>
    dplyr::group_by(age_floor, race) |>
    dplyr::summarize(
      deg_main = srvyr::survey_mean(deg_main, vartype = NULL),
      deg_casual = srvyr::survey_mean(deg_casual, vartype = NULL),
      deg_inst = srvyr::survey_mean(deg_inst_high, vartype = NULL)
    ) |>
    dplyr::mutate(type = "empirical")

  pred <- wsvy |>
    dplyr::mutate(age_floor = floor(age)) |>
    dplyr::group_by(age_floor, race) |>
    dplyr::summarize(
      deg_main = srvyr::survey_mean(deg_main2, vartype = NULL),
      deg_casual = srvyr::survey_mean(deg_casual2, vartype = NULL),
      deg_inst = srvyr::survey_mean(deg_inst2, vartype = NULL)
    ) |>
    dplyr::mutate(type = "predicted")

  all <- rbind(pred, emp)

  ### plot comparison
  library(ggplot2)
  all |>
    ggplot(aes(x = age_floor, y = deg_main, col = type)) +
    geom_point() +
    geom_smooth() +
    facet_wrap(~race)

  all |>
    ggplot(aes(x = age_floor, y = deg_casual, col = type)) +
    geom_point() +
    geom_smooth() +
    facet_wrap(~race)

  all |> ggplot(aes(x = age_floor, y = deg_inst, col = type)) +
    geom_point() +
    geom_smooth() +
    facet_wrap(~race)

  # this looks a lot better than using the empirical estimates

  # add to out
  out$main$nodefactor$age_race <- pred$deg_main
  out$casual$nodefactor$age_race <- pred$deg_casual
  out$inst$nodefactor$age_race <- pred$deg_inst
  out$inst$summary_time <- "year"

  # -----------------------------------------------------------------
  # Estimate Within-Partnership Characteristics
  # -----------------------------------------------------------------
  ## If estimate_type = "empirical", use weighted survey data to calculate estimates
  ## If estimate_type = "predicted", use glms based on weighted survey data to calculate smoothed estimates
  ## NOTE: for within partnership characteristcs, both methods produce similar results relative to above

  ## Mean absolute value of difference in sqrt of age between partners (for absdiff ERGM term)
  ## with female age adjusted to account for asymmmetry (they are usually younger than male partners)

  ## Age Mixing ######################################################
  # predict asymmetric different in age (on average females younger)
  adjs_glm <- survey::svyglm(asym_agediff ~ female, design = lsvy)
  lsvy$variables$pred_asym <- predict(adjs_glm, type = "response", newdata = lsvy$variables)

  # estimate asymmetric age adjustment for females
  pred_adj <- lsvy |>
    dplyr::summarize(mean = srvyr::survey_mean(pred_asym, vartype = NULL)) |>
    dplyr::pull()

  # calculate adjusted ages and diff_sqrt_age based on age_adj
  lsvy <- lsvy |>
    dplyr::mutate(
      age_adj1 = ifelse(female == 1, age + pred_adj, age),
      age_adj2 = ifelse(female == 0, alter_age + pred_adj, alter_age),
      diff_sqrt_age = abs(sqrt(age_adj1) - sqrt(age_adj2))
    )

  # fit and predict
  agemix_glm <- survey::svyglm(diff_sqrt_age ~ female + rel2, design = lsvy)
  lsvy$variables$agemix <- predict(agemix_glm, type = "response", newdata = lsvy$variables)

  agemix_per_rel <- lsvy |>
    dplyr::group_by(rel2) |>
    dplyr::summarize(mean = srvyr::survey_mean(agemix, vartype = NULL)) |>
    dplyr::pull()

  ## Update parameters
  out$pop$age$female_age_adj <- pred_adj # for pop setup
  out$main$absdiff_sqrt_age <- agemix_per_rel[1] # diff in sqrt age, main
  out$casual$absdiff_sqrt_age <- agemix_per_rel[2] # diff in sqrt age, casual
  ## no information about age mixing for non-current (one-time) partners
  ## use current casuals as proxy
  out$inst$absdiff_sqrt_age <- agemix_per_rel[2]

  ## Age Group & Race Matching (for nodematch ERGM terms) #########################
  race_glm <- survey::svyglm(race_match ~ race + age_group + rel2, design = lsvy, family = quasibinomial())
  lsvy$variables$race_match_pred <- predict(race_glm, newdata = lsvy$variables, type = "response")

  rm <- lsvy |>
    dplyr::group_by(rel2, race) |>
    dplyr::summarize(mean = srvyr::survey_mean(race_match_pred, vartype = NULL)) |>
    dplyr::pull()

  agegrp_glm <- survey::svyglm(agegrp_match ~ age_group + race + rel2, design = lsvy, family = quasibinomial())
  lsvy$variables$agegrp_match_pred <- predict(agegrp_glm, newdata = lsvy$variables, type = "response")

  am <- lsvy |>
    dplyr::group_by(rel2, age_group) |>
    dplyr::summarize(mean = srvyr::survey_mean(agegrp_match_pred, vartype = NULL)) |>
    dplyr::pull()

  ## Update parameters
  out$main$nodematch$race <- rm[1:4]
  out$casual$nodematch$race <- rm[5:8]
  out$main$nodematch$age_group <- am[1:7]
  out$casual$nodematch$age_group <- am[8:14]

  ## Mixing matrices for constraints ##################################
  alter_race_glm <- survey::svyolr(as.factor(alter_race) ~ race + age + rel2,
    design = lsvy
  )
  lsvy$variables[, c("B", "H", "O", "W")] <- predict(alter_race_glm, newdata = lsvy$variables, type = "probs")

  totrels <- lsvy |>
    dplyr::group_by(rel2, race) |>
    dplyr::summarize(
      tot = srvyr::survey_total(vartype = NULL),
      B = srvyr::survey_mean(B, vartype = NULL),
      H = srvyr::survey_mean(H, vartype = NULL),
      O = srvyr::survey_mean(O, vartype = NULL),
      W = srvyr::survey_mean(W, vartype = NULL)
    ) |>
    dplyr::mutate(
      B = B * tot,
      H = H * tot,
      O = O * tot,
      W = W * tot
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-race, -rel2, -tot)

  mainrels <- totrels[1:4, ]
  casrels <- totrels[5:8, ]
  mainprops <- mainrels / sum(mainrels)
  casprops <- casrels / sum(casrels)

  devtools::load_all("epimodel-sti")
  main_race_mixmat <- matrix_symmetrical(mainprops)
  casual_race_mixmat <- matrix_symmetrical(casprops)

  saveRDS(main_race_mixmat, here::here("params", "main_mixmat_predicted.rds"))
  saveRDS(casual_race_mixmat, here::here("params", "casual_mixmat_predicted.rds"))

  alter_age_glm <- survey::svyolr(as.factor(alter_age_group) ~ race + age + rel2, design = lsvy)
  lsvy$variables[, paste0("ag", 1:7)] <- predict(alter_age_glm, newdata = lsvy$variables, type = "probs")

  atotrels <- lsvy |>
    dplyr::group_by(rel2, age_group) |>
    dplyr::summarize(
      tot = srvyr::survey_total(vartype = NULL),
      ag1 = srvyr::survey_mean(ag1, vartype = NULL),
      ag2 = srvyr::survey_mean(ag2, vartype = NULL),
      ag3 = srvyr::survey_mean(ag3, vartype = NULL),
      ag4 = srvyr::survey_mean(ag4, vartype = NULL),
      ag5 = srvyr::survey_mean(ag5, vartype = NULL),
      ag6 = srvyr::survey_mean(ag6, vartype = NULL),
      ag7 = srvyr::survey_mean(ag7, vartype = NULL)
    ) |>
    dplyr::mutate(
      ag1 = ag1 * tot,
      ag2 = ag2 * tot,
      ag3 = ag3 * tot,
      ag4 = ag4 * tot,
      ag5 = ag5 * tot,
      ag6 = ag6 * tot,
      ag7 = ag7 * tot
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-age_group, -rel2, -tot)

  amainrels <- atotrels[1:7, ]
  acasrels <- atotrels[8:14, ]
  amainprops <- amainrels / sum(amainrels)
  acasprops <- acasrels / sum(acasrels)

  main_ag_mixmat <- matrix_symmetrical(amainprops)
  casual_ag_mixmat <- matrix_symmetrical(acasprops)

  saveRDS(main_ag_mixmat, here::here("params", "main_mixmat_ag_predicted.rds"))
  saveRDS(casual_ag_mixmat, here::here("params", "casual_mixmat_ag_predicted.rds"))


  ## For nodecov target, mean of combined age and agesq across all rels
  comb_age_glm <- survey::svyglm(comb_age ~ age + race + rel2, design = lsvy)
  comb_agesq_glm <- survey::svyglm(comb_agesq ~ age + race + rel2, design = lsvy)
  lsvy$variables$comb_age_pred <- predict(comb_age_glm, newdata = lsvy$variables, type = "response")
  lsvy$variables$comb_agesq_pred <- predict(comb_agesq_glm, newdata = lsvy$variables, type = "response")

  comb_age_tars <- lsvy |>
    dplyr::group_by(rel2) |>
    dplyr::summarize(
      ca = srvyr::survey_mean(comb_age_pred, vartype = NULL),
      casq = srvyr::survey_mean(comb_agesq_pred, vartype = NULL)
    )

  ## Update
  out$main$nodecov$age <- as.numeric(comb_age_tars[1, 2])
  out$main$nodecov$agesq <- as.numeric(comb_age_tars[1, 3])
  out$casual$nodecov$age <- as.numeric(comb_age_tars[2, 2])
  out$casual$nodecov$agesq <- as.numeric(comb_age_tars[2, 3])

  # --------------------------------------------------
  # Estimate Relationship Duration
  # --------------------------------------------------
  lsvy <- lsvy |>
    dplyr::mutate(
      agegrp_dyad = paste0(age_group, alter_age_group),
      agegrp_dyad = ifelse(agegrp_match == TRUE, agegrp_dyad, NA),
      race_dyad = paste0(race, alter_race),
      race_dyad = ifelse(race_match == TRUE, race_dyad, NA)
    )

  dur_glm <- survey::svyglm(log(partdur) ~ rel2 + age + race, design = lsvy)
  lsvy$variables$dur_pred <- predict(dur_glm, newdata = lsvy$variables, type = "response")
  lsvy$variables$dur <- exp(lsvy$variables$dur_pred)

  durs_main <- lsvy |>
    dplyr::filter(rel2 == "Marriage/Cohab") |>
    dplyr::group_by(agegrp_dyad) |>
    dplyr::summarize(med = srvyr::survey_median(dur * (365 / 12), vartype = NULL)) |>
    dplyr::pull()

  durs_cas <- lsvy |>
    dplyr::filter(rel2 == "Casual/Other") |>
    dplyr::group_by(agegrp_dyad) |>
    dplyr::summarize(mean = srvyr::survey_median(dur * (365 / 12), vartype = NULL)) |>
    dplyr::pull()

  durs_main_single <- lsvy |>
    dplyr::filter(rel2 == "Marriage/Cohab") |>
    dplyr::summarize(med = srvyr::survey_median(dur * (365 / 12), vartype = NULL)) |>
    dplyr::pull()

  durs_cas_single <- lsvy |>
    dplyr::filter(rel2 == "Casual/Other") |>
    dplyr::summarize(med = srvyr::survey_median(dur * (365 / 12), vartype = NULL)) |>
    dplyr::pull()

  ## Update parameters
  out$main$duration$duration$agegrp_match <- durs_main
  out$casual$duration$duration$agegrp_match <- durs_cas
  out$main$duration$duration$overall <- durs_main_single
  out$casual$duration$duration$overall <- durs_cas_single
  out$main$duration$metric <- "days"
  out$casual$duration$metric <- "days"
}


# save out as yaml
params_name <- paste0("nw_params_", estimate_type, ".yaml")
yaml::write_yaml(out, here::here("params", params_name))

# nolint end
