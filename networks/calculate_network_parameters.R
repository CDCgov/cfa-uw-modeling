# Here we set up list with information about desired population characteristics
# Will be used to generate starting network and ergm target statistics
# for network fitting process
# output: yaml config file

# -----------------------------------------------------------------
# Set up list and declare estimate type, load data
# -----------------------------------------------------------------
estimate_type <- "empirical" # , #"predicted" #alt:
out <- list()

w <- readRDS(here::here("data", "nsfg_wide.rds"))
l <- readRDS(here::here("data", "nsfg_long.rds"))
wsvy <- srvyr::as_survey_design(w, weights = weight)
lsvy <- srvyr::as_survey_design(l, weights = weight)

# -----------------------------------------------------------------
# Add enough information to parameterize starting nodal attributes
# -----------------------------------------------------------------
# requires some assumptions and previous knowledge of data
# required subset name: pop
# nolint start

out$pop$size <- 100000

out$pop$female$levels <- 0:1
out$pop$female$dist <- c(0.5, 0.5)

out$pop$age$min <- 15
out$pop$age$max <- 49
out$pop$age_group$levels <- 1:7
out$pop$age_group$dist <- c(rep(0.143, 6), 0.142) # even ages with slightly fewer 45-49

out$pop$race$levels <- c("B", "H", "O", "W")
out$pop$race$dist <- c(0.12, 0.19, 0.11, 0.58) # based on 2020 census

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
}

if (estimate_type == "predicted") {
  # -----------------------------------------------------------------
  # Estimate Partnership Frequencies
  # -----------------------------------------------------------------

  # First, we use the above empirical framework to get cross- and within-network concurrency rates
  # (using predicted values for these tend to need really detailed glms that don't smooth out activity in a desireable way)
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

  # Move on to smoothed predicted activity by age/race vals
  library(survey)
  wsvy <- wsvy |> dplyr::mutate(agesq = age^2)
  s <- wsvy$variables
  vars <- c("ego", "weight", "age", "agesq", "race", "female")
  testpop <- s |> dplyr::select(all_of(vars))
  testpop_18 <- testpop |> dplyr::filter(age < 19)
  testpop_19 <- testpop |> dplyr::filter(age >= 19)
  ## TRY WITH FITTING CASUAL FIRST

  cas_glm_18 <- svyglm(deg_casual ~ age + agesq + race + female,
    design = wsvy, subset = age < 19,
    family = quasipoisson()
  )

  cas_glm_19 <- svyglm(deg_casual ~ age + agesq + race + female, ,
    design = wsvy, subset = age >= 19,
    family = quasipoisson()
  )

  testpop_18$deg_casual_prob <- predict(cas_glm_18, newdata = testpop_18, type = "response")
  testpop_19$deg_casual_prob <- predict(cas_glm_19, newdata = testpop_19, type = "response")
  testpop2 <- rbind(testpop_18, testpop_19)

  main_glm <- svyglm(deg_main ~ age + agesq + race + female, design = wsvy, family = quasibinomial())
  testpop2$deg_main_prob <- predict(main_glm, newdata = testpop2, type = "response")

  inst_glm <- svyglm(deg_inst_high ~ age + agesq + race + female, design = wsvy, family = quasipoisson())
  testpop2$deg_inst_prob <- predict(inst_glm, newdata = testpop2, type = "response")

  fullpop <- testpop2 |>
    dplyr::mutate(
      deg_main = rbinom(nrow(testpop2), 1, deg_main_prob),
      deg_casual = rpois(nrow(testpop2), deg_casual_prob),
      deg_inst = rpois(nrow(testpop2), deg_inst_prob)
    )



  library(tidyverse)
  emp <- wsvy |>
    mutate(age_floor = floor(age)) |>
    group_by(age_floor, race) |>
    summarize(
      deg_main = srvyr::survey_mean(deg_main, vartype = NULL),
      deg_casual = srvyr::survey_mean(deg_casual, vartype = NULL),
      deg_inst = srvyr::survey_mean(deg_inst_high, vartype = NULL)
    ) |>
    mutate(type = "empirical")

  pred <- fullpop |>
    mutate(age_floor = floor(age)) |>
    group_by(age_floor, race) |>
    summarize(
      deg_main = mean(deg_main_prob),
      deg_casual = mean(deg_casual_prob),
      deg_inst = mean(deg_inst_prob)
    ) |>
    mutate(type = "predicted")

  all <- rbind(pred, emp)

  ### plot comparison
  all |>
    ggplot(aes(x = age_floor, y = deg_main, col = type)) +
    geom_point() +
    facet_wrap(~race)

  all |>
    ggplot(aes(x = age_floor, y = deg_casual, col = type)) +
    geom_point() +
    facet_wrap(~race)

  all |> ggplot(aes(x = age_floor, y = deg_inst, col = type)) +
    geom_point() +
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

  # Mean absolute value of difference in sqrt of age between partners (for absdiff ERGM term)
  # with female age adjusted to account for asymmmetry (they are usually younger than male partners)
  l <- lsvy$variables
  vars_long <- c("ego", "weight", "female", "age", "asym_agediff", "rel2", "race", "partdur")
  testpop_long <- l |>
    dplyr::select(all_of(vars_long)) |>
    dplyr::filter(!is.na(rel2)) # remove 1 rel where rel2 is NA


  adjs_glm <- svyglm(asym_agediff ~ female, design = lsvy)
  testpop_long$pred_asym <- predict(adjs_glm, type = "response", newdata = testpop_long)
  pred_adj <- testpop_long |>
    srvyr::as_survey_design(weights = weight) |>
    summarize(mean = srvyr::survey_mean(pred_asym, vartype = NULL)) |>
    pull()

  lsvy <- lsvy |>
    dplyr::mutate(
      age_adj1 = ifelse(female == 1, age + as.numeric(pred_adj), age),
      age_adj2 = ifelse(female == 0, alter_age + as.numeric(pred_adj), alter_age),
      diff_sqrt_age = abs(sqrt(age_adj1) - sqrt(age_adj2))
    )


  agemix_glm <- svyglm(diff_sqrt_age ~ female + rel2, design = lsvy)
  testpop_long$agemix <- predict(agemix_glm, type = "response", newdata = testpop_long)

  agemix_per_rel <- testpop_long |>
    srvyr::as_survey_design(weights = weight) |>
    group_by(rel2) |>
    summarize(mean = srvyr::survey_mean(agemix, vartype = NULL)) |>
    pull()

  # Update parameters
  out$pop$age$female_age_adj <- pred_adj # for pop setup
  out$main$absdiff_sqrt_age <- agemix_per_rel[1] # diff in sqrt age, main
  out$casual$absdiff_sqrt_age <- agemix_per_rel[2] # diff in sqrt age, casual
  # no information about age mixing for non-current (one-time) partners
  # use current casuals as proxy
  out$inst$absdiff_sqrt_age <- agemix_per_rel[2]

  # Race Matching (for nodematch ERGM term)
  race_glm <- svyglm(race_match ~ race + rel2, design = lsvy, family = quasibinomial())
  testpop_long$race_match <- predict(race_glm, newdata = testpop_long, type = "response")

  rm <- testpop_long |>
    srvyr::as_survey_design(weights = weight) |>
    group_by(rel2, race) |>
    summarize(mean = srvyr::survey_mean(race_match, vartype = NULL)) |>
    pull()

  # Update parameters
  out$main$nodematch$race <- rm[1:4]
  out$casual$nodematch$race <- rm[5:8]

  # --------------------------------------------------
  # Estimate Relationship Duration
  # --------------------------------------------------
  dur_glm <- svyglm(partdur ~ rel2 + age + race, design = lsvy)
  testpop_long$dur <- predict(dur_glm, newdata = testpop_long, type = "response")

  durs <- testpop_long |>
    srvyr::as_survey_design(weights = weight) |>
    group_by(rel2) |>
    summarize(mean = srvyr::survey_median(dur, vartype = NULL)) |>
    pull() * 4

  # Update parameters
  out$main$duration$duration <- durs[1]
  out$casual$duration$duration <- durs[2]
  out$main$duration$metric <- "weeks"
  out$casual$duration$metric <- "weeks"
}


# save out as yaml
params_name <- paste0("nw_params_", estimate_type, ".yaml")
yaml::write_yaml(out, here::here("params", params_name))

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

casual_race_mixmat <- matrix_symmetrical(c_mixmat_unadj)

saveRDS(main_race_mixmat, here::here("params", "main_mixmat.rds"))
saveRDS(casual_race_mixmat, here::here("params", "casual_mixmat.rds"))
# nolint end
