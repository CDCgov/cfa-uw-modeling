library(tidyverse)
library(epimodelcfa)
# load data
# location of survey data
input_folder <- file.path(here::here("input", "nsfg"))

w <- readRDS(file.path(input_folder, "nsfg_wide.rds"))
## change levels to force W to be ref category in glms
w <- w |>
  mutate(race = factor(race, levels = c("W", "B", "O", "H")))
l <- readRDS(file.path(input_folder, "nsfg_long.rds"))
l <- l |>
  mutate(race = factor(race, levels = c("W", "B", "O", "H")))

# destination folder for parameter yaml
param_folder <- here::here("local_tests", "input", "params")
if (!dir.exists(param_folder)) {
  dir.create(param_folder, recursive = TRUE)
}

# Set population parameters ---------------------------
out <- list()
out$pop$size <- 50000

out$pop$female$levels <- 0:1
out$pop$female$dist <- c(0.5, 0.5)

out$pop$age$min <- 15
out$pop$age$max <- 49
out$pop$age_group$levels <- 1:7
# 15-18, 19-24, then 5-yr ints
out$pop$age_group$dist <- c(0.114, 0.171, rep(0.143, 5))

out$pop$race$levels <- c("B", "H", "O", "W")
out$pop$race$dist <- c(0.12, 0.19, 0.11, 0.58) # based on 2020 census

# Rake surveys to population margins ---------------------------
us_pop <- 144764299 # from census
pop_agegroup <- data.frame(
  age_group = out$pop$age_group$levels,
  Freq = out$pop$age_group$dist * us_pop
)
pop_race <- data.frame(
  race = out$pop$race$levels,
  Freq = out$pop$race$dist * us_pop
)

pop_sex <- data.frame(
  female = out$pop$female$levels,
  Freq = out$pop$female$dist * us_pop
)

wsvy <- srvyr::as_survey_design(
  w,
  ids = secu,
  strata = sest,
  weights = weight,
  nest = TRUE
)
lsvy <- srvyr::as_survey_design(
  l,
  ids = secu,
  strata = sest,
  weights = weight,
  nest = TRUE
)

wsvy <- survey::rake(
  wsvy,
  list(~age_group, ~race, ~female),
  list(pop_agegroup, pop_race, pop_sex)
)

lsvy <- survey::rake(
  lsvy,
  list(~age_group, ~race, ~female),
  list(pop_agegroup, pop_race, pop_sex)
)


# Add age squared and cross-network ties ---------------------------
wsvy <- wsvy |>
  dplyr::mutate(
    agesq = age^2,
    agecb = age^3,
    cross_network = ifelse(deg_main > 0 & deg_casual > 0, 1, 0),
    conc = ifelse(deg_casual > 1, 1, 0)
  )

lsvy <- lsvy |>
  dplyr::mutate(
    absdiff_sqrt_age = abs(sqrt(age) - sqrt(alter_age))
  )

# Make Prediction Pop DF for degree- related GLMs ------------------
ages <- seq(out$pop$age$min, (out$pop$age$max))
races <- out$pop$race$levels

prediction_pop <- tibble(
  age = rep(ages, each = length(races)),
  agesq = age^2,
  race = rep(races, times = length(ages)),
  age_group = case_when(
    age < 50 & age >= 45 ~ 7,
    age < 45 & age >= 40 ~ 6,
    age < 40 & age >= 35 ~ 5,
    age < 35 & age >= 30 ~ 4,
    age < 30 & age >= 25 ~ 3,
    age < 25 & age >= 19 ~ 2,
    age < 19 ~ 1
  )
)

out$prediction_pop$age <- prediction_pop$age
out$prediction_pop$agesq <- prediction_pop$agesq
out$prediction_pop$age_group <- prediction_pop$age_group
out$prediction_pop$race <- prediction_pop$race

# -------------------------------------------------------------------------
# Fit GLMs for main, casual, inst degree EDGES
# -------------------------------------------------------------------------
m_edges_mod <- survey::svyglm(
  deg_main ~ 1,
  design = wsvy,
  family = quasipoisson()
)

c_edges_mod <- survey::svyglm(
  deg_casual ~ 1,
  design = wsvy,
  family = quasipoisson()
)

i_edges_mod <- survey::svyglm(
  deg_inst_high ~ 1,
  design = wsvy,
  family = quasipoisson()
)

out$main$edges <- exp(coef(m_edges_mod)[[1]])
out$casual$edges <- exp(coef(c_edges_mod)[[1]])
out$inst$edges <- exp(coef(i_edges_mod)[[1]])


## Older partners
op_glm <- survey::svyglm(
  olderpartner ~ age + race,
  design = wsvy,
  family = quasibinomial()
)

prediction_pop$olderpartner <- predict(
  op_glm,
  newdata = prediction_pop,
  type = "response"
) #nolint
out$main$olderpartner <- prediction_pop$olderpartner |> as.numeric()

# -----------------------------------------------------------------------------
# Nodecov (mean of combined age and agesq across all rels)
# -----------------------------------------------------------------------------
## MAIN
mod <- survey::svyglm(
  comb_age ~ 1,
  design = lsvy,
  subset = rel2 == "Marriage/Cohab",
  family = quasipoisson()
)
pred <- exp(coef(mod)[[1]])
out$main$nodecov$age <- pred

mod <- survey::svyglm(
  comb_agesq ~ 1,
  design = lsvy,
  subset = rel2 == "Marriage/Cohab",
  family = quasipoisson()
)
pred <- exp(coef(mod)[[1]])
out$main$nodecov$agesq <- pred

## CASUAL
mod <- survey::svyglm(
  comb_age ~ 1,
  design = lsvy,
  subset = rel2 == "Casual/Other",
  family = quasipoisson()
)
pred <- exp(coef(mod)[[1]])
out$casual$nodecov$age <- pred

mod <- survey::svyglm(
  comb_agesq ~ 1,
  design = lsvy,
  subset = rel2 == "Casual/Other",
  family = quasipoisson()
)

pred <- exp(coef(mod)[[1]])
out$casual$nodecov$agesq <- pred

# -----------------------------------------------------------------------------
# Nodematch (Age Group and Race, global and per-group)
# -----------------------------------------------------------------------------
## Age group - by group
agegrp_glm_each <- survey::svyglm(
  agegrp_match ~ rel2 + as.factor(age_group),
  design = lsvy
)

agegrp_data <- tibble(
  rel2 = rep(c("Marriage/Cohab", "Casual/Other"), each = 7),
  age_group = rep(as.factor(1:7), 2)
)

agegrp_data$prob <- predict(
  agegrp_glm_each,
  newdata = agegrp_data,
  type = "response"
)

out$main$nodematch$age_group$each <- agegrp_data |>
  filter(rel2 == "Marriage/Cohab") |>
  pull(prob) |>
  as.numeric()

out$casual$nodematch$age_group$each <- agegrp_data |>
  filter(rel2 == "Casual/Other") |>
  pull(prob) |>
  as.numeric()

## Age group - global
agegrp_glm <- survey::svyglm(
  agegrp_match ~ rel2,
  design = lsvy
)
rel_cats <- tibble(rel2 = c("Marriage/Cohab", "Casual/Other"))
rel_cats$ag_match <- predict(agegrp_glm, newdata = rel_cats, type = "response")

out$main$nodematch$age_group$global <- rel_cats |>
  filter(rel2 == "Marriage/Cohab") |>
  pull(ag_match) |>
  as.numeric()

out$casual$nodematch$age_group$global <- rel_cats |>
  filter(rel2 == "Casual/Other") |>
  pull(ag_match) |>
  as.numeric()

## Race - by group
race_glm_each <- survey::svyglm(
  race_match ~ rel2 + as.factor(race),
  design = lsvy
)

race_data <- tibble(
  rel2 = rep(c("Marriage/Cohab", "Casual/Other"), each = 4),
  race = rep(c("B", "H", "O", "W"), 2)
)

race_data$prob <- predict(
  race_glm_each,
  newdata = race_data,
  type = "response"
)

out$main$nodematch$race$each <- race_data |>
  filter(rel2 == "Marriage/Cohab") |>
  pull(prob) |>
  as.numeric()

out$casual$nodematch$race$each <- race_data |>
  filter(rel2 == "Casual/Other") |>
  pull(prob) |>
  as.numeric()

## Race - global
race_glm <- survey::svyglm(
  race_match ~ rel2,
  design = lsvy
)
rel_cats <- tibble(rel2 = c("Marriage/Cohab", "Casual/Other"))
rel_cats$race_match <- predict(race_glm, newdata = rel_cats, type = "response")

out$main$nodematch$race$global <- rel_cats |>
  filter(rel2 == "Marriage/Cohab") |>
  pull(race_match) |>
  as.numeric()

out$casual$nodematch$race$global <- rel_cats |>
  filter(rel2 == "Casual/Other") |>
  pull(race_match) |>
  as.numeric()

# -----------------------------------------------------------------
# Casual Concurrency ----------------------------------------------
# -----------------------------------------------------------------
## Concurrency in casual network
conc_glm <- survey::svyglm(
  conc ~ age + age_group + race,
  design = wsvy,
  family = quasibinomial()
)

predicted <- predict(conc_glm, newdata = prediction_pop, type = "response")
out$casual$concurrent <- mean(predicted)

# -----------------------------------------------------------------
# Absdiff Age ----------------------------------------------
# -----------------------------------------------------------------

# Absdiff (mean abs val in sqrt of age between partners)
mod <- survey::svyglm(
  absdiff_sqrt_age ~ 1,
  design = lsvy,
  subset = rel2 == "Marriage/Cohab",
  family = quasipoisson()
)
pred <- exp(coef(mod)[[1]])
out$main$absdiff_sqrt_age <- pred


## Absdiff (mean abs val in sqrt of age between partners)
mod <- survey::svyglm(
  absdiff_sqrt_age ~ 1,
  design = lsvy,
  subset = rel2 == "Casual/Other" & curr == 1,
  family = quasipoisson()
)
pred <- exp(coef(mod)[[1]])
out$casual$absdiff_sqrt_age <- pred

## Absdiff (mean abs val in sqrt of age between partners)
## no info in data on age of only inst partners, so use casual network
out$inst$absdiff_sqrt_age <- out$casual$absdiff_sqrt_age

# Relationship Durations -----------------------------------
# calc median rel dur from the survey data
dur <- lsvy |>
  dplyr::filter(rel2 == "Marriage/Cohab") |>
  dplyr::summarize(dur = srvyr::survey_median(partdur, na.rm = TRUE)) |>
  dplyr::pull(dur)

out$main$duration$overall <- dur * (365 / 12)
out$main$duration$metric <- "days"

dur <- lsvy |>
  dplyr::filter(rel2 == "Casual/Other" & curr == 1) |>
  dplyr::summarize(dur = srvyr::survey_median(partdur, na.rm = TRUE)) |>
  dplyr::pull(dur)

out$casual$duration$overall <- dur * (365 / 12)
out$casual$duration$metric <- "days"

# Mixing matrices for constraints ##################################
alter_race_glm <- survey::svyolr(
  as.factor(alter_race) ~ race + age + rel2,
  design = lsvy
)
lsvy$variables[, c("B", "H", "O", "W")] <- predict(
  alter_race_glm,
  newdata = lsvy$variables,
  type = "probs"
)

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

main_race_mixmat <- epimodelcfa::matrix_symmetrical(mainprops)
casual_race_mixmat <- epimodelcfa::matrix_symmetrical(casprops)

out$main$mixmat$race <- main_race_mixmat
out$casual$mixmat$race <- casual_race_mixmat

alter_age_glm <- survey::svyolr(
  as.factor(alter_age_group) ~ race + age + rel2,
  design = lsvy
)
lsvy$variables[, paste0("ag", 1:7)] <- predict(
  alter_age_glm,
  newdata = lsvy$variables,
  type = "probs"
)

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

out$main$mixmat$age_group <- main_ag_mixmat
out$casual$mixmat$age_group <- casual_ag_mixmat


# save out as yaml
params_name <- "nw_params.yaml"
yaml::write_yaml(out, file.path(param_folder, params_name))
