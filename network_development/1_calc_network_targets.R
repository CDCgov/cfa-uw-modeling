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
param_folder <- here::here("input", "params")
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
# Fit GLMs for main, casual, inst degree, olderpartner (nodefactor) -------
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

main_glm <- survey::svyglm(
  deg_main ~ age_group + sqrt(age_group) + race,
  design = wsvy,
  family = quasipoisson()
)

casual_glm <- survey::svyglm(
  deg_casual ~ age_group + sqrt(age_group) + race,
  design = wsvy,
  family = quasipoisson()
)

inst_glm <- survey::svyglm(
  deg_inst_high ~ age_group + sqrt(age_group) + race,
  design = wsvy,
  family = quasipoisson()
)

prediction_pop$deg_main <- predict(
  main_glm,
  newdata = prediction_pop,
  type = "response"
)

prediction_pop$deg_casual <- predict(
  casual_glm,
  newdata = prediction_pop,
  type = "response"
)

prediction_pop$deg_inst <- predict(
  inst_glm,
  newdata = prediction_pop,
  type = "response"
)

out$main$nodefactor <- prediction_pop$deg_main |> as.numeric()
out$casual$nodefactor <- prediction_pop$deg_casual |> as.numeric()
out$inst$nodefactor <- prediction_pop$deg_inst |> as.numeric()

## Older partners
op_glm <- survey::svyglm(
  olderpartner ~ age_group + sqrt(age_group) + race,
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
# Casual Concurrency & Cross Network ties--------------------------
# -----------------------------------------------------------------
## Concurrency in casual network
conc_glm <- survey::svyglm(
  conc ~ age + age_group + race,
  design = wsvy,
  family = quasibinomial()
)

predicted <- predict(conc_glm, newdata = prediction_pop, type = "response")
out$casual$concurrent <- mean(predicted)

m_cross <- survey::svyglm(
  cross_network ~ deg_main,
  design = wsvy
)

m_predicted <- predict(
  m_cross,
  newdata = tibble(deg_main = 1),
  type = "response"
)

c_cross <- survey::svyglm(
  cross_network ~ deg_casual > 0,
  design = wsvy
)

c_predicted <- predict(
  c_cross,
  newdata = tibble(deg_casual = 1),
  type = "response"
)

out$main$cross_network <- as.numeric(m_predicted)
out$casual$cross_network <- as.numeric(c_predicted)

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

## by age group
## to match formation model we need a vector with
## first entry: median of all unmatched rels and matched rels of ag7
## median estimates for matched rels in ags 1-6

lsvy <- lsvy |>
  dplyr::mutate(
    ag_cat = ifelse(
      age_group < 7 & agegrp_match,
      as.character(age_group),
      "0"
    )
  )

durs <- lsvy |>
  dplyr::group_by(rel2, ag_cat) |>
  dplyr::summarize(dur = srvyr::survey_median(partdur, na.rm = TRUE))

m_durs_age_group <- durs |>
  filter(rel2 == "Marriage/Cohab") |>
  pull(dur)

out$main$duration$by_age_group <- m_durs_age_group * (365 / 12)

c_durs_age_group <- durs |>
  filter(rel2 == "Casual/Other") |>
  pull(dur)

out$casual$duration$by_age_group <- c_durs_age_group * (365 / 12)

# save out as yaml
params_name <- "nw_params.yaml"
yaml::write_yaml(out, file.path(param_folder, params_name))
