library(tidyverse)
# load data
# location of survey data
input_folder <- file.path(here::here("input", "nsfg"))

w <- readRDS(file.path(input_folder, "nsfg_wide.rds"))
l <- readRDS(file.path(input_folder, "nsfg_long.rds"))

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
# even ages with slightly fewer 45-49
out$pop$age_group$dist <- c(rep(0.143, 6), 0.142)

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

ages <- seq(out$pop$age$min, out$pop$age$max)
races <- out$pop$race$levels
females <- out$pop$female$levels

prediction_pop <- tibble(
  age = rep(rep(ages, each = length(races)), times = length(females)),
  agesq = age^2,
  race = rep(rep(races, times = length(ages)), times = length(females)),
  female = rep(females, each = length(ages) * length(races))
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

# Fit GLMs for main, casual, inst degree ---------------------------
main_glm <- survey::svyglm(
  deg_main ~ age + agesq + race,
  design = wsvy,
  family = quasipoisson()
)

casual_glm <- survey::svyglm(
  deg_casual ~ I(age) + race + female,
  design = wsvy,
  family = quasipoisson()
)

inst_glm <- survey::svyglm(
  deg_inst_high ~ age + agesq + race + female,
  design = wsvy,
  family = quasipoisson()
)

casual_conc_glm <- survey::svyglm(
  conc ~ age + agesq + race + female,
  design = wsvy,
  family = quasibinomial()
)

observed <- wsvy |>
  group_by(age, agesq, age_group, race, female) |>
  reframe(
    deg_main = srvyr::survey_mean(deg_main, na.rm = TRUE),
    deg_casual = srvyr::survey_mean(deg_casual, na.rm = TRUE),
    deg_inst = srvyr::survey_mean(deg_inst_high, na.rm = TRUE),
    conc = srvyr::survey_mean(conc, na.rm = TRUE),
  )

prediction_pop <- observed |>
  select(age, agesq, age_group, race, female)

prediction_pop$deg_main <- predict(main_glm, newdata = prediction_pop, type = "response") #nolint
prediction_pop$deg_casual <- predict(casual_glm, newdata = prediction_pop, type = "response") #nolint
prediction_pop$deg_inst <- predict(inst_glm, newdata = prediction_pop, type = "response") #nolint
prediction_pop$conc <- predict(casual_conc_glm, newdata = prediction_pop, type = "response") #nolint

observed |>
  ggplot(aes(x = age, y = deg_casual)) +
  facet_wrap(~race + female) +
  geom_point() +
  geom_smooth() +
  geom_line(data = prediction_pop, aes(x = age, y = deg_casual), color = "red")

observed |>
  ggplot(aes(x = age, y = deg_main)) +
  facet_wrap(~race + female) +
  geom_point() +
  geom_smooth() +
  geom_line(data = prediction_pop, aes(x = age, y = deg_main), color = "red")

observed |>
  ggplot(aes(x = age, y = deg_inst)) +
  facet_wrap(~race + female) +
  geom_point() +
  geom_smooth() +
  geom_line(data = prediction_pop, aes(x = age, y = deg_inst), color = "red")

observed |>
  ggplot(aes(x = age, y = conc)) +
  facet_wrap(~race + female) +
  geom_point() +
  geom_smooth() +
  geom_line(data = prediction_pop, aes(x = age, y = conc), color = "red")

prediction_pop |>
  group_by(race) |>
  summarize(mean_casual_conc = mean(conc))

lsvy |>
  group_by(rel) |>

  summarize(mean_absdiff = srvyr::survey_mean(absdiff_sqrt_age, na.rm = TRUE))
