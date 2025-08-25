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
out$pop$age_group$dist <- c(rep(0.143, 6), 0.142) # even ages with slightly fewer 45-49

out$pop$race$levels <- c("B", "H", "O", "W")
out$pop$race$dist <- c(0.12, 0.19, 0.11, 0.58) # based on 2020 census

# Rake surveys to population margins ---------------------------
us_pop <- 144764299 # from census
pop_agegroup <- data.frame(age_group = out$pop$age_group$levels, Freq = out$pop$age_group$dist * us_pop)
pop_race <- data.frame(race = out$pop$race$levels, Freq = out$pop$race$dist * us_pop)

wsvy <- srvyr::as_survey_design(w, ids = secu, strata = sest, weights = weight, nest = TRUE)
lsvy <- srvyr::as_survey_design(l, ids = secu, strata = sest, weights = weight, nest = TRUE)

wsvy <- survey::rake(wsvy, list(~age_group, ~race), list(pop_agegroup, pop_race))
lsvy <- survey::rake(lsvy, list(~age_group, ~race), list(pop_agegroup, pop_race))

# Add age squared and cross-network ties ---------------------------
wsvy <- wsvy |> dplyr::mutate(
  agesq = age^2,
  agecb = age^3,
  cross_network = ifelse(deg_main > 0 & deg_casual > 0, 1, 0),
  conc = ifelse(deg_casual > 1, 1, 0)
)

lsvy <- lsvy |> dplyr::mutate(
  absdiff_sqrt_age = abs(sqrt(age) - sqrt(alter_age))
)

# Main Network Parameters -----------------------------------
## Edges
mod <- survey::svyglm(deg_main ~ 1, design = wsvy, family = quasipoisson())
pred <- exp(coef(mod)[[1]])
out$main$edges <- pred

## Older partners
mod <- survey::svyglm(olderpartner ~ age, design = wsvy, family = quasipoisson())
pred <- predict(mod, newdata = data.frame(age = 15:49), type = "response")
out$main$olderpartner <- as.numeric(pred)

## Cross-network ties (main-casual)
mod <- survey::svyglm(cross_network ~ 1, design = wsvy, family = quasipoisson())
pred <- exp(coef(mod)[[1]])
out$main$cross_network <- pred
out$casual$cross_network <- pred

## Nodefactor (activity by age group and race)
mod <- survey::svyglm(deg_main ~ age_group, design = wsvy, family = quasipoisson())
pred <- predict(mod, newdata = data.frame(age_group = out$pop$age_group$levels), type = "response")
out$main$nodefactor$age_group <- as.numeric(pred)

mod <- survey::svyglm(deg_main ~ race, design = wsvy, family = quasipoisson())
pred <- predict(mod, newdata = data.frame(race = out$pop$race$levels), type = "response")
out$main$nodefactor$race <- as.numeric(pred)

mod <- survey::svyglm(deg_main ~ age + agesq + race, design = wsvy, family = quasipoisson())
pred <- predict(mod, newdata = data.frame(
  age = rep(15:49, 4),
  agesq = (rep(15:49, 4))**2,
  race = rep(out$pop$race$levels, each = length(15:49))
), type = "response")
out$main$nodefactor$age_race <- as.numeric(pred)

## Nodecov (mean of combined age and agesq across all rels)
mod <- survey::svyglm(comb_age ~ 1, design = lsvy, subset = rel2 == "Marriage/Cohab", family = quasipoisson())
pred <- exp(coef(mod)[[1]])
out$main$nodecov$age <- pred

mod <- survey::svyglm(comb_agesq ~ 1, design = lsvy, subset = rel2 == "Marriage/Cohab", family = quasipoisson())
pred <- exp(coef(mod)[[1]])
out$main$nodecov$agesq <- pred

# Nodematch (by race, and overall)
mod <- survey::svyglm(race_match ~ race, design = lsvy, subset = rel2 == "Marriage/Cohab", family = quasipoisson())
pred <- predict(mod, newdata = data.frame(race = out$pop$race$levels), type = "response")
out$main$nodematch$race <- as.numeric(pred)

# Absdiff (mean abs val in sqrt of age between partners)
mod <- survey::svyglm(absdiff_sqrt_age ~ 1,
  design = lsvy,
  subset = rel2 == "Marriage/Cohab", family = quasipoisson()
)
pred <- exp(coef(mod)[[1]])
out$main$absdiff_sqrt_age <- pred

# Casual Relationships Parameters -------------------
## Edges
mod <- survey::svyglm(deg_casual ~ 1, design = wsvy, family = quasipoisson())
pred <- exp(coef(mod)[[1]])
out$casual$edges <- pred

## Concurrency
mod <- survey::svyglm(conc ~ 1, design = wsvy, family = quasipoisson())
pred <- exp(coef(mod)[[1]])
out$casual$concurrent <- pred

## Nodefactor (activity by age and race)
mod <- survey::svyglm(deg_casual ~ age + agesq + agecb + race, design = wsvy, family = quasipoisson())
pred <- predict(mod, newdata = data.frame(
  age = rep(15:49, 4),
  agesq = (rep(15:49, 4))**2,
  agecb = (rep(15:49, 4))**3,
  race = rep(out$pop$race$levels, each = length(15:49))
), type = "response")
out$casual$nodefactor$age_race <- as.numeric(pred)

## Nodematch (by race)
mod <- survey::svyglm(race_match ~ race, design = lsvy, subset = rel2 == "Casual/Other", family = quasipoisson())
pred <- predict(mod, newdata = data.frame(race = out$pop$race$levels), type = "response")
out$casual$nodematch$race <- as.numeric(pred)

mod <- survey::svyglm(agegrp_match ~ age_group, design = lsvy, subset = rel2 == "Casual/Other", family = quasipoisson())
pred <- predict(mod, newdata = data.frame(age_group = out$pop$age_group$levels), type = "response")
out$casual$nodematch$age_group <- as.numeric(pred)

## Absdiff (mean abs val in sqrt of age between partners)
mod <- survey::svyglm(absdiff_sqrt_age ~ 1,
  design = lsvy,
  subset = rel2 == "Casual/Other" & curr == 1, family = quasipoisson()
)
pred <- exp(coef(mod)[[1]])
out$casual$absdiff_sqrt_age <- pred

## Nodecov (mean of combined age and agesq across all rels)
mod <- survey::svyglm(comb_age ~ 1, design = lsvy, subset = rel2 == "Casual/Other" & curr == 1, family = quasipoisson())
pred <- exp(coef(mod)[[1]])
out$casual$nodecov$age <- pred

mod <- survey::svyglm(comb_agesq ~ 1,
  design = lsvy, subset = rel2 == "Casual/Other" & curr == 1,
  family = quasipoisson()
)
pred <- exp(coef(mod)[[1]])
out$casual$nodecov$agesq <- pred

# Instantaneous Relationships Parameters -------------------
out$inst$summary_time <- "year"
## Edges
mod <- survey::svyglm(deg_inst_high ~ 1, design = wsvy, family = quasipoisson())
pred <- exp(coef(mod)[[1]])
out$inst$edges <- pred

## Nodefactor (activity by age + race)
mod <- survey::svyglm(deg_inst_high ~ age + race, design = wsvy, family = quasipoisson())
pred <- predict(mod, newdata = data.frame(
  age = rep(15:49, 4),
  race = rep(out$pop$race$levels, each = length(15:49))
), type = "response")
out$inst$nodefactor$age_race <- as.numeric(pred)

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

main_race_mixmat <- epimodelcfa::matrix_symmetrical(mainprops)
casual_race_mixmat <- epimodelcfa::matrix_symmetrical(casprops)

out$main$mixmat$race <- main_race_mixmat
out$casual$mixmat$race <- casual_race_mixmat

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

out$main$mixmat$age_group <- main_ag_mixmat
out$casual$mixmat$age_group <- casual_ag_mixmat


# save out as yaml
params_name <- "nw_params.yaml"
yaml::write_yaml(out, file.path(param_folder, params_name))
