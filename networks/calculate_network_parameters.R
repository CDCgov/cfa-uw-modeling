# Here we set up list with information about desired population characteristics
# (and save as yaml config file)
# Will be used to generate starting network and ergm target statistics
# for network fitting process

# 1) set up list with population parameters
# And enough information to parameterize starting nodal attributes
# requires some assumptions and previous knowledge of data

# nolint start
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

# 2) load data from NSFG (long and wide versions)
# and convert to survey objects to use weights
w <- readRDS(here::here("data", "nsfg_wide.rds"))
wsvy <- srvyr::as_survey_design(w, weights = weight)

l <- readRDS(here::here("data", "nsfg_long.rds"))
lsvy <- srvyr::as_survey_design(l, weights = weight)

# 3) calculate population metrics of interest -----------------------------
## A. Proportion of respondents in a relationship, by main / casual (matrix)
## and overall prop for inst
## B. Proportion of respondents in a relationship, by type and age group
## C. Proportion of respondents in a relationship, by type and age group and female
## D. Proportion of respondents in a relationship, by type and race


## A.1. Proportion of respondents in a relationship, by main / casual (matrix)
degmat <- wsvy |>
  dplyr::group_by(deg_main, deg_casual) |>
  dplyr::summarize(sum = srvyr::survey_total(vartype = NULL)) |>
  dplyr::ungroup() |>
  dplyr::mutate(prop = sum / sum(sum)) |>
  dplyr::select(-sum) |>
  tidyr::pivot_wider(names_from = deg_casual, values_from = prop) |>
  dplyr::select(-deg_main) |>
  as.matrix()

### add overall concurrency types to out list
out$main$has_casual <- degmat[2, 2] + degmat[2, 3]
out$casual$has_main <- degmat[2, 2] + degmat[2, 3]
out$casual$concurrent <- sum(colSums(degmat, na.rm = TRUE)[3:4])

## A.2. Proportion of respondents with an inst rels
inst <- wsvy |>
  dplyr::group_by(deg_inst_high) |>
  dplyr::summarize(prop = srvyr::survey_prop(vartype = NULL))

inst_activity_year <-
  inst$prop[2] + 2 * inst$prop[3] + 3 * inst$prop[4] +
  4 * inst$prop[5] + 5 * inst$prop[6] + 7 * inst$prop[8]

### add to out list
out$inst$density <- inst_activity_year
out$inst$summary_time <- "year"

## B. For nodefactor terms
# Proportion of respondents in a relationship
### (sum of all rels in each group / pop)
### accounts for casual network having >1 active rel, and >1 inst in last year

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


## Activity (nodefactor) by Age, Race
out$main$nodefactor$age_race <- unname(unlist(
  target_props(wsvy, attr = age, attr2 = race, degtype = deg_main)$prop
))
out$casual$nodefactor$age_race <- unname(unlist(
  target_props(wsvy, attr = age, attr2 = race, degtype = deg_casual)$prop
))
out$inst$nodefactor$age_race <- unname(unlist(
  target_props(wsvy, attr = age, attr2 = race, degtype = deg_inst_high)$prop
))

## Activity (nodefactor) by Age Group
out$main$nodefactor$age_group <- unname(unlist(target_props(wsvy, attr = age_group, degtype = deg_main)[, 2]))
out$casual$nodefactor$age_group <- unname(unlist(target_props(wsvy, attr = age_group, degtype = deg_casual)[, 2]))
out$inst$nodefactor$age_group <- unname(unlist(target_props(wsvy, attr = age_group, degtype = deg_inst_high)[, 2]))

## Activity (nodefactor) by Race
out$main$nodefactor$race <- unname(unlist(
  target_props(wsvy, attr = race, degtype = deg_main)[, 2]
))
out$casual$nodefactor$race <- unname(unlist(
  target_props(wsvy, attr = race, degtype = deg_casual)[, 2]
))
out$inst$nodefactor$race <- unname(unlist(
  target_props(wsvy, attr = race, degtype = deg_inst_high)[, 2]
))

## Activity (nodefactor) by Age Group and Race
out$main$nodefactor$age_group_race <-
  target_props(wsvy, attr = age_group, attr2 = race, degtype = deg_main)$prop

out$casual$nodefactor$age_group_race <-
  target_props(wsvy, attr = age_group, attr2 = race, degtype = deg_casual)$prop

out$inst$nodefactor$age_group_race <-
  target_props(wsvy, attr = age_group, attr2 = race, degtype = deg_inst_high)$prop


# Dyad-based metrics of interest
# 1. absdiff sqrt age (by relationship type)
# 2. race mixing

## Absdiff in sqrt of age between partners
# with female age adjusted to account for asymmmetry (they are usually younger than male partners)
adjs <- lsvy |>
  dplyr::summarize(mean = srvyr::survey_mean(asym_agediff, na.rm = TRUE, vartype = NULL))

lsvy <- lsvy |>
  dplyr::mutate(
    age_adj1 = ifelse(female == 1, age + 1.78, age),
    age_adj2 = ifelse(female == 0, alter_age + 1.78, alter_age)
  )

agemix_vars <- lsvy |>
  dplyr::mutate(diff_sqrt_age = abs(sqrt(age_adj1) - sqrt(age_adj2))) |>
  dplyr::group_by(rel2) |>
  dplyr::summarize(abs_diff_sqrt_age = srvyr::survey_mean(diff_sqrt_age,
    na.rm = TRUE,
    vartype = NULL
  ))
out$main$absdiff_sqrt_age <- as.numeric(agemix_vars[1, 2])
out$casual$absdiff_sqrt_age <- as.numeric(agemix_vars[2, 2])
out$pop$age$female_age_adj <- 1.78

# no information about age mixing for non-current (one-time) partners
# use current casuals as proxy
out$inst$absdiff_sqrt_age <- as.numeric(agemix_vars[2, 2])


## Race Matching aka nodematch
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

### duration
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


# save out as yaml
yaml::write_yaml(out, here::here("params", "nw_params.yaml"))

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
