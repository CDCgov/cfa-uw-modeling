# Make "wide" data (cleaned survey data) "long"
# where each row is an ego/relationship pair
# of current rels for marcoh and casual
# and inst partners

# nolint start

source(here::here("data", "nsfg_scripts", "egodata_functions.R"))

w <- readRDS(here::here("data", "nsfg_wide.rds"))

# pull out ego-specific variables
ego_vars <- c(
  "ego", "weight", "wgt2017_2019", "secu", "sest", "female", "race",
  "age", "age_group", "sqrtage", "ageF", "ageM",
  "hadsex", "parts1yr", "parts1yr_capped",
  "deg", "deg_mar", "deg_cohab", "deg_casual",
  "deg_inst", "deg_main", "deg_inst_high",
  "sex4wk", "cond4wk", "sex1wk",
  "p_cond_month", "condom_use", "condom_use2"
)

egos <- w |> select(all_of(ego_vars))

# define alter variables in wide data
alter_vars <- c(
  "ego", "weight", "race", "age_group", "age", "female", # ego chars
  "once1", "rel1", "page1", "curr1", "partdur1", "prace1",
  "once2", "rel2", "page2", "curr2", "partdur2", "prace2",
  "once3", "rel3", "page3", "curr3", "partdur3", "prace3"
)

# filter wide data to only those respondents who have had sex
alters <- w |>
  dplyr::filter(hadsex == 1) |>
  dplyr::select(all_of(alter_vars))
# next function doesn't play nice with tibbles
alters <- as.data.frame(alters)

# make long
a <- reshape_edgelist(alters)

# reverse sex variable
# (currently reflects sex of associated ego, not the alter)
# rename "page" "age", prace, etc
a2 <- a |>
  dplyr::mutate(
    alter_female = ifelse(female == 0, 1, ifelse(female == 1, 0, NA)),
    alter_race = prace,
    alter_age = page,
    rel2 = ifelse(rel %in% 1:2, 1, rel),
    partdur = ifelse((is.na(partdur) & once == 1) | partdur == 0,
      0.5, partdur
    )
  ) |>
  dplyr::mutate(
    rel = factor(rel,
      levels = c(1:3),
      labels = c(
        "Marriage", "Cohabitation",
        "Casual/Other"
      )
    ),
    rel2 = factor(rel2,
      levels = c(1, 3),
      labels = c("Marriage/Cohab", "Casual/Other")
    )
  ) |>
  dplyr::select(-c("page", "prace"))

# filter to current marcoh / casual
a3 <- a2 |> dplyr::filter(curr == 1)

# For best ergm estimates, we need to input alter_age (missing a lot of under19 and casual alter_ages)
partner_age_glm <- glm(alter_age ~ age + race + female + rel2, data = a3, weights = weight)
agemiss <- which(is.na(a3$alter_age))

testdat <- a3[agemiss, ]
a3$alter_age[agemiss] <- predict(partner_age_glm, newdata = testdat, type = "response")


# Return to cleaning

a4 <- a3 |>
  dplyr::mutate(alter_age = ifelse(alter_age == 50, 49.9, alter_age)) |>
  dplyr::mutate(alter_age = ifelse(alter_age > 50, NA, alter_age)) |>
  dplyr::mutate(alter_age = ifelse(alter_age < 15, 15, alter_age)) |>
  dplyr::mutate( # race cat
    alter_race = ifelse(alter_race == 1, "H",
      ifelse(alter_race == 2, "W",
        ifelse(alter_race == 3, "B",
          ifelse(alter_race == 4, "O", NA)
        )
      )
    )
  )

a4$alter_age_group <- cut(round(a4$alter_age), 7)

a5 <- a4 |>
  dplyr::mutate(
    agegrp_match = ifelse(age_group == alter_age_group, TRUE, FALSE),
    race_match = ifelse(race == alter_race, TRUE, FALSE)
  ) |>
  dplyr::mutate(
    asym_agediff = ifelse(female == 1,
      alter_age - age,
      age - alter_age
    ),
    asym_agediff_sqrt = ifelse(female == 1,
      sqrt(alter_age) - sqrt(age),
      sqrt(age) - sqrt(alter_age)
    ),
    comb_age = age + alter_age,
    comb_agesq = age^2 + alter_age^2
  )


# finally, remove any remaining rows with NA reponses in key variables
a6 <- a5 |>
  dplyr::filter(!is.na(alter_age), !is.na(rel2), !is.na(alter_race))

# save out
saveRDS(a6, here::here("data", "nsfg_long.rds"))

# nolint end
