#' @title Modules for vital dynamics
#'
#' @description Handles node aging, departure, and arrivals
#'
#' @inheritParams initialize_mgen
#' @importFrom EpiModel get_attr set_attr get_epi set_epi get_param append_core_attr append_attr
#' get_edgelist apportion_lr
#'
#' @export

mod_aging <- function(dat, at) {
  # Calc Updated Age Attributes
  age <- get_attr(dat, "age")
  units <- get_param(dat, "units_per_year")
  age <- age + (1 / units)
  age_group <- dplyr::case_when(
    age < 20 ~ 1,
    age >= 20 & age < 25 ~ 2,
    age >= 25 & age < 30 ~ 3,
    age >= 30 & age < 35 ~ 4,
    age >= 35 & age < 40 ~ 5,
    age >= 40 & age < 45 ~ 6,
    age >= 45 ~ 7
  )

  # Update Attributes
  dat <- set_attr(dat, "age", age)
  dat <- set_attr(dat, "agesq", age^2)
  dat <- set_attr(dat, "age_group", age_group)

  ## Summary statistics ##
  dat <- set_epi(dat, "meanAge", at, mean(age, na.rm = TRUE))
  dat
}

# Departures Module ----------------------------------------------------

mod_departures <- function(dat, at) {
  ## Attributes
  active <- get_attr(dat, "active")
  exitTime <- get_attr(dat, "exitTime")
  age <- get_attr(dat, "age")
  status <- get_attr(dat, "status")

  ## Parameters
  exitAge <- get_param(dat, "exitAge")
  ## if we had ASMR we would add that here

  ## Query alive but past simulation age range
  ## this setup a little odd make it easier to include ASMR later
  idsElig <- which(active == 1 & ceiling(age) >= exitAge + 1)
  nElig <- length(idsElig)
  nDepts <- 0

  if (nElig > 0) {
    idsDept <- idsElig
    nDepts <- length(idsDept)
    ## Update nodal attributes
    active[idsDept] <- 0
    exitTime[idsDept] <- at

    ## If departing nodes are in main relationship,
    ## Set partner's olderpartner attr to 1
    el <- get_edgelist(dat = dat, network = 1) # get edgelist
    rels1 <- which(el[, 1] %in% idsDept) # any departing id in row 1
    rels2 <- which(el[, 2] %in% idsDept) # any departing id in row 2
    allrels <- c(rels1, rels2) # combine

    if (length(allrels) > 0) {
      allParts <- el[allrels, ] # get all IDs (departing and partner)
      idsParts <- setdiff(allParts, idsDept) # extract partner IDs
      # Update Partner Attribute
      olderpartner <- get_attr(dat, "olderpartner")
      olderpartner[idsParts] <- 1
      dat <- set_attr(dat, "olderpartner", olderpartner)
    }
  }

  ## Reset attr
  dat <- set_attr(dat, "active", active)
  dat <- set_attr(dat, "exitTime", exitTime)

  ## Summary statistics
  dat <- set_epi(dat, "d.flow", at, nDepts)
  dat <- set_epi(dat, "edges_main", at, nrow(dat$run$el[[1]]))
  dat <- set_epi(dat, "edges_casual", at, nrow(dat$run$el[[2]]))
  dat <- set_epi(dat, "edges_inst", at, nrow(dat$run$el[[3]]))

  dat
}


# Arrivals Module ----------------------------------------------------

mod_arrivals <- function(dat, at) {
  ## Parameters
  n <- sum(get_attr(dat, "active") == 1)
  aType <- get_param(dat, "arrivalType")
  entryAge <- get_param(dat, "entryAge")
  femaleAgeAdj <- get_param(dat, "entryFemaleAgeAdj")
  femaleProb <- get_param(dat, "entryFemaleProb")
  raceNames <- get_param(dat, "entryRaceNames")
  raceProbs <- get_param(dat, "entryRaceProbs")

  nArrivals <- 0

  if (!aType %in% c("rate", "departures")) {
    stop("Arrival Type must be either 'rate' or 'departures'")
  }

  if (aType == "rate") {
    a_rate <- get_param(dat, "arrival.rate")

    ## Process
    nArrivalsExp <- n * a_rate
    nArrivals <- stats::rpois(1, nArrivalsExp)
  }

  if (aType == "departures") {
    nArrivals <- get_epi(dat, "d.flow", at)
  }

  if (nArrivals > 0) {
    ## Determine sex, race
    arrivalSex <- stats::rbinom(nArrivals, 1, femaleProb)
    arrivalRace <- apportion_lr(nArrivals, raceNames, raceProbs)

    ## Update attributes
    dat <- append_core_attr(dat, at, nArrivals)
    dat <- append_attr(dat, "status", "s", nArrivals)
    dat <- append_attr(dat, "infTime", NA, nArrivals)
    dat <- append_attr(dat, "age", entryAge, nArrivals)
    dat <- append_attr(dat, "agesq", entryAge^2, nArrivals)
    dat <- append_attr(dat, "age_group", 1, nArrivals)
    dat <- append_attr(dat, "race", arrivalRace, nArrivals)
    dat <- append_attr(dat, "female", arrivalSex, nArrivals)
    dat <- append_attr(dat, "olderpartner", 0, nArrivals)
  }

  ## Summary statistics
  dat <- set_epi(dat, "a.flow", at, nArrivals)

  dat
}
