#' @title Model Mgen Infection between Partners
#'
#' @description Disease transmission
#'
#' @inheritParams initialize_mgen
#' @importFrom EpiModel get_attr set_attr set_epi get_param set_transmat discord_edgelist
#'
#' @export

mgen_infection <- function(dat, at) {
  # Notes
  ## for now:
  ## - no age/rel difference in act rates
  ## - no condoms
  ## leaving infection-stage inf prob rate functionality in but unused atm

  # Variables ---------------------------------------------------------------

  active <- get_attr(dat, "active")
  infTime <- get_attr(dat, "infTime")
  status <- get_attr(dat, "status")
  female <- get_attr(dat, "female")

  infProbMTF <- get_param(dat, "infProbMTF")
  infProbFTM <- get_param(dat, "infProbFTM")
  act.rate <- get_param(dat, "act.rate")

  inter.eff <- get_param(dat, "inter.eff")
  inter.start <- get_param(dat, "inter.start")

  # Vector of infected and susceptible IDs
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)

  # Initialize vectors
  # G2 = female postscript (female attr == 1)
  nInf <- nInfG2 <- totInf <- 0

  # Process -----------------------------------------------------------------
  # If some infected AND some susceptible, then proceed
  if (nElig > 0 && nElig < nActive) {
    # Get discordant edgelist
    del_list <- lapply(seq_len(dat$num.nw), discord_edgelist, dat = dat, at = at, include.network = TRUE)
    del <- dplyr::bind_rows(del_list)

    # If some discordant edges, then proceed
    if (NROW(del) > 0) {
      # Infection duration to at
      del$infDur <- at - infTime[del$inf]
      del$infDur[del$infDur == 0] <- 1

      # Calculate infection-stage transmission rates
      linf.prob <- length(infProbMTF)
      if (is.null(infProbFTM)) {
        del$transProb <- ifelse(del$infDur <= linf.prob,
          infProbMTF[del$infDur],
          infProbMTF[linf.prob]
        )
      } else {
        # FLAG
        del$transProb <- ifelse(female[del$sus] == 1,
          ifelse(del$infDur <= linf.prob,
            infProbMTF[del$infDur],
            infProbMTF[linf.prob]
          ),
          ifelse(del$infDur <= linf.prob,
            infProbFTM[del$infDur],
            infProbFTM[linf.prob]
          )
        )
      }

      # Interventions
      if (!is.null(inter.eff) && at >= inter.start) {
        del$transProb <- del$transProb * (1 - inter.eff)
      }

      # Calculate infection-stage act/contact rates
      lact.rate <- length(act.rate)
      del$actRate <- ifelse(del$infDur <= lact.rate,
        act.rate[del$infDur],
        act.rate[lact.rate]
      )

      # Calculate final transmission probability per timestep
      del$finalProb <- 1 - (1 - del$transProb)^del$actRate

      # Randomize transmissions and subset df
      transmit <- stats::rbinom(nrow(del), 1, del$finalProb)
      del <- del[which(transmit == 1), ]

      # Set new infections vector
      idsNewInf <- unique(del$sus)
      status[idsNewInf] <- "i"
      dat <- set_attr(dat, "status", status)
      infTime[idsNewInf] <- at
      dat <- set_attr(dat, "infTime", infTime)
      nInf <- sum(female[idsNewInf] == 1)
      nInfG2 <- sum(female[idsNewInf] == 2)
      totInf <- nInf + nInfG2
    } # end some discordant edges condition
  } # end some active discordant nodes condition


  # Output ------------------------------------------------------------------

  # Save transmission matrix
  if (totInf > 0) {
    dat <- set_transmat(dat, del, at)
  }

  ## Save incidence vector
  dat <- set_epi(dat, "si.flow.female0", at, nInf)
  dat <- set_epi(dat, "si.flow.female1", at, nInfG2)
  return(dat) # nolint
}
