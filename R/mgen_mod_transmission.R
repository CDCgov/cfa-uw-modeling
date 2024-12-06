#' @title Model Mgen Infection between Partners
#'
#' @description mgen transmission
#'
#' @param $$
#'
#' @export
#'

transmit_mgen <- function(dat, at) {
  # Variables ---------------------------------------------------------------
  active <- get_attr(dat, "active")
  infTime <- get_attr(dat, "infTime")
  status <- get_attr(dat, "status")
  group <- get_attr(dat, "group")

  inf.prob <- get_param(dat, "inf.prob")
  inf.prob.g2 <- get_param(dat, "inf.prob.g2")
  act.rate <- get_param(dat, "act.rate")

  # Vector of infected and susceptible IDs
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)

  # Initialize vectors
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
      linf.prob <- length(inf.prob)
      if (is.null(inf.prob.g2)) {
        del$transProb <- ifelse(del$infDur <= linf.prob,
          inf.prob[del$infDur],
          inf.prob[linf.prob]
        )
      } else {
        # FLAG
        del$transProb <- ifelse(group[del$sus] == 1,
          ifelse(del$infDur <= linf.prob,
            inf.prob[del$infDur],
            inf.prob[linf.prob]
          ),
          ifelse(del$infDur <= linf.prob,
            inf.prob.g2[del$infDur],
            inf.prob.g2[linf.prob]
          )
        )
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
      transmit <- rbinom(nrow(del), 1, del$finalProb)
      del <- del[which(transmit == 1), ]

      # Set new infections vector
      idsNewInf <- unique(del$sus)
      status[idsNewInf] <- "i"
      dat <- set_attr(dat, "status", status)
      infTime[idsNewInf] <- at
      dat <- set_attr(dat, "infTime", infTime)
      nInf <- sum(group[idsNewInf] == 1)
      nInfG2 <- sum(group[idsNewInf] == 2)
      totInf <- nInf + nInfG2
    } # end some discordant edges condition
  } # end some active discordant nodes condition


  # Output ------------------------------------------------------------------

  # Save transmission matrix
  if (totInf > 0) {
    dat <- set_transmat(dat, del, at)
  }

  ## Save incidence vector
  dat <- set_epi(dat, "si.flow", at, nInf)
  dat <- set_epi(dat, "si.flow.g2", at, nInfG2)

  return(dat)
}
