#' @title Network Resimulation for Time Steps 2+
#'
#' @description Resimulation of dynamic networks
#'
#' @inheritParams initialize_mgen
#'
#' @export

resim_nets_mgen <- function(dat, at) {
  # nolint start: cyclocomp_linter, object_name_linter
  # Calculate active attribute
  active <- get_attr(dat, "active")
  idsActive <- which(active == 1)
  anyActive <- ifelse(length(idsActive) > 0, TRUE, FALSE)
  if (dat$param$groups == 2) {
    group <- get_attr(dat, "group")
    groupids.1 <- which(group == 1) # nolint
    groupids.2 <- which(group == 2) # nolint
    nActiveG1 <- length(intersect(groupids.1, idsActive))
    nActiveG2 <- length(intersect(groupids.2, idsActive))
    anyActive <- ifelse(nActiveG1 > 0 & nActiveG2 > 0, TRUE, FALSE)
  }

  # Network resimulation, with dat.updates interspersed
  if (anyActive == TRUE && get_control(dat, "resimulate.network") == TRUE) {
    ## Edges Correction
    dat <- edges_correct_mpox(dat, at)

    ## network resimulation
    dat.updates <- NVL(get_control(dat, "dat.updates"), function(dat, ...) dat) # nolint
    dat <- dat.updates(dat = dat, at = at, network = 0L)

    # only update all networks every week, otherwise just inst
    # (and don't need dat.updates, which updates main/cas degree attrs)
    if (at %% 7 == 0) {
      for (network in seq_len(dat$num.nw)) {
        dat <- simulate_dat(dat = dat, at = at, network = network)
        dat <- dat.updates(dat = dat, at = at, network = network)
      }
    } else {
      dat <- simulate_dat(dat = dat, at = at, network = 3)
    }
  }

  # Cummulative edgelist
  truncate.el.cuml <- get_control(dat, "truncate.el.cuml")
  for (network in seq_len(dat$num.nw)) {
    dat <- update_cumulative_edgelist(dat, network, truncate.el.cuml)
  }

  return(dat)

  # nolint end
}
