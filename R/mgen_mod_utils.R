# nolint start: object_name_linter

#' @title Epidemic Model Initial Conditions
#'
#' @description Sets the initial conditions for a stochastic epidemic models
#'              simulated with [`netsim`].
#'
#' @return
#' A list object of class `init.net`, which can be passed to EpiModel function [`netsim`].
#'
#' @export
#'
init_mgen <- function(ninfs = 5,
                      ...) {
  p <- get_args(
    formal.args = formals(sys.function()),
    dot.args = list(...)
  )

  class(p) <- "init.net"
  return(p)
}

#' @title Updates to Network Attributes During Mpox Resimulation
#'
#' @description This function is updates the network-related degree attributes on the three-layer
#' MSM sexual network that account for the resimulated network structure.
#'
#' @param network Integer for network number (values 1, 2, or 3 for main, casual, and one-time
#'                networks).
#'
#' @details
#' This function is called between network resimulations in [`EpiModel::resim_nets`], passed into
#' [`control_msm`] through the `dat.updates` argument. This implementation updates degree attributes
#' calculated from the current network snapshot for use as ERGM terms in the other network layers
#' (e.g., degree in the casual network is a function of the degree in the main network). See the
#' general documentation for `dat.updates` at [`EpiModel::control.net`].
#'
#'
#' @export
#'
resimnet_updates <- function(dat, at, network) {
  if (network == 0L) {
    dat <- set_attr(dat, "deg.pers", EpiModel::get_degree(dat$el[[2]]))
  } else if (network == 1L) {
    dat <- set_attr(dat, "deg.main", EpiModel::get_degree(dat$el[[1]]))
  } else if (network == 2L) {
    dat <- set_attr(dat, "deg.pers", EpiModel::get_degree(dat$el[[2]]))
    dat <- set_attr(dat, "deg.tot", pmin(
      get_attr(dat, "deg.main") +
        EpiModel::get_degree(dat[["el"]][[2]]),
      3
    ))
  }
  return(dat)
}

# nolint end
