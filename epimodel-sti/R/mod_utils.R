#' @title Updates to Network Attributes During Mpox Resimulation
#'
#' @description This function is updates the network-related degree attributes on the three-layer
#' MSM sexual network that account for the resimulated network structure.
#'
#' @inheritParams mod_aging
#' @param network Integer for network number (values 1, 2, or 3 for main, casual, and one-time
#'                networks).
#'
#' @details
#' This function is called between network resimulations in [`EpiModel::resim_nets`], passed into
#' [`EpiModel::control.net`] through the `dat.updates` argument. This implementation updates degree attributes
#' calculated from the current network snapshot for use as ERGM terms in the other network layers
#' (e.g., degree in the casual network is a function of the degree in the main network). See the
#' general documentation for `dat.updates` at [`EpiModel::control.net`].
#'
#'
#' @export
#'
resimnet_updates_sti <- function(dat, at, network) {
  if (network == 0L) {
    dat <- EpiModel::set_attr(dat, "deg_casual", EpiModel::get_degree(dat$run$el[[2]]))
  } else if (network == 1L) {
    dat <- EpiModel::set_attr(dat, "deg_main", EpiModel::get_degree(dat$run$el[[1]]))
  } else if (network == 2L) {
    dat <- EpiModel::set_attr(dat, "deg_casual", EpiModel::get_degree(dat$run$el[[2]]))
    dat <- EpiModel::set_attr(
      dat, "deg_tot",
      EpiModel::get_attr(dat, "deg_main") + EpiModel::get_degree(dat$run$el[[2]])
    )
  }
  dat
}
