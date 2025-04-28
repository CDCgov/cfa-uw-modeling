# nolint start: object_name_linter

#' @title Epidemic Model Parameters (Mgen)
#'
#' @description Sets the epidemic parameters for stochastic network models
#'              simulated with [`netsim`]
#'
#' @param ... Additional arguments passed to the function.
#'
#' @return
#' A list object of class `param.net`, which can be passed to EpiModel function [`netsim`].
#'
#' @export
#'
param_mgen <- function(
    # placeholder sexual act rates
    act_rate_main = 5.63 / 7,
    act_rate_casual = 1.64 / 7,
    act_rate_instant = 1,
    # placeholder natural history
    latent_period = 7,
    infectious_period = 30 * 4,
    screening_rate_females = 1 / 365,
    screening_rate_males = NA,
    vital = FALSE,
    ...) {
  p <- get_args(
    formal.args = formals(sys.function()),
    dot.args = list(...)
  )

  class(p) <- "param.net"
  return(p)
}

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


#' @title Epidemic Model Control Settings
#'
#' @description Sets the controls for stochastic network models simulated with [`netsim`].
#'
#' @param simno Unique ID for the simulation run, used for file naming purposes if used
#'        with the `EpiModelHPC` package.
#' @param nsims Number of simulations to run overall.
#' @param ncores Number of parallel cores per run, if parallel processing is used.
#' @param nsteps Number of time steps per simulation.
#' @param start Starting time step for simulation, with default to `1` to run new simulation. This
#'        may also be set to 1 greater than the final time step of a previous simulation to resume
#'        the simulation with different parameters.
#' @param cumulative.edgelist If `TRUE`, models tracks a cumulative edgelist (a historical list of
#'        all partners within each network over time), useful for partner services models.
#' @param truncate.el.cuml Number of time steps of the cumulative edgelist to retain. See help file
#'        for `EpiModel::update_cumulative_edgelist` for details
#' @param initialize.FUN Module function to use for initialization of the epidemic model.
#' @param aging.FUN Module function for aging.
#' @param departure.FUN Module function for general and disease-related departures.
#' @param arrival.FUN Module function for entries into the sexually active population.
#' @param partident.FUN Module function for partner identification process.
#' @param resim_nets.FUN Module function for network resimulation at each time
#'        step.
#' @param summary_nets.FUN Module function for network statistic extraction at
#'        each time step.
#' @param acts.FUN Module function to simulate the number of sexual acts within
#'        partnerships.
#' @param condoms.FUN Module function to simulate condom use within acts.
#' @param prev.FUN Module function to calculate prevalence summary statistics.
#' @param verbose.FUN Module function to print model progress to the console or
#'        external text files.
#' @param cleanup.FUN Module function to tidy the `netsim_dat` object at the end
#'        of each step.
#' @param save.nwstats Calculate and save network statistics as defined in the
#'        `simnet` modules.
#' @param nwstats.formula Right-hand side formula for network statistics to
#'        monitor, with default `"formation"` equal to the model formula.
#'        Supports \code{\link{multilayer}} specification.
#' @param tergmLite Logical indicating usage of `tergmLite` (non-modifiable
#'        currently for `EpiModelHIV`).
#' @param save.network If `TRUE`, networkLite object is saved at simulation end.
#' @param tergmLite.track.duration If `TRUE`, track duration information
#'        for models in `tergmLite` simulations. Supports
#'        \code{\link{multilayer}} specification.
#' @param save.other A character vector of elements on the main `netsim_dat`
#'        data object to save out after each simulation.
#' @param verbose If `TRUE`, print out simulation progress to the console
#'        if in interactive mode or text files if in batch mode.
#' @param set.control.ergm Control arguments passed to `ergm`'s
#'        `simulate_formula.network`. Supports \code{\link{multilayer}}
#'        specification.
#' @param set.control.tergm Control arguments passed to `tergm`'s
#'        `simulate_formula.network`. Supports \code{\link{multilayer}}
#'        specification.
#' @param dat.updates See the general documentation for `dat.updates` at [`EpiModel::control.net`].
#'        The standard implementation for the current MSM model is found in
#'        [`standard_resimnet_updates`].
#' @param ... Additional arguments passed to the function.
#'
#' @return
#' A list object of class `control.net`, which can be passed to the EpiModel function [`netsim`].
#'
#' @export
#'
control_mgen <- function(
    simno = 1,
    nsims = 1,
    ncores = 1,
    nsteps = 10,
    start = 1,
    cumulative.edgelist = FALSE,
    save.nwstats = TRUE,
    nwstats.formula = "formation",
    tergmLite = TRUE,
    truncate.el.cuml = 0,
    # modules
    # amr.FUN = amr_mgen,
    # initialize.FUN = initialize_mgen,
    initialize.FUN = initialize.net,
    # infection.FUN = transmit_mgen,
    infection.FUN = infection.net,
    recovery.FUN = recovery.net,
    # progress.FUN = progress_mgen,
    # prevalence.FUN = trackers_mgen,
    prevalence.FUN = prevalence.net,
    resim_nets.FUN = resim_nets,
    # resim_nets.FUN = resim_nets_mgen
    summary_nets.FUN = summary_nets,
    # tx.FUN = tx_mgen,
    verbose.FUN = verbose.net,
    # cleanup.FUN = cleanup_msm,
    # end modules
    save.network = FALSE,
    tergmLite.track.duration = FALSE,
    set.control.ergm = control.simulate.formula(
      MCMC.burnin = 2e5
    ),
    set.control.tergm = control.simulate.formula.tergm(
      MCMC.burnin.min = 5000
    ),
    save.other = c("attr", "temp", "el", "net_attr"),
    verbose = FALSE,
    dat.updates = resimnet_updates_sti,
    ...) {
  formal.args <- formals(sys.function())
  dot.args <- list(...)
  p <- get_args(formal.args, dot.args)

  p$skip.check <- TRUE
  p$save.transmat <- FALSE

  bi.mods <- grep(".FUN", names(formal.args), value = TRUE)
  bi.mods <- bi.mods[which(sapply(bi.mods, function(x) !is.null(eval(parse(text = x))),
    USE.NAMES = FALSE
  ) == TRUE)]
  p$bi.mods <- bi.mods
  p$user.mods <- grep(".FUN", names(dot.args), value = TRUE)
  p[["f.names"]] <- c(p[["bi.mods"]], p[["user.mods"]])

  if (is.null(p$verbose.int)) {
    p$verbose.int <- 1
  }

  p$resimulate.network <- TRUE
  p$save.diss.stats <- FALSE

  p <- set.control.class("control.net", p)
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
resimnet_updates_sti <- function(dat, at, network) {
  if (network == 0L) {
    dat <- set_attr(dat, "deg_casual", EpiModel::get_degree(dat$run$el[[2]]))
  } else if (network == 1L) {
    dat <- set_attr(dat, "deg_main", EpiModel::get_degree(dat$run$el[[1]]))
  } else if (network == 2L) {
    dat <- set_attr(dat, "deg_casual", EpiModel::get_degree(dat$run$el[[2]]))
    dat <- set_attr(
      dat, "deg_tot",
      get_attr(dat, "deg_main") + EpiModel::get_degree(dat$run$el[[2]])
    )
  }
  return(dat)
}

# nolint end
