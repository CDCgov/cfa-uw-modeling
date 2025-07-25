#' @title Get Target Mean Degrees for Main and Casual Networks by Age and Race
#' @description Extracts target degrees for main and casual networks from a YAML file. Currently assumes the YAML file
#' includes sexual activity terms under "nodefactor" for each network stratified by a joint set of attributes,
#' age and race, and description of attributes distributions under a sublist named "pop".
#' @param target_yaml_file Path to the YAML file containing network targets.
#' @param nets Character vector specifying which networks to extract targets for. Default is c("main", "casual").
#' @param joint_attrs Character vector specifying the two joint attribute to extract targets for.
#' Default is c("age", "race").
#' @return A data frame where each cell represents the mean degree for that age, race, and network type.
#' @export
get_target_degrees_age_race <- function(target_yaml_file, nets = c("main", "casual"), joint_attrs = c("age", "race")) {
  # load network targets from yaml file
  x <- yaml::read_yaml(target_yaml_file)

  # Validate inputs, currently only supports main and casual networks, age and race joint attributes
  if (nets != c("main", "casual")) {
    stop("Currently only 'main' and 'casual' networks are supported, in that order.")
  }

  if (joint_attrs != c("age", "race")) {
    stop("Currently only race and race are supported as joint attributes, in that order.")
  }

  # Make sure the specified networks and attributes exist in the YAML file
  if (!all(nets %in% names(x))) {
    stop("Specified networks not found in the YAML file.")
  }

  if (length(joint_attrs) != 2) {
    stop("joint_attrs must be a character vector of length 2.")
  }
  if (!attr %in% names(x[[nets[1]]][["nodefactor"]])) {
    stop("Specified attribute not found in the YAML file.")
  }

  joint_name <- paste0(joint_attrs[1], "_", joint_attrs[2])

  if (!joint_name %in% intersect(names(x[[nets[1]]]$nodefactor), names(x[[nets[2]]]$nodefactor))) {
    stop("Joint attribute not found in the YAML file.")
  }

  if (is.null(x$pop)) {
    stop("Population attributes not found in the YAML file.")
  }

  # Extract the attribute values for the specified joint attributes
  ages <- x$pop$age$min:x$pop$age$max
  races <- x$pop$race$levels


  dat <- data.frame(
    main = x[[nets[1]]]$nodefactor[[joint_name]],
    casual = x[[nets[1]]]$nodefactor[[joint_name]],
    age = rep(ages, (length(races))),
    race = rep(races, each = length(ages))
  )

  # Summarize mean degrees and reshape to long format
  dat |>
    dplyr::group_by(age, race) |>
    dplyr::summarize(
      main = mean(main),
      casual = mean(casual),
      .groups = "drop"
    ) |>
    dplyr::mutate(data = "targets") |>
    tidyr::pivot_longer(
      cols = c("main", "casual"),
      names_to = "type",
      values_to = "degree"
    )
}

#' @title Get Edges History from Simulation
#' @description Extracts edges history as the number of edges per time step from a simulation object,
#' calculating the difference between the simulated edges and the target edges for both main and casual networks.
#' NOTE: Assumes edge history is tracked in the simulation object in the "epi" list, not via setting the control
#' parameter 'nwstats' TRUE.
#' @inheritParams get_target_degrees_age_race
#' @param sim A simulation object of class `EpiModel::netsim`.
#' @return A data frame with time, simulation identifier, edges for main and casual networks,
#' the difference from target edges, and the percentage difference from target edges.
#' @importFrom rlang .data
#' @export
get_edges_history <- function(sim, nets = c("main", "casual")) {
  edges <- paste0("edges_", nets)

  if (!all(edges %in% names(sim$epi))) {
    stop("Edges history not found in the simulation object.")
  }

  # Extract edges target for given network
  # always stored as first element in target.stats for each network
  # assumes main network is first in the list
  target_main <- sim$nwparam[[1]]$target.stats[1]
  target_casual <- sim$nwparam[[2]]$target.stats[1]

  # Extract edges history from simulation object
  sim |>
    as.data.frame() |>
    dplyr::select(.data$time, .data$sim, dplyr::all_of(edges)) |>
    dplyr::mutate(
      main_diff = edges[1] - .data$target_main,
      casual_diff = edges[2] - .data$target_casual,
      main_diff_perc = (main_diff) / .data$target_main * 100,
      casual_diff_perc = (casual_diff) / .data$target_casual * 100
    ) |>
    dplyr::group_by(time) |>
    dplyr::mutate(
      mean_main_diff_perc = mean(main_diff_perc, na.rm = TRUE),
      mean_casual_diff_perc = mean(casual_diff_perc, na.rm = TRUE)
    )
}


# frequency of rels by age in networks at end of simulation
summarize_final_degrees <- function(sim) {
  simdat <- NULL

  for (i in seq_len(nsims)) {
    this_sim <- paste0("sim", i)
    d <- data.frame(
      age = floor(sim[["network"]][[this_sim]][[1]] %v% "age"),
      race = sim[["network"]][[this_sim]][[1]] %v% "race",
      deg_main = get_degree(sim[["network"]][[this_sim]][[1]]),
      deg_cas = get_degree(sim[["network"]][[this_sim]][[2]]),
      sim = this_sim
    )
    simdat <- rbind(simdat, d)
  }

  sim_summary <- simdat |>
    dplyr::group_by(sim, age, race) |>
    dplyr::summarize( # calc mean degree by sim, age, race
      main = mean(deg_main),
      casual = mean(deg_cas),
      .groups = "drop"
    ) |>
    tidyr::pivot_longer(
      cols = c("main", "casual"),
      names_to = "type",
      values_to = "deg"
    ) |>
    dplyr::group_by(age, race, type) |>
    dplyr::summarize( # mean and sd of mean degree across sims
      degree = mean(deg),
      IQR1 = quantile(deg, 1 / 4),
      IQR3 = quantile(deg, 3 / 4),
      .groups = "drop"
    ) |>
    dplyr::mutate(data = "simulated")
}

# mean rel durs at end (may not match targets if simulation is not long enough)
get_mean_durations <- function(sim) {
  main_durs <- NULL
  casual_durs <- NULL
  nsims <- sim$control$nsims
  nsteps <- sim$control$nsteps

  for (i in seq_len(nsims)) {
    this_sim <- paste0("sim", i)
    m <- sim[["network"]][[this_sim]][[1]] %n% "lasttoggle"
    mdurs <- nsteps - m[, 3]
    main_durs <- c(main_durs, mean(mdurs, na.rm = TRUE))

    c <- sim[["network"]][[this_sim]][[2]] %n% "lasttoggle"
    cdurs <- nsteps - c[, 3]
    casual_durs <- c(casual_durs, mean(cdurs, na.rm = TRUE))
  }

  data.frame(
    type = c("main", "casual"),
    mean_duration = c(mean(main_durs), mean(casual_durs)),
    sd_duration = c(sd(main_durs), sd(casual_durs))
  )
}
