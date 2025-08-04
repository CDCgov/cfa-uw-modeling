#' @title Get Target Mean Degrees for Main and Casual Networks by Age and Race
#' @description Extracts target degrees for main and casual networks from a YAML file. Currently assumes the YAML file
#' includes sexual activity terms under "nodefactor" for each network stratified by a joint set of attributes,
#' age and race, and description of attributes distributions under a sublist named "pop".
#' @param target_yaml_file Path to the YAML file containing network targets.
#' @param nets Character vector specifying which networks to extract targets for. Default is c("main", "casual").
#' @param joint_attrs Character vector specifying the two joint attribute to extract targets for.
#' Default is c("age", "race").
#' @return A data frame where each cell represents the mean degree for that age, race, and network type.
#' @importFrom yaml read_yaml
#' @importFrom dplyr mutate group_by summarize
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data
#' @export
get_target_degrees_age_race <- function(yaml_params_loc, nets = c("main", "casual"), joint_attrs = c("age", "race")) {
  # load network targets from yaml file
  x <- read_yaml(yaml_params_loc)

  # Validate inputs, currently only supports main and casual networks, age and race joint attributes
  if (sum(nets == c("main", "casual")) != 2) {
    stop("Currently only 'main' and 'casual' networks are supported, in that order.")
  }

  if (sum(joint_attrs == c("age", "race")) != 2) {
    stop("Currently only race and race are supported as joint attributes, in that order.")
  }

  # Make sure the specified networks and attributes exist in the YAML file
  if (!all(nets %in% names(x))) {
    stop("Specified networks not found in the YAML file.")
  }

  if (length(joint_attrs) != 2) {
    stop("joint_attrs must be a character vector of length 2.")
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
    casual = x[[nets[2]]]$nodefactor[[joint_name]],
    age = rep(ages, (length(races))),
    race = rep(races, each = length(ages))
  )

  # Summarize mean degrees and reshape to long format
  dat |>
    group_by(.data$age, .data$race) |>
    summarize(
      main = mean(.data$main),
      casual = mean(.data$casual),
      .groups = "drop"
    ) |>
    mutate(data = "targets") |>
    pivot_longer(
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
#' @importFrom dplyr select rename_with mutate filter group_by ungroup all_of
#' @importFrom tidyr pivot_longer
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
    select(.data$time, .data$sim, all_of(edges)) |>
    rename_with(~ gsub("edges_", "", .), all_of(edges)) |>
    pivot_longer(cols = all_of(nets), names_to = "net", values_to = "edges") |>
    mutate(
      target = ifelse(.data$net == "main", target_main, target_casual),
      absolute = .data$edges - .data$target,
      percent = (.data$absolute / .data$target) * 100
    ) |>
    pivot_longer(
      cols = c("absolute", "percent", "edges"),
      names_to = "diff_type",
      values_to = "diff"
    ) |>
    group_by(.data$time, .data$net, .data$diff_type) |>
    mutate(mean = mean(.data$diff, na.rm = TRUE)) |>
    ungroup() |>
    mutate(target = ifelse(.data$diff_type == "edges", .data$target, 0))
}

#' @title Plot Edges History
#' @description Plots the edges history for a specified network and type of difference (absolute, percent, or edges).
#' @param edges_df A data frame containing the edges history,
#' typically obtained from `get_edges_history()`.
#' @param network A character string specifying the network type, either "main" or "casual".
#' @param type A character string specifying the type of difference to plot, either "percent", "absolute", or "edges".
#' @return A ggplot object showing the edges history over time for the specified network and type.
#' @importFrom ggplot2 ggplot aes geom_line geom_hline labs
#' @importFrom dplyr filter pull
#' @importFrom rlang .data
#' @export
plot_edges_history <- function(x, network, type) {
  if (!class(x) %in% c("netsim", "data.frame")) {
    stop("x must be a netsim object or a data frame.")
  }
  if (!type %in% c("percent", "absolute", "edges")) {
    stop("type must be one of 'percent', 'absolute', or 'edges'.")
  }
  if (!network %in% c("main", "casual")) {
    stop("network must be either 'main', or 'casual'.")
  }

  if (class(x) == "netsim") {
    edges_df <- get_edges_history(x, nets = network)
  } else {
    edges_df <- x
  }

  if (!all(c("time", "sim", "net", "target", "diff_type", "diff", "mean") %in% names(edges_df))) {
    stop("edges_df must contain the columns: time, sim, net, target, diff_type, diff, and mean.")
  }

  target_val <- edges_df |>
    filter(net == network, diff_type == type) |>
    pull(target) |>
    unique()

  edges_df |>
    filter(.data$net == network, .data$diff_type == type) |>
    ggplot(aes(x = .data$time, y = .data$diff, color = .data$sim)) +
    geom_line() +
    geom_line(aes(y = .data$mean), color = "black", linewidth = 1) +
    geom_hline(aes(yintercept = target_val)) +
    labs(
      title = paste("Edges history for ", network, " network (", type, ")", sep = ""),
      y = paste(type, "difference"),
      x = "time"
    )
}

#' @title Summarize Final Degrees from Simulation
#' @description Summarizes the final degrees of individuals in the main and casual networks
#' at the end of the simulation and calculates the mean degree for each age and race combination.
#' @param sim A simulation object of class `EpiModel::netsim`.
#' @return A data frame summarizing the mean degree, interquartile range (IQR), and data source
#' for each age and race combination
#' @export

# frequency of rels by age in networks at end of simulation
summarize_final_degrees <- function(sim) {
  simdat <- NULL
  nsims <- sim$control$nsims

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

#' @title Plot Final Degrees for Main and Casual Networks
#' @description Plots the final degrees of individuals in the main and casual networks summarized across simulations
#' and compares them to target degrees extracted from a YAML file.
#' @param sim A simulation object of class `EpiModel::netsim`.
#' @param network A character string specifying the network type, either "main" or "casual".
#' @param yaml_params_loc Path to the YAML file containing target degrees.
#' @return A ggplot object showing the final degrees for the specified network type,
#' comparing simulated degrees to target degrees.
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar facet_wrap
#' @importFrom dplyr filter mutate
#' @importFrom rlang .data
#' @export
plot_final_degrees <- function(sim, network, yaml_params_loc) {
  if (!network %in% c("main", "casual")) {
    stop("network must be either 'main' or 'casual'.")
  }

  s <- summarize_final_degrees(sim)
  t <- get_target_degrees_age_race(yaml_params_loc) |>
    dplyr::mutate(IQR1 = degree, IQR3 = degree) # targets do not have IQRs

  y <- rbind(s, t)

  y |>
    dplyr::filter(.data$type == network) |>
    ggplot2::ggplot(ggplot2::aes(x = .data$age, y = .data$degree, color = .data$data)) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$IQR1, ymax = .data$IQR3), width = 0.2) +
    ggplot2::facet_wrap(~ .data$race)
}


#' @title Get Mean Durations of Relationships at End of Simulation
#' @description Calculates the mean durations of relationships in the main and casual networks at the end
#' of the simulation, comparing them to target durations specified in a YAML file.
#' @param sim A simulation object of class `EpiModel::netsim`.
#' @param nets A character vector specifying the networks to calculate durations for, default is c("main", "casual").
#' @param yaml_params_loc Path to the YAML file containing target durations.
#' @return A data frame summarizing the target and simulated mean durations for each network,
#' along with the standard deviation of the simulated durations.
#' @importFrom yaml read_yaml
#' @export
get_mean_durations <- function(sim, nets = c("main", "casual"), yaml_params_loc) {
  x <- read_yaml(yaml_params_loc)
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
    network = c("main", "casual"),
    target = c(
      x[[nets[1]]]$duration$overall,
      x[[nets[2]]]$duration$overall
    ),
    sim_mean = c(mean(main_durs), mean(casual_durs)),
    sim_sd = c(sd(main_durs), sd(casual_durs))
  )
}
