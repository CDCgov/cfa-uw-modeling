#' @title Generate Starting Network for ERGM Estimation
#'
#' @description Generates empty network of specified size and adds nodal attributes.
#' Currently these attributes are specific to this project:
#' a sex variable, race, age group, age, and an adjusted age variable used to capture
#' reported age asymmetry in heterosexual partnerships in NSFG
#' (females on average slightly younger than their male partners)
#'
#' @param params a list of parameters with a sublist called "pop" with the size of the network
#' and the named attributes and their distribution in the population
#' @param seed optional, used to set a seed to maintain the same distribution of nodal attributes
#' every time function is called
#' @param assign_deg_casual default = FALSE, if TRUE, uses additional parameters information in "params" to set
#' a plausible casual degree for each node (to help fit main/long-term partnership network,
#' which is usually fit first in workflow)
#'
#' @return an object of class network
#'
#' @export
#'

generate_init_network <- function(params, seed = NULL, assign_deg_casual = FALSE) {
  if (is.null(seed)) {
    warning("No seed specified")
  } else {
    set.seed(seed)
  }

  # desired population size
  if (!is.numeric(params$pop$size)) {
    stop("Population size must be numeric")
  }
  num <- as.integer(params$pop$size)

  # Generate empty network of desired population size
  nw <- network::network.initialize(num, directed = FALSE)

  ## Create initial attribute vectors
  ### sex variable
  # ensure that dist sums to 1
  if (sum(params$pop$female$dist) != 1) {
    stop("Distrution of sex attribute must sum to 1")
  }
  female <- rep(params$pop$female$levels, params$pop$female$dist * num)

  ### race variable
  race <- sample(EpiModel::apportion_lr(
    num,
    params$pop$race$levels,
    params$pop$race$dist
  ))
  ### age variables
  age_group <- sample(EpiModel::apportion_lr(
    num,
    params$pop$age_group$levels,
    params$pop$age_group$dist
  ))

  # assign age based on age group
  age <- rep(NA, num)
  # calc age group width
  width <- ((params$pop$age$max + 1) - params$pop$age$min) / length(params$pop$age_group$levels)
  # possible ages
  ages <- seq(from = params$pop$age$min, to = params$pop$age$max + 1, by = 1 / 52)

  # assign age
  for (i in seq_along(params$pop$age_group$levels)) {
    low <- params$pop$age$min + width * (i - 1)
    high <- params$pop$age$min + width * i
    ageopts <- ages[which(ages >= low & ages < high)]
    indices <- which(age_group == i)
    age[indices] <- sample(ageopts, length(indices), replace = TRUE)
  }

  # make adjusted age for female asymmetry in partnerships
  age_adj <- age
  age_adj[which(female == 1)] <- age[which(female == 1)] + params$pop$age$female_age_adj
  age_adj[which(age_adj > 50)] <- 50

  if (isTRUE(assign_deg_casual) && !is.null(params$casual)) {
    int_age_range <- params$pop$age$min:params$pop$age$max
    ## set obeserved casual degree from casual params by age & race
    ref <- data.frame(
      age = rep(int_age_range, each = length(params$pop$race$levels)),
      race = rep(params$pop$race$levels, times = length(int_age_range)),
      prob = params$casual$nodefactor$age_race
    )

    ref$comb <- paste0(ref$age, ref$race)

    popvec <- data.frame(comb = paste0(floor(age), race)) |>
      dplyr::left_join(ref, by = "comb")

    deg_casual <- rpois(num, popvec$prob)

    ## make attr lists
    attr_names <- c("female", "race", "age_group", "deg_casual", "age", "age_adj", "agesq")
    attr_values <- list(female, race, age_group, deg_casual, age, age_adj, age^2)
  }

  if (isTRUE(assign_deg_casual) && is.null(params$casual)) {
    warning("assign_deg_casual = TRUE, but there are no casual parameters in yaml.
      Not setting deg_casual in network attributes.")

    attr_names <- c("female", "race", "age_group", "age", "age_adj", "agesq")
    attr_values <- list(female, race, age_group, age, age_adj, age^2)
  }
  if (!assign_deg_casual) {
    attr_names <- c("female", "race", "age_group", "age", "age_adj", "agesq")
    attr_values <- list(female, race, age_group, age, age_adj, age^2)
  }

  # Set attributes on network
  nw <- network::set.vertex.attribute(nw, attr_names, attr_values)

  return(nw)
}

#' @title Calculate Network Target Statistics
#'
#' @description Calculates network targets given, at minimum, a starting network,
#' a list of parameters calculated from data, plus desired ERGM term and attribute.
#' Function currently supports calculations for edges, nodefactor, nodematch, absdiff,
#' and concurrent ERGM terms. Edges, nodefactor, and nodematch terms additionally require
#' a joint distribution (default = age and race) of parameters for calculations,
#' although these can be summed to output the marginal targets (e.g. nodefactor(~age))
#'
#' @param nw the network object outputted from "generate_init_network"
#' (usually already in environment during workflow)
#' @param params the parameter list object with parameters information for
#' all networks estimated from data
#' @param rel string, which relationship/network sub-lists to pull parameters from.
#' @param count_type string, which ERGM term to generate targets for. Currently can only be
#' in the form "edges", "nodefactor", "nodematch", "absdiff_sqrt_age", or "concurrent".
#' @param attr_name string of length 1, used for nodefactor and nodematch targets. If NULL, produces
#' targets based on the two joint_attrs specifed
#' @param joint_attrs string of length 2, the two attrs used to calculate joint distribution
#' of estimates calculated from empirical data
#' @param inst_correct default = FALSE, if TRUE, adjust year-long cumulative reporting of
#' one-time partnerships to daily or weekly counts
#' @param time_unit default = "weeks", the desired time unit for inst reporting conversion
#' @param level additional statification for nodedov target calculation function
#' @param attr_squared for nodecov target calculation, use squared version of attribute? (usually, age)
#' @export
#'

calc_targets <- function(nw, params, rel, count_type,
                         attr_name = NULL, joint_attrs = c("age", "race"),
                         inst_correct = FALSE, level = NULL, attr_squared = FALSE, time_unit = "weeks") {
  # Test that conditions are met to calculate specified target -----------------
  check_conditions(nw, params, rel, count_type, attr_name, joint_attrs)

  # Calculate Targets -------------------------------
  ## (currently assumes we use joint distribution to make edges calculations)
  ## get pop size, attribute vectors
  num <- network::network.size(nw)
  attrs <- EpiModelSTI::get_nw_attr_vecs(nw)
  joint_name <- paste0(joint_attrs[1], "_", joint_attrs[2])

  # calculate expected nodefactor targets
  nf_joint_counts <- calc_joint_nodefactor(params, attrs, joint_attrs, joint_name, rel)

  # calculate expected total edges based on sum of nodefactor
  edges <- calc_edges(nf_joint_counts)

  if (count_type == "edges") {
    final_targets <- edges
  }

  if (count_type == "nodecov") {
    final_targets <- calc_nodecov_age(params, rel, attr_name, edges, level, joint_attrs, nf_joint_counts, attr_squared)
  }

  # if true, calc target for absdiff sqrt age
  if (count_type == "absdiff_sqrt_age") {
    final_targets <- calc_absdiff(params, rel, count_type, edges)
  }

  if (count_type == "concurrent") {
    # calc number of people in rels
    # number of people with > 1 partners
    final_targets <- calc_concurrent(params, rel, num)
  }

  if (count_type == "cross_network") {
    # calc number of people in rels
    # number of people with > 1 partners
    final_targets <- calc_cross_network(params, rel, num)
  }

  # if true, calc nodefactor and then if true nodematch
  if (count_type %in% c("nodefactor", "nodematch")) {
    if (is.null(attr_name)) {
      # targets for full joint distribution
      attr_targets <- as.numeric(nf_joint_counts)
    } else {
      attr_targets <- calc_single_attr_nodefactor(params, attr_name, joint_attrs, nf_joint_counts)
    }

    # if nodefactor, leave targets as-is
    if (count_type == "nodefactor") {
      final_targets <- attr_targets
    }
    # if nodematch, use nodefactor targets using nodefactor info
    if (count_type == "nodematch") {
      final_targets <- calc_nodematch(params, attr_name, attr_targets, rel)
    }
  }

  # Check that these are reasonable targets
  check_targets(edges, final_targets, count_type)

  # Instantaneous Rel Correction
  ## correct for survey data reflecting number of one-times in last year
  if (rel == "inst" && isTRUE(inst_correct)) {
    final_targets <- inst_correction(final_targets, time_unit)
  }

  # Output targets
  return(final_targets)
}

#' @title Extract Nodal Attribute Vectors from Network
#'
#' @description Outputs nodal attribute vectors as list from network object
#'
#' @param nw a network object, in this workflow, usually the network generated from "generate_init_network"
#'
#' @return A list of nodal attribute vectors
#' @export
#'

get_nw_attr_vecs <- function(nw) {
  if (!"network" %in% class(nw)) {
    stop("input must be a network object")
  }

  n <- network::list.vertex.attributes(nw)

  attrs <- list()

  for (i in n) {
    if (i == "age") {
      attrs[[i]] <- floor(nw %v% i)
    } else {
      attrs[[i]] <- nw %v% i
    }
  }

  return(attrs)
}

#' @title Make Empirical Mixing Matrix Symmetrical
#'
#' @description Convert empirical mixing matrix to symmetrical matrix based on mean of upper/lower sections
#'
#' @param mat empirical mixing matrix
#'
#' @return matrix
#' @export

matrix_symmetrical <- function(mat) {
  if (dim(mat)[1] != dim(mat)[2]) {
    stop("Matrix must be square.")
  }

  ncats <- dim(mat)[1]

  newmat <- matrix(NA, ncats, ncats)
  for (i in 1:ncats) {
    for (j in 1:ncats) {
      newmat[i, j] <- mean(c(mat[i, j], mat[j, i]))
    }
  }
  return(newmat)
}

#' @title Target correction for instantaneous network
#'
#' @description In NSFG (and other surveys), information about the frequency of one-time (instantaneous)
#' partnerships are reported as the cumulative number over the last 12 months. This function converts the
#' empirical yearly target to per-week or per-day counts.
#'
#' @param targets a numeric vector of the ergm target statistics estimated from cumulative year data
#' @param time_unit either "days" or "weeks", the unit of time represented by each discrete step in simulation model
#'
#' @return A numeric vector
#' @export

inst_correction <- function(targets, time_unit = NULL) {
  # if available time unit not specifed, return targets unmodified
  # assumes default target time_unit is a year
  if (is.null(time_unit) || !time_unit %in% c("weeks", "days")) {
    warning("Specified time_unit not available, returning umodified targets.")
    unit_correction <- 1
  } else {
    if (time_unit == "weeks") {
      unit_correction <- 52
    } else {
      if (time_unit == "days") {
        unit_correction <- 365
      }
    }
  }

  return(targets / unit_correction)
}

#' @title Target Stats Calculation Helpers
#'
#' @description Small helper functions used as part of calc_targets()
#'
#' @inheritParams calc_targets
#' @param joint_name Name of joint attribute as found in parameter input, calculated in calc_targets()
#'
#'
#'
#' @name targets
NULL


#' @rdname targets
#' @param attrs which two attributes?
#' @export
calc_joint_nodefactor <- function(params, attrs, joint_attrs, joint_name, rel) {
  # pop counts by joint attr dist
  counts <- xtabs(~ attrs[[joint_attrs[1]]] + attrs[[joint_attrs[2]]])

  # shape joint nodefactor input probs into same shape as joint pop counts
  nf_joint_probs <- matrix(params[[rel]][["nodefactor"]][[joint_name]],
    nrow = nrow(counts), byrow = TRUE
  )

  # calculate expected nodefactor targets
  nf_joint_counts <- counts * nf_joint_probs

  return(nf_joint_counts)
}

#' @rdname targets
#' @param nf_joint_counts output from calc_joint_nodefactor
#' @export
calc_single_attr_nodefactor <- function(params, attr_name, joint_attrs, nf_joint_counts) {
  if (attr_name == joint_attrs[1]) {
    attr_targets <- rowSums(nf_joint_counts)
  }
  if (attr_name == joint_attrs[2]) {
    attr_targets <- colSums(nf_joint_counts)
  }

  # condense targets from age into age_group if requested
  if (attr_name == "age_group" && joint_attrs[1] == "age") {
    attr_targets <- rowSums(nf_joint_counts)
    age_range <- ((params$pop$age$max + 1) - params$pop$age$min)
    grp_width <- age_range / length(params$pop$age_group$levels)
    ngrps <- age_range / grp_width
    dat <- data.frame(
      grps = rep(1:ngrps, each = grp_width),
      tar = attr_targets
    )
    summed_targets <- dat |>
      dplyr::group_by(grps) |>
      dplyr::summarize(tars = sum(tar)) |>
      as.data.frame()
    attr_targets <- as.numeric(summed_targets[, "tars"])
  }

  return(attr_targets)
}

calc_nodematch <- function(params, attr_name, attr_targets, rel) {
  if (!attr_name %in% names(params[[rel]][["nodematch"]])) {
    stop("Attr name not available for nodematch statistic.")
  }
  attr_probs_nodematch <- params[[rel]][["nodematch"]][[attr_name]]
  # further refine nodefactor targets to get nodematch
  # how many edges of the estimated above activity are matching
  final_targets <- attr_probs_nodematch * attr_targets / 2

  return(final_targets)
}

#' @rdname targets
#' @param nf_joint_counts output from calc_joint_nodefactor()
#' @export
calc_edges <- function(nf_joint_counts) {
  return(sum(nf_joint_counts) / 2)
}

#' @rdname targets
#' @param edges output from calc_edges()
#' @export
calc_absdiff <- function(params, rel, count_type, edges) {
  avg <- params[[rel]][[count_type]]
  return(avg * edges)
}

#' @rdname targets
#' @export
calc_concurrent <- function(params, rel, num) {
  return(num * params[[rel]][["concurrent"]])
}

#' @rdname targets
#' @export
calc_cross_network <- function(params, rel) {
  return(params$pop$size * params[[rel]][["cross_network"]])
}

#' @rdname targets
#' @param nf_joint_counts output from calc_joint_nodefactor()
#' @param edges output from calc_edges()
#' @export
calc_nodecov_age <- function(params, rel, attr_name, edges, level = NULL,
                             joint_attrs, nf_joint_counts, attr_squared) {
  # first check if cutoff exists
  cutoff <- params[[rel]][["nodecov"]][["cutoff"]]

  if (is.null(cutoff)) {
    # get mean of value for attr(i) + attr(j) from data
    if (attr_squared) {
      attr_name <- paste0(attr_name, "sq")
    }
    attr_mean <- params[[rel]][["nodecov"]][[attr_name]]
  }

  if (!is.null(cutoff)) {
    # Need to re-calculate edges for that cutoff group
    nf_counts <- calc_single_attr_nodefactor(params, attr_name = "age", joint_attrs, nf_joint_counts)
    if (level == "low") {
      age_range <- params$pop$age$min:(cutoff - 1)
    }
    if (level == "high") {
      age_range <- cutoff:params$pop$age$max
    }
    counts <- nf_counts[which(names(nf_counts) %in% age_range)]
    edges <- sum(counts) / 2

    # if cutoff exists, need level to clarify nodecov target
    # get mean of value for attr(i) + attr(j) from data
    if (attr_squared) {
      attr_name <- paste0(attr_name, "sq")
    }
    attr_mean <- params[[rel]][["nodecov"]][[level]][[attr_name]]
  }

  # nodecov target is the attr_mean multiplied by the number of edges
  return(attr_mean * edges)
}

#' @rdname targets
#' @param edges output from calc_edges()
#' @param final_targets vector, final ergm term targets to be checked before output
#' @param threshold default = 0.01, proportion of expected activity based on edges that
#' calulated target is allowed within (+/- threshold)
#' @export
check_targets <- function(edges, final_targets, count_type, threshold = 0.01) {
  # check that all targets are positive
  if (sum(final_targets < 0) > 0) {
    stop("All targets must be >= 0")
  }

  expected_activity <- edges * 2

  if (count_type == "nodefactor") {
    high_threshold <- expected_activity * (1 + threshold)
    low_threshold <- expected_activity * (1 - threshold)
    if ((sum(final_targets) > high_threshold) || (sum(final_targets) < low_threshold)) {
      stop("Sum of nodefactor targets do not match expected activity,
        check distribiton of attribute and activity levels of attribute")
    }
  }

  if (count_type == "nodematch") {
    if (sum(final_targets) > expected_activity) {
      stop("Sum of nodematch targets cannot exceed expected activity level")
    }
  }

  if (count_type %in% c("concurrent", "absdiff_sqrt_age", "nodecov", "cross_network")) {
    if (length(final_targets) > 1) {
      stop("target must be of length 1 for nodecov, concurrent, cross_network and absdiff_sqrt_age, targets")
    }
  }
}

#' @rdname targets
#' @export
# nolint start
check_conditions <- function(nw, params, rel, count_type, attr_name, joint_attrs) {
  if (!"network" %in% class(nw)) {
    stop("inputted initial network must be a network object")
  }
  if (!"list" %in% class(params)) {
    stop("Parameters must be in the form of a list")
  }
  if (!rel %in% names(params)) {
    stop("Specified relationship type does not appear in parameters list")
  }
  if (!count_type %in% c(
    "edges", "nodefactor", "nodematch", "absdiff_sqrt_age",
    "concurrent", "nodecov", "cross_network"
  )) {
    stop("This function is only designed to estimate targets
                for edges, nodefactor, nodematch, absdiff by sqrt age, nodecov, cross_network or concurrent ergm terms")
  }
  if (count_type != "edges" && !count_type %in% names(params[[rel]])) {
    stop("Specified count type does not appear in parameter list for this relationship type")
  }
  if (is.null(joint_attrs)) {
    stop("calculating edges without joint attr distribution not currently supported")
  }
  if (!is.null(attr_name) && !attr_name %in% c(names(params[[rel]][[count_type]]), joint_attrs)) {
    if ((attr_name == "age_group" && "age" %in% joint_attrs) && count_type == "nodefactor") {

    } else {
      stop("Specified attribute name does not appear in parameter list for this relationship and count type")
    }
  }

  joint_name <- paste0(joint_attrs[1], "_", joint_attrs[2])

  if (
    count_type %in% c("nodefactor", "nodematch") &&
      !joint_name %in% names(params[[rel]][["nodefactor"]])
  ) {
    stop("Joint attributes either do not appear in nodefactor list for this relationship
    or are specifed in the wrong order")
  }
}
# nolint end
