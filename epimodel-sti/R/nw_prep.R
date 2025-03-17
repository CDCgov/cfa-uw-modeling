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
#' @param deg_casual default = FALSE, if TRUE, uses additional parameters information in "params" to set
#' a plausible casual degree for each node (to help fit main/long-term partnership network,
#' which is usually fit first in workflow)
#'
#' @return an object of class network
#'
#' @export
#'

generate_init_network <- function(params, deg_casual = FALSE, seed = NULL) {
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

  if (deg_casual && !is.null(params$casual)) {
    ## set obeserved casual degree from casual params by age group & race
    ref <- data.frame(
      age_group = rep(params$pop$age_group$levels, each = length(params$pop$race$levels)),
      race = rep(params$pop$race$levels, times = length(params$pop$age_group$levels)),
      prob = params$casual$nodefactor$age_group_race
    )

    ref$comb <- paste0(ref$age_group, ref$race)

    popvec <- data.frame(comb = paste0(age_group, race)) |>
      dplyr::left_join(ref, by = "comb")

    deg_casual <- rbinom(num, 1, popvec$prob)

    ## make attr lists
    attr_names <- c("female", "race", "age_group", "deg_casual", "age", "age_adj")
    attr_values <- list(female, race, age_group, deg_casual, age, age_adj)
  }

  if (deg_casual && is.null(params$casual)) {
    warning("deg_casual = TRUE, but there are no casual parameters in yaml.
      Not setting deg_casual in network attributes.")

    attr_names <- c("female", "race", "age_group", "age", "age_adj")
    attr_values <- list(female, race, age_group, age, age_adj)
  }
  if (!deg_casual) {
    attr_names <- c("female", "race", "age_group", "age", "age_adj")
    attr_values <- list(female, race, age_group, age, age_adj)
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
#' @export
#'

# nolint start
calc_targets <- function(nw, params, rel, count_type,
                         attr_name = NULL, joint_attrs = c("age", "race"),
                         inst_correct = FALSE, time_unit = "weeks") {
  # Test that conditions are met to calculate specified target -----------------

  if (!"network" %in% class(nw)) {
    stop("inputted initial network must be a network object")
  }
  if (!"list" %in% class(params)) {
    stop("Parameters must be in the form of a list")
  }
  if (!rel %in% names(params)) {
    stop("Specified relationship type does not appear in parameters list")
  }
  if (!count_type %in% c("edges", "nodefactor", "nodematch", "absdiff_sqrt_age", "concurrent")) {
    stop("This function is only designed to estimate targets
                for edges, nodefactor, nodematch, absdiff by sqrt age, or concurrent ergm terms")
  }
  if (count_type != "edges" && !count_type %in% names(params[[rel]])) {
    stop("Specified count type does not appear in parameter list for this relationship type")
  }
  if (is.null(joint_attrs)) {
    stop("calculating edges without joint attr distribution not currently supported")
  }
  if (!is.null(attr_name) && !attr_name %in% c(names(params[[rel]][[count_type]]), joint_attrs)) {
    stop("Specified attribute name does not appear in parameter list for this relationship and count type")
  }

  joint_name <- paste0(joint_attrs[1], "_", joint_attrs[2])

  if (count_type %in% c("nodefactor", "nodematch") &&
    !joint_name %in% names(params[[rel]][["nodefactor"]])) {
    stop("Joint attributes either do not appear in nodefactor list for this relationship
    or are specifed in the wrong order")
  }

  # Calculate Edges (use in all ERGM target calcs) -------------------------------
  ## (currently assumes we use joint distribution to make this calculations)
  ## get pop size, attribute vectors
  num <- network::network.size(nw)
  attrs <- EpiModelSTI::get_nw_attr_vecs(nw)

  # pop counts by joint attr dist
  counts <- xtabs(~ attrs[[joint_attrs[1]]] + attrs[[joint_attrs[2]]])

  # shape joint nodefactor input probs into same shape as joint pop counts
  nf_joint_probs <- matrix(params[[rel]][["nodefactor"]][[joint_name]],
    nrow = nrow(counts), byrow = TRUE
  )

  # calculate expected nodefactor targets
  nf_joint_counts <- counts * nf_joint_probs

  # calculate expected total edges based on sum of nodefactor
  edges <- sum(nf_joint_counts) / 2

  if (count_type == "edges") {
    final_targets <- edges
  }

  # if true, calc target for absdiff sqrt age
  if (count_type == "absdiff_sqrt_age") {
    avg <- params[[rel]][[count_type]]
    final_targets <- avg * edges
  }

  if (count_type == "concurrent") {
    # calc number of people in rels
    # number of people with > 1 partners
    final_targets <- num * params[[rel]][["concurrent"]]
  }

  # if true, calc nodefactor and then if true nodematch
  if (count_type %in% c("nodefactor", "nodematch")) {
    if (is.null(attr_name)) {
      attr_targets <- as.numeric(nf_joint_counts)
    } else {
      if (attr_name == joint_attrs[1]) {
        attr_targets <- rowSums(nf_joint_counts)
      }
      if (attr_name == joint_attrs[2]) {
        attr_targets <- colSums(nf_joint_counts)
      }

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
    }
    # Check that these are reasonable targets
    expected_activity <- edges * 2

    # test should be slightly different for nodematch
    if ((sum(attr_targets) > expected_activity * 1.01) ||
      (sum(attr_targets) < expected_activity * 0.99)) {
      stop("Sum of target does not match expected activity,
        check distribiton of attribute and activity levels of attribute")
    }

    if (count_type == "nodefactor") {
      final_targets <- attr_targets
    }

    if (count_type == "nodematch") {
      if (!attr_name %in% names(params[[rel]][["nodematch"]])) {
        stop("Attr name not available for nodematch statistic.")
      }
      attr_probs_nodematch <- params[[rel]][["nodematch"]][[attr_name]]
      # further refine nodefactor targets to get nodematch
      # how many edges of the estimated above activity are matching
      final_targets <- attr_probs_nodematch * attr_targets / 2
    }
  }

  # Instantaneous Rel Correction
  ## correct for survey data reflecting number of one-times in last year
  if (rel == "inst" && isTRUE(inst_correct)) {
    if (time_unit == "weeks") {
      unit_correction <- 52
    } else if (time_unit == "days") {
      unit_correction <- 365
    }
    final_targets <- final_targets / unit_correction
  }

  # Output targets
  return(final_targets)
}
# nolint end

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
