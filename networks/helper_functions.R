# Helper functions for target stats
# nolint start
calc_targets <- function(net = nw, params = x, rel, count_type, attr_name = NULL,
                         joint_attrs = c("age", "race"),
                         dummy_var = FALSE, dummy_var_position = NULL,
                         pop_attr_name = "pop", dist_name = "dist",
                         female_dist = NULL,
                         inst_correct = FALSE, time_unit = "weeks") {
  if (!count_type %in% c("edges", "nodefactor", "nodematch", "absdiff_sqrt_age", "concurrent")) {
    stop("This function is only designed to estimate targets
                for edges, nodefactor, nodematch, absdiff, or concurrent ergm terms")
  }

  if (!is.null(attr_name) && !attr_name %in% c("age", "age_group", "race")) {
    stop("This function is only designed to estimate targets
                  for age, age_group, race attributes")
  }

  if (is.null(joint_attrs)) {
    stop("calculating edges without joint attr distribution not currently supported")
  }

  # get pop size, attribute vectors
  num <- network::network.size(net)
  n <- network::list.vertex.attributes(net)

  attrs <- list()
  for (i in n) {
    if (i == "age") {
      attrs[[i]] <- floor(net %v% i)
    } else {
      attrs[[i]] <- net %v% i
    }
  }

  # calc edges based on joint attr distribution
  # this also gets us joint nodefactor counts

  # pop counts by joint attr dist
  counts <- xtabs(~ attrs[[joint_attrs[1]]] + attrs[[joint_attrs[2]]])
  nf_name <- paste0(joint_attrs[1], "_", joint_attrs[2])

  # shape joint nodefactor input probs into same shape as joint pop counts
  nf_joint_probs <- matrix(params[[rel]][["nodefactor"]][[nf_name]],
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
        age_range <- ((x$pop$age$max + 1) - x$pop$age$min)
        grp_width <- age_range / length(x$pop$age_group$levels)
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


    if (count_type == "nodematch") {
      attr_probs_nodematch <- params[[rel]][["nodematch"]][[attr_name]]
      # further refine nodefactor targets to get nodematch
      # how many edges of the estimated above activity are matching
      attr_targets <- attr_probs_nodematch * attr_targets / 2
    }

    if (dummy_var) {
      if (is.null(dummy_var_position)) {
        final_targets <- attr_targets[-1]
      } else {
        final_targets <- attr_targets[-dummy_var_position]
      }
    } else {
      final_targets <- attr_targets
    }
  }

  # correct for survey data reflecting number of one-times in last year
  if (inst_correct) {
    if (time_unit == "weeks") {
      unit_correction <- 52
    } else if (time_unit == "days") {
      unit_correction <- 365
    }
    final_targets <- final_targets / unit_correction
  }

  # output targets
  return(final_targets)
}

# make matrix symmetrical
# (mean of race-alter_rate and alter-race-race reports)
matrix_symmetrical <- function(mat, ncats) {
  newmat <- matrix(NA, ncats, ncats)
  for (i in 1:ncats) {
    for (j in 1:ncats) {
      newmat[i, j] <- mean(c(mat[i, j], mat[j, i]))
    }
  }
  return(newmat)
}

generate_init_network <- function(x, seed, deg_casual = TRUE) {
  set.seed(seed)
  num <- x$pop$size
  nw <- network::network.initialize(num, directed = FALSE)

  ## Create initial attribute vectors
  ### sex variable
  female <- rep(x$pop$female$levels, x$pop$female$dist * num)

  ### race variable
  race <- sample(EpiModel::apportion_lr(
    num,
    x$pop$race$levels,
    x$pop$race$dist
  ))
  ### age variables
  age_group <- sample(EpiModel::apportion_lr(
    num,
    x$pop$age_group$levels,
    x$pop$age_group$dist
  ))

  # assign age based on age group
  age <- rep(NA, num)
  # calc age group width
  width <- ((x$pop$age$max + 1) - x$pop$age$min) / length(x$pop$age_group$levels)
  # possible ages
  ages <- seq(from = x$pop$age$min, to = x$pop$age$max + 1, by = 1 / 52)

  # assign age
  for (i in seq_along(x$pop$age_group$levels)) {
    low <- x$pop$age$min + width * (i - 1)
    high <- x$pop$age$min + width * i
    ageopts <- ages[which(ages >= low & ages < high)]
    indices <- which(age_group == i)
    age[indices] <- sample(ageopts, length(indices), replace = TRUE)
  }

  # make adjusted age for female asymmetry in partnerships
  age_adj <- age
  age_adj[which(female == 1)] <- age[which(female == 1)] + x$pop$age$female_age_adj
  age_adj[which(age_adj > 50)] <- 50

  if (deg_casual) {
    ## Third, set obeserved casual degree from casual params by age group & race
    ref <- data.frame(
      age_group = rep(x$pop$age_group$levels, each = length(x$pop$race$levels)),
      race = rep(x$pop$race$levels, times = length(x$pop$age_group$levels)),
      prob = x$casual$nodefactor$age_group_race
    )

    ref$comb <- paste0(ref$age_group, ref$race)

    popvec <- data.frame(comb = paste0(age_group, race)) |>
      dplyr::left_join(ref, by = "comb")

    deg_casual <- rbinom(num, 1, popvec$prob)

    ## Fourth, assign attrs to nodes on network
    attr_names <- c("female", "race", "age_group", "deg_casual", "age", "age_adj")
    attr_values <- list(female, race, age_group, deg_casual, age, age_adj)
  } else {
    ## Fourth, assign attrs to nodes on network
    attr_names <- c("female", "race", "age_group", "age", "age_adj")
    attr_values <- list(female, race, age_group, age, age_adj)
  }

  nw <- network::set.vertex.attribute(nw, attr_names, attr_values)
  return(nw)
}
# nolint end
