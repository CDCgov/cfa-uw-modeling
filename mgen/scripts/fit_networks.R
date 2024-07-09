# This script generates an example set of network objects for testing purposes

library(EpiModel)

# Define population size
num <- 10000

# Generate empty network
nw <- network::network.initialize(num, directed = FALSE)

# Set initial attributes
sex <- rep(c(0, 1), each = num / 2)
nw <- network::set.vertex.attribute(nw, "sex", sex)

# Network targets (density, duration in days)
main_tar <- c(2000, 365 * 3)
casual_tar <- c(1000, 365)
once_tar <- c(30, NA)

# Constraints
constraints <-
  ~ blocks(attr = ~sex, levels2 = diag(TRUE, 2)) + bd(maxout = 1) + sparse


# 1. Main Model -----------------------------------------------------------

# Formulas
fit_m <- netest(nw,
  formation = ~edges,
  target.stats = main_tar[1],
  coef.diss = dissolution_coefs(
    dissolution = ~ offset(edges),
    duration = main_tar[2],
    d.rate = 0
  ),
  constraints = constraints
)

# 2. Casual Model ---------------------------------------------------------

# Initialize network
nw_pers <- nw

# Assign main degree from above fitted network
nw_pers %v% "deg.main" <- get_degree(fit_m$newnetwork) # nolint

# Fit model
fit_p <- netest(nw_pers,
  formation = ~edges,
  target.stats = casual_tar[1],
  coef.diss = dissolution_coefs(
    dissolution = ~ offset(edges),
    duration = casual_tar[2],
    d.rate = 0
  ),
  constraints = constraints
)

# Fit inst model ----------------------------------------------------------

# Initialize network
nw_inst <- nw_pers

# Update pers degree
nw_inst %v% "deg.pers" <- get_degree(fit_p$newnetwork) # nolint

# Fit model
fit_i <- netest(nw_inst,
  formation = ~edges,
  target.stats = once_tar[1],
  coef.diss = dissolution_coefs(~ offset(edges), 1),
  constraints = constraints
)

# save out all three networks as list
# this becomes input for epimodel code
est <- list(fit_m, fit_p, fit_i)
save(est, file = paste0("mgen/data/fit.rda"))
