# functions used to reshape alters and egos list as prep for egodata object
# these were written by JKB

# nolint start

###########################################################
# reshape_edgelist
###########################################################
# TO DO: Add pre-post reshape checks to the data

#' Reshape an ego file into an ego-alter edgelist
#' @param dat Ego data frame
#' @param n Number of alter variables per ego
#' @param delete_empty Vector of required variables. If ANY are NA, the row will be deleted.
#' @param all Set to TRUE to require that all vars be NA in order to delete; if FALSE,
#' row will be deleted if any of the vars are NA
#' @examples
#' data(nsfg1)
#' nrow(reshape_edgelist(nsfg1, delete_empty = c("active", "rel", "len")))
#' nrow(reshape_edgelist(nsfg1, delete_empty = c("active", "rel", "len"), all = FALSE))
#' @export
reshape_edgelist <- function(dat, n = 3, delete_empty = NULL, all = FALSE) {
  # Variables
  vars <- colnames(dat)
  # Create ID if you need to
  if (!"ego" %in% vars) {
    warning("In reshape_edgelist, creating ego ID variable using row number")
    dat$ego <- 1:nrow(dat)
  }
  # Identify varying vars as all that have numeric suffixes
  # This will probably not be generalizable - will need to get nuanced?
  varythese <- vars[grepl("[0-9]$", vars)]
  varynames <- varnx(varythese, reverse = TRUE, allx = TRUE)
  varylist <- lapply(varynames, function(x) varnx(x))

  dfl <- reshape(dat,
    idvar = c("ego"),
    varying = varylist,
    v.names = varynames,
    times = 1:n,
    timevar = "alter",
    direction = "long",
    sep = ""
  )

  if (!is.null(delete_empty)) {
    if (1 == 0) { # DELETE THIS CODE AFTER I'VE CONFIRMED THAT I DON'T NEED IT
      # Indicator of NA for the cols of interest
      nacols <- apply(dfl[, delete_empty], 2, is.na)
      # Delete if ANY are missing
      delrows <- rowSums(nacols) != 0
      # Deprecated: code to delete if ALL rows are empty

      # Delete rows that are NA for ALL of the variables
      # specified in delete_empty
      delrows <- rowSums(nacols) == ncol(nacols)
    }
    dfl <- delete_empty_rows(dfl, delete_empty, all)
  }
  return(dfl)
}

###########################################################
# var_nx
###########################################################

#' Return variable names with suffixes n to x
#' @param var Vector of variable names
#' @param n Start number for suffixes
#' @param x End number for suffixes
#' @param reverse Logical: set to TRUE to return unique prefixes
#' from a vector of variables with suffixes
#' @param allx Logical: set to TRUE to return only those prefixes in
#' var that have all x variables (presumes counting starts at n=1)
#' @examples
#' data(nsfg1)
#' varnx("len")
#' varnx(c("len1", "len1"), x = 2, reverse = TRUE, allx = TRUE)
#' varnx(c("len1", "len2"), reverse = TRUE, allx = FALSE)
#' @export
varnx <- function(var, n = 1, x = 3, reverse = FALSE, allx = FALSE) {
  if (!reverse) {
    c(sapply(var, function(v) paste0(v, n:x)))
  } else {
    # Strip suffixes
    allprefix <- gsub("[0-9]$", "", var)
    uniques <- unique(allprefix)
    if (!allx) {
      return(uniques)
    } else {
      counts <- sapply(uniques, function(p) sum(allprefix == p))
      return(uniques[counts == x])
    }
  }
}

###########################################################
# rdelete_because_NA
###########################################################

#' Delete rows in a data frame based on NAs in specified variables
#'
#' @param data Data frame
#' @param vars Character vector of variable names in the data to evaluate for NA
#' @param all Set to TRUE to require that all vars be NA in order to delete; if FALSE,
#' row will be deleted if any of the vars are NA
#' @param return The original data frame minus the deleted rows
#' @examples
#' # Load data and use varnx to define variables 1 to 3 of active, rel, len
#' data(nsfg1)
#' nrow(nsfg1)
#' thesevars <- c(varnx("active"), varnx("rel"), varnx("len"))
#' # Compare nrows after deletions
#' nrow(delete_empty_rows(nsfg1, thesevars, all = TRUE))
#' @export
delete_empty_rows <- function(data, vars, all) {
  if (!is.logical(all)) stop("all parameter is not logical")
  # Indicator of NA for the cols of interest
  nacols <- apply(data[, vars], 2, is.na)
  # Deprecated: code to delete if ALL rows are empty
  if (all) {
    # Delete rows that are NA for ALL of the variables
    # specified in delete_empty
    delrows <- rowSums(nacols) == ncol(nacols)
  } else {
    # Delete if ANY are missing
    delrows <- rowSums(nacols) != 0
  }
  data <- data[!delrows, ]
  return(data)
}

# nolint end
