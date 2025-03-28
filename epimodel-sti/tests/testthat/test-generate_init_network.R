# Tests for generate_init_network()
test_that("Initial nw generation sets attributes correctly", {
  # import parameters
  x <- yaml::read_yaml(test_path("input", "nw_params_for_test.yaml"))

  # Test functionality with good params (testing function options & object creation)
  ## expect warning when seed not specifed
  expect_warning(generate_init_network(x))

  ## test that nw gets created properly when deg_casual = FALSE (default)
  expect_no_warning(nw <- generate_init_network(x, seed = 123))
  expect_s3_class(nw, "network")

  ## test that none of the nodal attrs are NA
  attrs1 <- network::list.vertex.attributes(nw)
  for (i in seq_along(attrs1)) {
    vec <- nw %v% i
    expect_equal(sum(is.na(vec)), 0)
  }

  ## expect warning that deg_casual not getting set as attribute if no casual params in yaml
  expect_warning(nw_no_cas <- generate_init_network(x, seed = 123, assign_deg_casual = TRUE))

  ## test that nw gets created and doesn't contain deg_casual attribute
  expect_s3_class(nw_no_cas, "network")
  attrs <- network::list.vertex.attributes(nw_no_cas)
  expect_false("deg_casual" %in% attrs)

  # Test bad parameter inputs
  x$pop$size <- "string"
  expect_error(generate_init_network(x, seed = 123))

  ## error if sex attribute dist != 1
  x$pop$size <- 100
  x$pop$female$dist <- c(0.6, 0.5)
  expect_error(generate_init_network(x, seed = 123))

  ## error if sum of any of the categorical dists != 1 or don't match n categories
  ## (EpiModel apportion_lr error)
  x$pop$female$dist <- c(0.5, 0.5)
  x$pop$race$dist <- rep(0.1, 4)
  expect_error(generate_init_network(x, seed = 123))

  x$pop$race$dist <- c(rep(0.1, 4), 0.6)
  expect_error(generate_init_network(x, seed = 123))
})

test_that("Attribute extraction only works for network objects", {
  # Generate network
  x <- yaml::read_yaml(test_path("input", "nw_params_for_test.yaml"))
  net <- generate_init_network(x, seed = 123)

  # no errors / warnings expect here
  expect_no_warning(attrs <- get_nw_attr_vecs(net))
  expect_no_error(get_nw_attr_vecs(net))

  # check that list contains all network attributes
  attr_names <- network::list.vertex.attributes(net)
  expect_equal(names(attrs), attr_names)

  # check that error occurs when using wrong input
  expect_error(get_nw_attr_vecs(attrs))
})
