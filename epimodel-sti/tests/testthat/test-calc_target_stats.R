test_that("Calculate target stats", {
  # load data & generate init network
  params <- yaml::read_yaml(test_path("input", "full_nw_params_for_test.yaml"))
  nw <- EpiModelSTI::generate_init_network(params, seed = 123)

  # Test that we get no errors in this simple case
  expect_no_error(calc_targets(nw = nw, params = params, rel = "main", count_type = "edges"))

  # Tests for names that don't exist in parameters / wrong objects
  expect_error(calc_targets(nw = data.frame())) ## net should be a network object
  expect_error(calc_targets(nw = nw, params = data.frame())) # #x should be a list
  expect_error(calc_targets(nw = nw, params = params, rel = "marriage")) ## needs count_type
  expect_error(calc_targets(nw = nw, params = params, rel = "main", count_type = "test")) # unavailable count_type
  ## assumes no concurrent rels in main net
  expect_error(calc_targets(nw = nw, params = params, rel = "main", count_type = "concurrent"))
  ## unavailable attr names
  expect_error(calc_targets(nw = nw, params = params, rel = "main", count_type = "nodefactor", attr_name = "test"))
  expect_error(calc_targets(
    nw = nw, params = params, rel = "main", count_type = "nodefactor",
    joint_attrs = c("test", "test2")
  ))
  ## nm only avail for race
  expect_error(calc_targets(
    nw = nw, params = params, rel = "main", count_type = "nodematch",
    attr_name = "age", diff = TRUE
  ))

  # currently need joint attrs dist for all count types
  expect_error(calc_targets(nw = nw, params = params, rel = "main", count_type = "edges", joint_attrs = NULL))
  expect_error(calc_targets(nw = nw, params = params, rel = "main", count_type = "nodefactor", joint_attrs = NULL))
  expect_error(calc_targets(nw = nw, params = params, rel = "main", count_type = "nodematch", joint_attrs = NULL))
  expect_error(calc_targets(
    nw = nw, params = params, rel = "main", count_type = "absdiff_sqrt_age",
    joint_attrs = NULL
  ))
  expect_error(calc_targets(nw = nw, params = params, rel = "main", count_type = "concurrent", joint_attrs = NULL))

  # expect correct length of returned results
  rels <- names(params)[-1]
  for (rel in rels) {
    # edges
    edge_tar <- calc_targets(nw = nw, params = params, rel = rel, count_type = "edges")
    expect_equal(length(edge_tar), 1)
    # nodefactor
    nf_tar <- calc_targets(nw = nw, params = params, rel = rel, count_type = "nodefactor")
    nf_data <- params[[rel]][["nodefactor"]][["age_race"]]
    expect_equal(length(nf_tar), length(nf_data))
    # nodematch (currently only available for race in main/casual nets)
    if (!rel %in% "inst") {
      nm_tar <- calc_targets(
        nw = nw, params = params, rel = rel, count_type = "nodematch",
        attr_name = "race", diff = TRUE
      )
      nm_data <- params[[rel]][["nodematch"]][["race"]]
      expect_equal(length(nm_tar), length(nm_data))
    }
  }

  # concurrent (only used for casual)
  expect_equal(length(calc_targets(nw = nw, params = params, rel = "casual", count_type = "concurrent")), 1)
})
