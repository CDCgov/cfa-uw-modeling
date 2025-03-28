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
  expect_error(calc_targets(nw = nw, params = params, rel = "main", count_type = "nodematch", attr_name = "age"))

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
      nm_tar <- calc_targets(nw = nw, params = params, rel = rel, count_type = "nodematch", attr_name = "race")
      nm_data <- params[[rel]][["nodematch"]][["race"]]
      expect_equal(length(nm_tar), length(nm_data))
    }
  }

  # concurrent (only used for casual)
  expect_equal(length(calc_targets(nw = nw, params = params, rel = "casual", count_type = "concurrent")), 1)
})

# Testing of calc_targets sub-functions
test_that("correction for instantanous rel data collection modifies target correctly", {
  # set targets for evalulation
  cuml_rels <- 500
  per_week_target <- cuml_rels / 52
  per_day_target <- cuml_rels / 365

  cuml_rels_vec <- c(500, 1000, 10)
  per_week_target_vec <- cuml_rels_vec / 52
  per_day_target_vec <- cuml_rels_vec / 365

  # test that function gives correct answer with expected inputs
  expect_equal(inst_correction(cuml_rels, "weeks"), per_week_target)
  expect_equal(inst_correction(cuml_rels, "days"), per_day_target)
  expect_equal(inst_correction(cuml_rels_vec, "weeks"), per_week_target_vec)
  expect_equal(inst_correction(cuml_rels_vec, "days"), per_day_target_vec)

  # expect warnings with unexpected time_unit inputs and returns unmodified target
  expect_warning(t1 <- inst_correction(cuml_rels))
  expect_warning(t2 <- inst_correction(cuml_rels, time_unit = "seconds"))
  expect_warning(t3 <- inst_correction(cuml_rels, time_unit = NA))
  expect_warning(t4 <- inst_correction(cuml_rels, time_unit = 0))
  expect_identical(t1, cuml_rels)
  expect_identical(t2, cuml_rels)
  expect_identical(t3, cuml_rels)
  expect_identical(t4, cuml_rels)

  # expect error when no/wrong targets specified
  expect_error(inst_correction(time_unit = "days"))
  expect_error(inst_correction("test", time_unit = "days"))
  bad_input <- list(c(1, 2, 3), c("test"))
  bad_input2 <- list(c(1, 2, 3), c(5, 6, 7))
  expect_error(inst_correction(bad_input, time_unit = "days"))
  expect_error(inst_correction(bad_input2, time_unit = "days"))
})
