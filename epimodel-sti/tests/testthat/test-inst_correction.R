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
