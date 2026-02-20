testthat::test_that(
  "estimate_ML_binary returns expected result",
  {

    load(file = testthat::test_path("testdata", "binary_test_data.rda"))
    estimate_ML_binary_expected <- readRDS(testthat::test_path("testdata", "estimate_ML_binary_expected.rds"))

    new_version <- estimate_ML_binary(binary_test_data$generated_data, tol = 1e-7, max_iter = 1000, init = list(prev_1 = 0.2, se_1 = rep(.75, 4), sp_1 = rep(0.75, 4)))

    testthat::expect_equal(
      new_version@results,
      estimate_ML_binary_expected@results

    )
  }
)

