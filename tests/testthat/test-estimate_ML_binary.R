testthat::test_that(
  "estimate_ML_binary returns expected result",
  {

    load(file = testthat::test_path("testdata", "binary_test_data.rda"))
    estimate_ML_binary_expected <- readRDS(testthat::test_path("testdata", "estimate_ML_binary_expected.rds"))

    estimate_ML_binary_expected <- setFreqs(estimate_ML_binary_expected)

    new_version <- estimate_ML_binary(binary_test_data$generated_data, tol = 1e-7, max_iter = 1000, init = list(prev_1 = 0.2, se_1 = rep(.75, 4), sp_1 = rep(0.75, 4)))

    unique_binary <- unique_obs_summary(binary_test_data$generated_data)

    new_version_fast <- estimate_ML_binary(unique_binary$unique_obs, freqs = unique_binary$obs_freqs, tol = 1e-7, max_iter = 1000, init = list(prev_1 = 0.2, se_1 = rep(.75, 4), sp_1 = rep(0.75, 4)))

    testthat::expect_equal(
      getResults(new_version),
      getResults(estimate_ML_binary_expected)
    )

    testthat::expect_equal(
      getResults(new_version_fast)[1:3],
      getResults(estimate_ML_binary_expected)[1:3]
    )

  }
)

