testthat::test_that(
  "estimate_ML_ordinal returns expected result",
  {

    load(file = testthat::test_path("testdata", "ordinal_test_data.rda"))
    estimate_ML_ordinal_expected <- readRDS(testthat::test_path("testdata", "estimate_ML_ordinal_expected.rds"))
    init_ordinal <- readRDS(testthat::test_path("testdata", "init_ordinal.rds"))

    estimate_ML_ordinal_expected <- setFreqs(estimate_ML_ordinal_expected)

    new_version <- estimate_ML_ordinal(ordinal_test_data, tol = 1e-7, max_iter = 1000, init = init_ordinal)

    unique_ordinal <- unique_obs_summary(ordinal_test_data)
    new_version_fast <- estimate_ML_ordinal(unique_ordinal$unique_obs, unique_ordinal$obs_freq, tol = 1e-7, max_iter = 1000, init = init_ordinal)

    testthat::expect_equal(
      getResults(new_version),
      getResults(estimate_ML_ordinal_expected)
    )

    testthat::expect_equal(
      getResults(new_version_fast)[1:4],
      getResults(estimate_ML_ordinal_expected)[1:4]
    )


  }
)

