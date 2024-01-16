#' Create data sets which simulate paired measurements of multiple methods
#'
#' @inheritDotParams generate_multimethod_binary
#' @inheritDotParams generate_multimethod_ordinal
#' @inheritDotParams generate_multimethod_continuous
#' @param type A string specifying the data type of the methods under evaluation
#' @param n_method An integer representing the number of methods to simulate
#' @param n_obs An integer representing the number of observations to simulate
#' @param prev Optional value between 0-1 which represents the proportion of positive results in the target population
#' @param D Optional binary vector representing the true classification of each observation
#' @param method_names Optional vector of names used to identify each method
#' @param obs_names Optional vector of names used to identify each observation
#' @param ...
#'
#' @return A list containing a simulated data set and the parameters used to create it
#' @export

generate_multimethod_data <-
  function(
    type = c("binary", "ordinal", "continuous"),
    n_method = 3,
    n_obs = 100,
    prev = 0.5,
    D = NULL,
    method_names = NULL,
    obs_names = NULL,
    ...){
    type <- match.arg(type)
    switch (type,
            binary = generate_multimethod_binary(n_method, n_obs, prev, D, method_names = method_names, obs_names = obs_names, ...),
            ordinal = generate_multimethod_ordinal(n_method, n_obs, prev, D, method_names = method_names, obs_names = obs_names, ...),
            continuous = generate_multimethod_continuous(n_method, n_obs, prev, D, method_names = method_names, obs_names = obs_names, ...)
    )
  }

#' Estimate maximum likelihood accuracy statistics by expectation maximization
#'
#' @inheritDotParams estimate_ML_ordinal
#' @param type A string specifying the data type of the methods under evaluation
#' @param data An n_obs by n_method matrix containing the observed values for each method. If dimensions are named, row names will be used for obs_names, and column names will be used for method_names
#' @param init An optional list of initial values used to seed the EM algorithm. If initial values are not provided, a function will be called on the data to estimate starting values.
#' @param max_iter The maximum number of EM algorithm iterations to compute before reporting a result.
#' @param tol The minimum change in statistic estimates needed to continue iterating the EM algorithm.
#' @param save_progress A logical indication whether to save interim calculations used in the EM algorithm.
#' @param ...
#'
#' @return a MultiMethodMLEstimate class S4 object containing ML accuracy statistics by EM
#' @export

estimate_ML <-
  function(
    type = c("binary", "ordinal", "continuous"),
    data,
    init = list(NULL),
    max_iter = 1000,
    tol = 1e-7,
    save_progress = TRUE,
    ...){
    type <- match.arg(type)
    switch (type,
            binary = estimate_ML_binary(data = data, init = init, max_iter = max_iter, tol = tol, save_progress = save_progress, ... = ...),
            ordinal = estimate_ML_ordinal(data = data, init = init, max_iter = max_iter, tol = tol, save_progress = save_progress, ... = ...),
            continuous = estimate_ML_continuous(data = data, init = init, max_iter = max_iter, tol = tol, save_progress = save_progress, ... = ...)
    )
  }

#' Create plots of ML estimation process
#'
#' @param type A string specifying the data type of the methods under evaluation
#' @param ML_est A MultiMethodEstimate class object
#' @param params An optional list of population parameters
#'
#' @return a list of ggplot2 plots
#' @export

plot_ML <-
  function(
    ML_est,
    params = NULL){
    type <- ML_est@type
    switch (type,
            binary = plot_ML_binary(ML_est, params),
            ordinal = plot_ML_ordinal(ML_est, params),
            continuous = plot_ML_continuous(ML_est, params)
    )
  }


#' Bootstrap Accuracy Statistic Estimation for Multimethod Data
#'
#' @inheritDotParams estimate_ML
#' @param n_boot number of bootstrap estimates to compute
#' @param seed optional seed for RNG
#'
#' @return a list containing accuracy estimates `v` and the parameters used
#' @export
#'
boot_ML <-
  function(
    type = c("binary", "ordinal", "continuous"),
    data,
    n_boot = 100,
    max_iter = 1000,
    tol = 1e-7,
    seed = NULL){

  if(!is.null(seed)) set.seed(seed)

    type <- match.arg(type)

    v_0 <- estimate_ML(type, data, save_progress = FALSE)

    n_obs <- nrow(data)

    v_star <-
    lapply(1:n_boot, function(b){
      tmp <- data[sample(n_obs, n_obs, replace = TRUE), ]
      estimate_ML(type, tmp, save_progress = FALSE)
    })

    return(
      list(
        v_0 = v_0,
        v_star = v_star,
        params = list(
          data = data,
          n_boot = n_boot,
          max_iter = max_iter,
          tol = tol,
          seed = seed
        )

      )
    )

  }


