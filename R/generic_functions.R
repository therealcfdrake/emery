#' @title Create data sets which simulate paired measurements of multiple methods
#' @description
#' `generate_multimethod_data()` is a general function for creating a data set which
#' simulates the results one might see when using several different methods to measure a set
#' of objects.
#' @order 1
#' @param type A string specifying the data type of the methods being simulated.
#' @param n_method An integer representing the number of methods to simulate.
#' @param n_obs An integer representing the number of observations to simulate.
#' @param prev A value between 0-1 which represents the proportion of
#' "positive" results in the target population.
#' @param D Optional binary vector representing the true classification of
#'  each observation.
#' @param method_names Optional vector of names used to identify each method.
#' @param obs_names Optional vector of names used to identify each observation.
#' @param ... Additional parameters
#' @returns A list containing a simulated data set and the parameters used to create it
#' @details
#' The function supports binary measurement methods, e.g., Pass/Fail;
#' ordinal measurement methods, e.g., the Likert scale; and continuous measurement
#' methods, e.g., height. The data are generated under the assumption that the
#' underlying population consists of a mixture of two groups. The primary
#' application of this is to simulate a sample from a population which has some
#' prevalence of disease.
#'
#' @export
#' @example man/examples/ML_example.R

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

#' @title Estimate maximum likelihood accuracy statistics by expectation maximization
#' @description
#' `estimate_ML()` is a general function for estimating the maximum likelihood accuracy
#' statistics for a set of methods with no known reference value, i.e. "truth", or
#' "gold standard".
#' @param type A string specifying the data type of the methods under evaluation.
#' @param data An `n_obs` by `n_method` \code{\link{matrix}} containing the
#' observed values for each method. If the dimensions are named, row names will
#' be used to name each observation (`obs_names`) and column names will be used
#' to name each measurement method (`method_names`).
#' @param init An optional list of initial values used to seed the EM algorithm.
#'  If initial values are not provided, the `pollinate_ML()` function will be
#'  called on the data to estimate starting values. It is recommended to try several
#'  sets of starting parameters to ensure that the algorithm converges to the same
#'  results. This is to verify that the result does not represent a local extrema.
#' @param max_iter The maximum number of EM algorithm iterations to compute before
#' reporting a result.
#' @param tol The minimum change in statistic estimates needed to continue
#' iterating the EM algorithm.
#' @param save_progress A logical indication of whether to save interim
#' calculations used in the EM algorithm.
#' @param ... Additional arguments
#' @order 1
#' @returns `estimate_ML()` returns an S4 object of class "MultiMethodMLEstimate"
#' containing the maximum likelihood accuracy statistics calculated by EM.
#' @details
#' The lack of an infallible reference method is referred to
#' as an imperfect gold standard (GS). Accuracy statistics which rely on a GS
#' method, such as sensitivity, specificity, and AUC,
#' can be estimated using imperfect gold standards by iteratively estimating the
#' maximum likelihood values of these statistics while the conditional independence
#' assumption holds. `estimate_ML()` relies on a collection of expectation maximization (EM) algorithms
#' to achieve this. The EM algorithms used in this function are based on those presented in
#' *Statistical Methods in Diagnostic Medicine, Second Edition*
#' \insertCite{Zhou_Obuchowski_McClish_2011}{emery} and have been validated on
#' several examples therein. Additional details about these algorithms can be found
#' for binary \insertCite{Walter1988-oq}{emery}, ordinal \insertCite{Zhou2005-gk}{emery},
#'  and continuous \insertCite{Hsieh_Su_Zhou_2011}{emery} methods.
#'  Minor changes to the literal calculations have been
#' made for efficiency, code readability, and the like, but the underlying steps
#' remain functionally unchanged.
#' @importFrom Rdpack reprompt
#' @export
#' @example man/examples/ML_example.R
#' @references
#' \insertRef{Zhou_Obuchowski_McClish_2011}{emery}
#'
#' \insertRef{Walter1988-oq}{emery}
#'
#' \insertRef{Zhou2005-gk}{emery}
#'
#' \insertRef{Hsieh_Su_Zhou_2011}{emery}
#'

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

#' @title Create plots visualizing the ML estimation process and results.
#' @description
#' `plot_ML()` is a general function for visualizing results generated by `estimate_ML()`.
#' @param ML_est A MultiMethodMLEstimate class object
#' @param params A list of population parameters. This is primarily used to evaluate
#' results from a simulation where the target parameters are known, but can be used to visualize
#' results with respect to some True value.
#' @order 1
#' @returns A list of ggplot2 plots.
#' @export
#' @example man/examples/ML_example.R

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

#' @title Generate seed values for EM algorithm
#' @description
#' `pollinate_ML()` is a general helper function which can be used to generate starting
#' values, i.e. seeds, for the `estimate_ML()` function from a multi-method data set.
#' @inheritParams estimate_ML
#' @param ... Additional arguments
#' @order 1
#' @returns a list of EM algorithm initialization values
#' @export

pollinate_ML <-
  function(
    type = c("binary", "ordinal", "continuous"),
    data,
    ...){
    type <- match.arg(type)
    switch (type,
            binary = pollinate_ML_binary(type, data, ...),
            ordinal = pollinate_ML_ordinal(type, data, ...),
            continuous = pollinate_ML_continuous(type, data, ...)
    )
  }


#' @title Bootstrap ML accuracy statistic estimation for multi-method data
#' @description
#' `boot_ML()` is a function used to generate bootstrap estimates of results generated
#' by `estimate_ML()` primarily for use in creating nonparametric confidence intervals.
#'
#' @inheritParams estimate_ML
#'
#' @param n_boot number of bootstrap estimates to compute
#' @param n_study sample size to select for each bootstrap estimate
#' @param seed optional seed for RNG
#' @returns a list containing accuracy estimates, `v`, and the parameters used.
#' \item{v_0}{result from original data}
#' \item{v_star}{list containing results from each bootstrap resampling}
#' \item{params}{list containing the parameters used}
#' @export
#' @example man/examples/bootstrap_example.R

boot_ML <-
  function(
    type = c("binary", "ordinal", "continuous"),
    data,
    n_boot = 100,
    n_study = NULL,
    max_iter = 1000,
    tol = 1e-7,
    n_obs = NULL,
    seed = NULL,
    ...){

  if(!is.null(seed)) set.seed(seed)

    type <- match.arg(type)

    v_0 <- estimate_ML(type, data, save_progress = FALSE)

    if(is.null(n_obs)) n_obs <- nrow(data)
    if(is.null(n_study)) n_study <- nrow(data)

    v_star <-
    lapply(1:n_boot, function(b){
      tmp <- data[sample(n_obs, n_study, replace = TRUE), ]
      estimate_ML(type, tmp, save_progress = FALSE)@results
    })

    return(
      list(
        v_0 = v_0,
        v_star = v_star,
        params = list(
          data = data,
          n_boot = n_boot,
          n_study = n_study,
          max_iter = max_iter,
          tol = tol,
          n_obs = n_obs,
          seed = seed
        )
      )
    )
  }


#' @description
#' `aggregate_boot_ML()` rearranges the bootstrap results from `boot_ML()` by
#' statistic instead of bootstrap iteration.
#'
#' @param boot_ML_result a list returned by `bootML()`.
#' @returns a named list of long format data frames containing aggregated statistic estimates.
#' \item{boot_id}{index of bootstrap sample which resulted in value}
#' \item{col_id}{value identifier}
#' \item{row_id}{optional value identifier used when the result has more than 1 dimension}
#' \item{value}{statistic value}
#' @export
#' @importFrom stats setNames
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr any_of
#' @example man/examples/bootstrap_example.R

aggregate_boot_ML <-
  function(
    boot_ML_result = NULL
  ){
    lapply(names(boot_ML_result$v_0@results), function(w){
      lapply(1:length(boot_ML_result$v_star), function(l){
        cbind("boot_id" = l, do.call(cbind, boot_ML_result$v_star[[l]][w]))
      }) |>
        do.call(what = rbind) |>
        {\(.) if(!is.null(rownames(.))) cbind(data.frame(., row.names = NULL), row_id = rownames(.)) else data.frame(.)}() |>
        tidyr::pivot_longer(-tidyr::any_of(c("boot_id", "row_id")), names_to = "col_id")
    }) |> stats::setNames(names(boot_ML_result$v_0@results))
  }

#' @description
#' `plot_boot_ML()` creates univariate plots of bootstrap results from `boot_ML()`.
#' @inheritParams estimate_ML
#' @returns a named list of named plots.
#' @export
#' @import ggplot2

 plot_boot_ML <-
   function(boot_ML_result, probs = c(0.05, 0.50, 0.95)){
     agg_results <- aggregate_boot_ML(boot_ML_result)
     stats_to_plot <- names(agg_results)[names(agg_results) %in% c("prev_est", "se_est", "sp_est", "A_i_est", "A_j_est")]

     lapply(stats_to_plot, function(x){

     q_summary <-
       dplyr::reframe(agg_results[[x]],
                      value = stats::quantile(value, probs, na.rm = TRUE),
                      quantile = probs,
                      .by = any_of(c("col_id", "row_id"))) #|>
       # dplyr::summarize(upper = max(value), lower = min(value), .by = any_of(c("col_id", "row_id")))

     # print(q_summary)

       ggplot2::ggplot(agg_results[[x]], aes(x = value, color = if("row_id" %in% colnames(agg_results[[x]])) row_id else "Group")) +
         ggplot2::geom_histogram(bins = 100, boundary = 0, position = "identity", aes(fill = after_scale(alpha(colour, 0.5)))) +
         ggplot2::geom_vline(data = q_summary, lty = 2, aes(xintercept = value, color = if("row_id" %in% colnames(agg_results[[x]])) row_id else "Group")) +
         ggplot2::facet_grid(col_id ~ .) +
         ggplot2::scale_x_continuous(x, limits = c(0, 1), expand = ggplot2::expansion(add = 0.01), breaks = seq(0, 1, 0.1)) +
         ggplot2::scale_color_brewer("ID", palette = "Set1") +
         ggplot2::theme(panel.background = ggplot2::element_blank(),
                        panel.grid = ggplot2::element_line(color = "gray80"),
                        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                        axis.title.y = ggplot2::element_blank(),
                        axis.text.y = ggplot2::element_blank(),
                        axis.ticks.y = ggplot2::element_blank(),
                        legend.position = "bottom") +
         ggtitle("Bootstrap Distribution")
     })

   }
















