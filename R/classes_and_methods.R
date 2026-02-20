
#' S4 object containing the results of multi-method ML accuracy estimates
#' @slot results a list of estimated accuracy statistics
#' @slot data a matrix containing the raw data used for estimation
#' @slot freqs a vector containing the frequencies at which each row in `data` was observed
#' @slot names a list containing vectors of names of various dimensions
#' @slot data a copy of the data used to generate the estimated values
#' @slot iter an integer number of iterations needed for the EM algorithm to converge
#' @slot prog a list containing the values calculated during each iteration of
#' the EM algorithm
#' @slot type a string describing the data type
#'
#' @import methods

setClass(
  "MultiMethodMLEstimate",
  slots = c(
    results = "list",
    data = "matrix",
    freqs = "numeric",
    names = "list",
    iter = "numeric",
    prog = "list",
    type = "character"
  ),
  prototype = list(
    results = list(),
    data = matrix(),
    freqs = c(),
    names = list(),
    iter = 0,
    prog = list(),
    type = NA_character_
  )
)

#' Show a MultiMethodMLEstimate S4 object
#' @description
#' Print the accuracy statistic estimates stored in a MultiMethodMLEstimate object.
#'
#' @param object An object of class MultiMethodMLEstimate.
#' @return A list containing relevant accuracy statistic estimates. This is a
#' subset of the list stored in `results` slot of the MultiMethodMLEstimate object.
#' @export
#'
setMethod(
  "show",
  signature(object = "MultiMethodMLEstimate"),
    function(object){
      print(object@results[!grepl("q|z", names(object@results))])
    }
  )

#' Create plots from a MultiMethodMLEstimate object
#' @description
#' Create a list of plots visualizing the expectation maximization process and
#' resulting accuracy statistics stored in a MultiMethodMLEstimate object.
#' @inheritDotParams plot_ML params
#' @param x a MultiMethodMLEstimate S4 object
#' @param y not used
#' @param ... Additional arguments
#' @return A list of ggplot2 plots
#' @export
#'
setMethod(
  "plot",
  signature(x = "MultiMethodMLEstimate"),
  function(x, ...){
    plot_ML(x, ...)
  })

#' Create new boot_ML class object
#' @description
#' Wrapper for creating boot_ML class object.
#' @param n_obs Number of observations in data
#' @inheritParams boot_ML
#' @inheritParams estimate_ML
#' @param v_0 MultiMethodMLEstimate S4 object
#' @param v_star results slot of bootstrapped MultiMethodMLEstimate objects
#' @return a boot_ML object
#'
new_boot_ML <-
  function(v_0, v_star, data, freqs = NULL, n_boot, n_study, max_iter, tol, n_obs, seed){
    structure(
      .Data = list(
        v_0 = v_0,
        v_star = v_star,
        params = list(
          data = data,
          freqs = freqs,
          n_boot = n_boot,
          n_study = n_study,
          max_iter = max_iter,
          tol = tol,
          n_obs = n_obs,
          seed = seed
        )
      ),
      class = "boot_ML"
    )
  }

setGeneric("getResults", function(x){
  standardGeneric("getResults")
  })

setMethod("getResults", signature("MultiMethodMLEstimate"), function(x){
  x@results
})

setGeneric("getNames", function(x, name){
  standardGeneric("getNames")
  })

setMethod("getNames", signature("MultiMethodMLEstimate"), function(x, name){
  x@names[name]
})

