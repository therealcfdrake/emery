
#' S4 object containing the results of multi-method ML accuracy estimates
#' @slot results a list of estimated accuracy statistics
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
    names = "list",
    iter = "numeric",
    prog = "list",
    type = "character"
  ),
  prototype = list(
    results = list(),
    data = matrix(),
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
