
#' S4 object for storing the results of multi-method ML accuracy estimates
#' @slot results a list of estimated accuracy statistics
#' @slot names a list containing the names of various dimensions
#' @slot data a copy of the data used to generate the estimate
#' @slot iter an integer number of iterations used by EM algorithm
#' @slot prog a list of interim statistic estimates
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

#' Print final results stored in MultiMethodMLEstimate S4 object
#'
#' @param object MultiMethodMLEstimate.
#'
#' @export
#'
setMethod(
  "show",
  signature(object = "MultiMethodMLEstimate"),
    function(object){
      print(object@results)
    }
  )

#' Call `plot_ML` on MultiMethodMLEstimate S4 objects
#'
#' @param x MultiMethodMLEstimate S4 object
#' @param y not used
#' @param ... Additional arguments
#'
#' @return A list of ggplot2 plots
#' @export
#'
setMethod(
  "plot",
  signature(x = "MultiMethodMLEstimate"),
  function(x, ...){
    plot_ML(x, ...)
  })
