#' @slot results a list of estimated accuracy statistics
#'
#' @slot iter an integer number of iterations used by EM algorithm
#' @slot prog a list of interim statistic estimates
#' @slot type a string describing the data type
#'
#' @import methods

setClass(
  "MultiMethodMLEstimate",
  slots = c(
    results = "list",
    iter = "numeric",
    prog = "list",
    type = "character"
  ),
  prototype = list(
    results = list(),
    iter = 0,
    prog = list(),
    type = NA_character_
  )
)

#' Print Final Results Stored in MultiMethodMLEstimate S4 Objects
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

#' Plot Results Stored in MultiMethodMLEstimate S4 Objects
#'
#' @param x MultiMethodMLEstimate.
#'
#' @return a list of ggplot2 plots
#' @export
#'
setMethod(
  "plot",
  signature(x = "MultiMethodMLEstimate"),
  function(x, ...){
    plot_ML(x, ...)
  })
