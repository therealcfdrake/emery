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

setMethod(
  "show",
  signature(object = "MultiMethodMLEstimate"),
    function(object){
      print(object@results)
    }
  )

setMethod(
  "plot",
  signature(x = "MultiMethodMLEstimate"),
  function(x, ...){
    plot_ML(x, ...)
  })
