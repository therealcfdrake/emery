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
      print(object@iter)
    }
  )

setMethod(
  "plot",
  signature(x = "MultiMethodMLEstimate"),
  function(x, ...){
    switch(
      x@type,
      binary = plot_ML_binary(x, ...),
      ordinal = plot_ML_binary(x, ...),
      continuous = plot_ML_binary(x, ...))
  })
