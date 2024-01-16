
#' Create unique names for a set of things
#'
#' @param thing a string that describes the set of items to name
#' @param n an integer number of unique names to create
#'
#' @return a vector of unique names
#'
name_thing <-
  function(thing = "", n = 1){

    paste0(thing, stringr::str_pad(1:n, width = nchar(n), side = "left", pad = "0"))

  }

#' Define the True disease state of a simulated sample
#'
#' @param D True disease state
#' @param n_obs integer number of observations
#' @param prev prevalence of positive observations in population
#'
#' @return A list of features defining the true disease status of each observation
#'
define_disease_state <-
  function(D = NULL, n_obs = NULL, prev = NULL){

    if(is.null(D)){
      if(!is.numeric(prev)|!is.numeric(n_obs)){
        stop("Invalid disease state specified", call. = FALSE)
      }else{
        pos <- round(n_obs * prev, 0)
        neg <- n_obs - pos
        D <- c(rep(1, pos), rep(0, neg))
      }
    }else{
      n_obs <- length(D)
      pos <- sum(D)
      neg <- sum(1 - D)
      prev <- mean(D)
    }

    return(
      list(
        D = D,
        n_obs = n_obs,
        prev = prev,
        pos = pos,
        neg = neg
      )
    )

  }


