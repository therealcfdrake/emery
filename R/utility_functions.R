
#' Create unique names for a set of things
#'
#' @param thing a string that describes the set of items to name
#' @param n an integer number of unique names to create
#'
#' @return a vector of unique names
#' @importFrom stringr str_pad
#' @export

name_thing <-
  function(thing = "", n = 1){

    paste0(thing, stringr::str_pad(1:n, width = nchar(n), side = "left", pad = "0"))

  }

#' Define the True disease state of a simulated sample
#'
#' @inheritParams generate_multimethod_data
#'
#' @return A list of features defining the true disease status of each observation

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

#' Censor data randomly rowwise
#'
#' @inheritParams generate_multimethod_data

censor_data <- function(
    n_obs = dis$n_obs,
    first_reads_all = first_reads_all,
    n_method_subset = n_method_subset,
    n_method = n_method){

  lapply(1:n_obs, function(i)
    if(!first_reads_all){
      sample(c(rep(1, n_method_subset), rep(NA, n_method - n_method_subset)), n_method, replace = FALSE)
    }else{
      c(1, sample(c(rep(1, n_method_subset - 1), rep(NA, n_method - n_method_subset)), n_method - 1, replace = FALSE))
    }
  ) |> do.call(what = rbind, args = _)

}

#' @title Reduce Data by Summarizing Observations
#' @description
#' Extract unique rows from a data frame. Collect row names of duplicate observations and frequencies of each.
#'
#' @inheritParams estimate_ML
#' @return list containing matrix of unique rows with names and frequencies of duplicates.
#' @export

unique_obs_summary <- function(data){

  unique_obs <- unique(as.matrix(data))
  duplicate_obs <- lapply(1:nrow(unique_obs),
                          function(r){rownames(data[sapply(1:nrow(data), function(s){
                            identical(data[s, ], unique_obs[r, ])}), , drop = FALSE])})
  obs_freq <- vapply(duplicate_obs, length, 1)

  list(
    unique_obs = unique_obs,
    duplicate_obs = duplicate_obs,
    obs_freq = obs_freq
  )

}


#' @title Multivariate Normal Densities
#' @description
#' Return the density of a point in a multivariate normal distribution
#'
#' @param x matrix of observations
#' @param mu vector of method means
#' @param sigma method covariance matrix

dmvnorm <-
  function(x, mu, sigma){
    x_minus_mu <- t(x) - as.vector(mu)
    k <- ncol(x)
    v <- exp(-1 / 2 * t(x_minus_mu) %*% solve(sigma) %*% x_minus_mu) /
      sqrt((2*pi)^k * det(sigma))
    return(diag(v))
  }


# `%mmult%` <-
#   function(A, B){
#     Ar <- nrow(A)
#     apply(B, 2, function(Bc)
#       apply(A, 1, function(Ar)
#         sum(Ar * Bc, na.rm = TRUE)))
#   }
