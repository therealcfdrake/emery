
# d disease status
# i methods
# j response
# k observations


### Build Sample Data Set

#' Title
#'
#' @param n_method number of methods to generate simulated data for
#' @param n_obs number of observations to simulate
#' @param prev true prevalence of positives in target population
#' @param D binary vector representing the true disease status of the observations
#' @param mu_i1 vector of method mean values for positive observations
#' @param sigma_i1 covariance matrix of method positive observations
#' @param mu_i0 vector of method mean values for negative observations
#' @param sigma_i0 covariance matrix of method negative observations
#' @param method_names vector of names of each method
#'
#' @return list containing generated data and input parameters
#' @export
#' @importFrom mvtnorm rmvnorm
#' @examples
#'

gen_multi_con <-
  function(
    n_method = 2,
    n_obs = 100,
    prev = 0.5,
    D = NULL,
    mu_i1 = rep(12, n_method),
    sigma_i1 = diag(n_method),
    mu_i0 = rep(10, n_method),
    sigma_i0 = diag(n_method),
    method_names = NULL
    ){

  pos <- round(n_obs * prev, 0)
  neg <- n_obs - pos

  X <- mvtnorm::rmvnorm(n = pos, mean = mu_i1, sigma = sigma_i1)
  Y <- mvtnorm::rmvnorm(n = neg, mean = mu_i0, sigma = sigma_i0)

  generated_data <- rbind(X, Y)

  D <- c(rep(1, pos), rep(0, neg))

  return(
    list(generated_data = generated_data,
         n_method = n_method,
         n_obs = n_obs,
         D = D,
         mu_i1 = mu_i1,
         sigma_i1 = sigma_i1,
         mu_i0 = mu_i0,
         sigma_i0 = sigma_i0,
         method_names = method_names
         )
  )
}

### Internal Functions

calc_l_cond_continuous <- function(prev_m, z_k1_m, z_k0_m, f_X_m, f_Y_m){

  sum(
    (z_k1_m *
       log((prev_m) * #paper is inconsistent on whether to use log p or log of product
             pmax(f_X_m, 1e-300))) + #small minimum value, algorithm fails is this reaches 0
      (z_k0_m *
         log((1 - prev_m) * #paper is inconsistent on whether to use log p or log of product
               pmax(f_Y_m, 1e-300))) #small minimum value, algorithm fails is this reaches 0
  )

}

calc_z_kd <- function(prev_m, f_X_m, f_Y_m, D){

  #combined numerators into single function. D determines which part is used.
  ((prev_m * f_X_m) * D + ((1 - prev_m) * f_Y_m) * (1 - D)) /
    ((prev_m * f_X_m) + ((1 - prev_m) * f_Y_m))

}

calc_next_prev <- function(z_k1_m){

  mean(z_k1_m)

}

calc_next_mu_id <- function(t_k, z_kd_m){

  #returns vector of mu_id(m+1) estimates

  colSums((z_kd_m * t_k)) /
    sum(z_kd_m)

}

calc_next_sigma_d <- function(prev_m, t_k, z_kd, mu_id){

  # returns covariance matrix sigma_d(m+1) estimate

  n_method <- length(mu_id)
  sigma_d <- matrix(nrow = n_method, ncol = n_method)

  for(i in 1:n_method){

    for(j in 1:i){

      sigma_d[i, j] <-
        sum(z_kd * (t_k[, i] - mu_id[i]) * (t_k[, j] - mu_id[j])) /
        sum(z_kd)

      sigma_d[j, i] <- sigma_d[i, j]

    }

  }
  return(sigma_d)
}

eta_j <- function(mu_i1, sigma_i1, mu_i0, sigma_i0){

  (mu_i1 - mu_i0) / sqrt(diag(sigma_i1) + diag(sigma_i0))

}

A_j <- function(eta){

  test_names <- if(is.null(names(eta))){paste0("test_", 1:length(eta))}else{names(eta)}
  data.frame(test = test_names, A_hat = pnorm(eta))

}

#' Title
#'
#' @param data A matrix with observations as rows and methods as columns.
#' @param init A named list of initial values for prevalence, location, and dispresion of disease groups.
#' @param iter An integer specifying the maximum number of iterations
#'
#' @return
#' @export
#'
#' @examples
estimate_EM_continuous <-
  function(data,
           init = list(
             prev_1 = NULL,
             mu_i1_1 = NULL,
             sigma_i1_1 = NULL,
             mu_i0_1 = NULL,
             sigma_i0_1 = NULL),
           iter = 100){

  if(any(sapply(init, is.null))){init <- pollinate_EM_continuous(data)}

  t_k <- data
  prev_m <- init$prev_1
  mu_i1_m <- init$mu_i1_1
  sigma_i1_m <- init$sigma_i1_1
  mu_i0_m <- init$mu_i0_1
  sigma_i0_m <- init$sigma_i0_1

  list_prev <- list()
  list_mu_i1 <- list()
  list_sigma_i1 <- list()
  list_mu_i0 <- list()
  list_sigma_i0 <- list()

  list_z_k1 <- list()
  list_z_k0 <- list()
  list_l_cond <- list()

  for(j in 1:iter){

    f_X_m <- mvtnorm::dmvnorm(t_k, mean = mu_i1_m, sigma = sigma_i1_m)
    f_Y_m <- mvtnorm::dmvnorm(t_k, mean = mu_i0_m, sigma = sigma_i0_m)

    z_k1_m <- calc_z_kd(prev_m, f_X_m, f_Y_m, D = 1)
    z_k0_m <- calc_z_kd(prev_m, f_X_m, f_Y_m, D = 0)

    l_cond_m <- calc_l_cond_continuous(prev_m, z_k1_m, z_k0_m, f_X_m, f_Y_m)
    list_prev[[j]] <- prev_m
    list_mu_i1[[j]] <- mu_i1_m
    list_sigma_i1[[j]] <- sigma_i1_m
    list_mu_i0[[j]] <- mu_i0_m
    list_sigma_i0[[j]] <- sigma_i0_m

    list_z_k1[[j]] <- z_k1_m
    list_z_k0[[j]] <- z_k0_m
    list_l_cond[[j]] <- l_cond_m

    if(j > 1){if(abs(list_l_cond[[j]] - list_l_cond[[j - 1]]) < .000001){break}}

    prev_m <- calc_next_prev(z_k1_m)
    mu_i1_m <- calc_next_mu_id(t_k, z_k1_m)
    sigma_i1_m <- calc_next_sigma_d(prev_m, t_k, z_k1_m, mu_i1_m)
    mu_i0_m <- calc_next_mu_id(t_k, z_k0_m)
    sigma_i0_m <- calc_next_sigma_d(prev_m, t_k, z_k0_m, mu_i0_m)

  }

  list(list_prev = list_prev,
       list_mu_i1 = list_mu_i1,
       list_sigma_i1 = list_sigma_i1,
       list_mu_i0 = list_mu_i0,
       list_sigma_i0 = list_sigma_i0,
       list_l_cond = list_l_cond,
       list_calc_z_k1 = list_z_k1,
       list_calc_z_k0 = list_z_k0,
       iter = j)

}

### Starting Value Generator

pollinate_EM_continuous <- function(t_k, prev = 0.5, q_seeds = c((1 - prev) / 2, 1 - (prev / 2)), high_pos = TRUE){

  # adjust seeds depending on whether high (default) or low values are associated with "positive" diagnosis
  q_seeds <- sort(q_seeds, decreasing = high_pos)

  muD_mat <- apply(t_k, MARGIN = 2, quantile, q_seeds)

  mu_i1_1 <- unlist(muD_mat[1, ])
  mu_i0_1 <- unlist(muD_mat[2, ])

  sigma_i1_1 <- diag(mu_i1_1 ^ 2)
  sigma_i0_1 <- diag(mu_i0_1 ^ 2)

  list(prev_1 = prev, mu_i1_1 = mu_i1_1, sigma_i1_1 = sigma_i1_1, mu_i0_1 = mu_i0_1, sigma_i0_1 = sigma_i0_1)

}
