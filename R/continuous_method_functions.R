
# d disease status
# i methods
# j response
# k observations

### Build Sample Data Set

#' Generate data set for multiple continuous methods measuring the same binary variable
#'
#' @param n_method number of methods to generate simulated data for
#' @param n_obs number of observations to simulate
#' @param prev true prevalence of positives in target population
#' @param D binary vector representing the true disease status of the observations
#' @param mu_i1 vector of method mean values for positive observations
#' @param sigma_i1 covariance matrix of method positive observations
#' @param mu_i0 vector of method mean values for negative observations
#' @param sigma_i0 covariance matrix of method negative observations
#' @param obs_names vector of names of each observation
#'
#' @return list containing generated data and input parameters
#' @export
#' @importFrom mvtnorm rmvnorm
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
    method_names = NULL,
    obs_names = NULL
    ){

  dis <- define_disease_state(D, n_obs, prev)

  X <- mvtnorm::rmvnorm(n = dis$pos, mean = mu_i1, sigma = sigma_i1)
  Y <- mvtnorm::rmvnorm(n = dis$neg, mean = mu_i0, sigma = sigma_i0)

  generated_data <- rbind(X, Y)

  params <-
    list(
      n_method = n_method,
      n_obs = dis$n_obs,
      prev = dis$prev,
      D = dis$D,
      mu_i1 = mu_i1,
      sigma_i1 = sigma_i1,
      mu_i0 = mu_i0,
      sigma_i0 = sigma_i0,
      method_names = method_names,
      obs_names = obs_names
    )

  return(
    list(generated_data = generated_data,
         params = params
         )
  )
}

#' EM for continuous methods
#'
#' @param data A matrix with observations as rows and methods as columns.
#' @param max_iter Maximum number of iterations to stop seeking convergence of l estimate.
#' @param tol Change in l below which l is said to converge.
#' @param save_progress Logical value denoting whether to output the calculated values of each iteration.
#' @param init A named list of initial values for prevalence, location, and dispersion of disease groups.
#'
#' @return A list containing the final calculated values.
#' @export
#'
estimate_EM_continuous <-
  function(data,
           init = list(
             prev_1 = NULL,
             mu_i1_1 = NULL,
             sigma_i1_1 = NULL,
             mu_i0_1 = NULL,
             sigma_i0_1 = NULL),
           max_iter = 100,
           tol = 1e-7,
           save_progress = FALSE){

  calc_l_cond_continuous <- function(){
      sum(
        (z_k1_m *
           log((prev_m) * #paper is inconsistent on whether to use log p or log of product
                 pmax(f_X_m, 1e-300))) + #small minimum value, algorithm fails is this reaches 0
          (z_k0_m *
             log((1 - prev_m) * #paper is inconsistent on whether to use log p or log of product
                   pmax(f_Y_m, 1e-300))) #small minimum value, algorithm fails is this reaches 0
      )
    }
  calc_z_kd <- function(d){
      #combined numerators into single function. d determines which part is used.
      ((prev_m * f_X_m) * d + ((1 - prev_m) * f_Y_m) * (1 - d)) /
        ((prev_m * f_X_m) + ((1 - prev_m) * f_Y_m))
    }
  calc_next_prev <- function(){
      mean(z_k1_m)
    }
  calc_next_mu_id <- function(z_kd_m){
      #returns vector of mu_id(m+1) estimates
      colSums((z_kd_m * t_k)) /
        sum(z_kd_m)
    }
  calc_next_sigma_d <- function(z_kd_m, mu_id_m){
      # returns covariance matrix sigma_d(m+1) estimate
      n_method <- length(mu_id_m)
      sigma_d <- matrix(nrow = n_method, ncol = n_method)
      for(i in 1:n_method){
        for(j in 1:i){
          sigma_d[i, j] <-
            sum(z_kd_m * (t_k[, i] - mu_id_m[i]) * (t_k[, j] - mu_id_m[j])) /
            sum(z_kd_m)
          sigma_d[j, i] <- sigma_d[i, j]
        }
      }
      return(sigma_d)
    }
  calc_eta_j <- function(){
    (mu_i1_m - mu_i0_m) / sqrt(diag(sigma_i1_m) + diag(sigma_i0_m))
  }
  calc_A_j <- function(){
    setNames(pnorm(eta_j_m), method_names)
  }

  t_k <- data
  n_method <- ncol(t_k)
  method_names <- if(is.null(colnames(t_k))){name_thing("method", n_method)}else{colnames(t_k)}

  if(any(sapply(init, is.null))){init <- pollinate_EM_continuous(t_k)}
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

  list_eta_j <- list()
  list_A_j <- list()
  list_z_k1 <- list()
  list_z_k0 <- list()
  list_l_cond <- list()

  for(iter in 1:max_iter){

    f_X_m <- mvtnorm::dmvnorm(t_k, mean = mu_i1_m, sigma = sigma_i1_m)
    f_Y_m <- mvtnorm::dmvnorm(t_k, mean = mu_i0_m, sigma = sigma_i0_m)

    eta_j_m <- calc_eta_j()
    A_j_m <- calc_A_j()
    z_k1_m <- calc_z_kd(d = 1)
    z_k0_m <- calc_z_kd(d = 0)
    l_cond_m <- calc_l_cond_continuous()

    list_prev <- c(list_prev, list(prev_m))
    list_mu_i1 <- c(list_mu_i1, list(mu_i1_m))
    list_sigma_i1 <- c(list_sigma_i1, list(sigma_i1_m))
    list_mu_i0 <- c(list_mu_i0, list(mu_i0_m))
    list_sigma_i0 <- c(list_sigma_i0, list(sigma_i0_m))
    list_eta_j <- c(list_eta_j, list(eta_j_m))
    list_A_j <- c(list_A_j, list(A_j_m))
    list_z_k1 <- c(list_z_k1, list(z_k1_m))
    list_z_k0 <- c(list_z_k0, list(z_k0_m))
    list_l_cond <- c(list_l_cond, list(l_cond_m))

    if(iter > 1){if(abs(list_l_cond[[iter]] - list_l_cond[[iter - 1]]) < tol){break}}

    prev_m <- calc_next_prev()
    mu_i1_m <- calc_next_mu_id(z_k1_m)
    sigma_i1_m <- calc_next_sigma_d(z_k1_m, mu_i1_m)
    mu_i0_m <- calc_next_mu_id(z_k0_m)
    sigma_i0_m <- calc_next_sigma_d(z_k0_m, mu_i0_m)

  }

  output <-
    new("MultiMethodMLEstimate",
        results = list(
          prev_est = prev_m,
          mu_i1_est = mu_i1_m,
          sigma_i1_est = sigma_i1_m,
          mu_i0_est = mu_i0_m,
          sigma_i0_est = sigma_i0_m,
          eta_j_est = eta_j_m,
          A_j_est = A_j_m,
          z_k1_est = z_k1_m,
          z_k0_est = z_k0_m),
        iter = iter,
        type = "continuous")

  if(save_progress){
    output@prog <-
      list(
        prev = list_prev,
        mu_i1 = list_mu_i1,
        sigma_i1 = list_sigma_i1,
        mu_i0 = list_mu_i0,
        sigma_i0 = list_sigma_i0,
        eta_j = list_eta_j,
        A_j = list_A_j,
        z_k1 = list_z_k1,
        z_k0 = list_z_k0,
        l_cond = list_l_cond)
  }

  return(output)

  }

### Starting Value Generator

#' Create initialization values for EM estimation
#'
#' @param t_k A data set of `n_obs` rows and `n_method` columns
#' @param prev A double between 0-1 representing the proportion of positives in the population
#' @param q_seeds A vector of length 2 representing the quantiles at which the two groups are assumed to be centered
#' @param high_pos A logical indicating whether larger values are considered "positive"
#'
#' @return
#' @export
#'
#' @examples
pollinate_EM_continuous <-
  function(t_k,
           prev = 0.5,
           q_seeds = c((1 - prev) / 2, 1 - (prev / 2)),
           high_pos = TRUE){

  # adjust seeds depending on whether high (default) or low values are associated with "positive" diagnosis
  q_seeds <- sort(q_seeds, decreasing = high_pos)

  muD_mat <- apply(t_k, MARGIN = 2, quantile, q_seeds)

  mu_i1_1 <- unlist(muD_mat[1, ])
  mu_i0_1 <- unlist(muD_mat[2, ])

  sigma_i1_1 <- diag(mu_i1_1 ^ 2)
  sigma_i0_1 <- diag(mu_i0_1 ^ 2)

  return(
    list(
      prev_1 = prev,
      mu_i1_1 = mu_i1_1,
      sigma_i1_1 = sigma_i1_1,
      mu_i0_1 = mu_i0_1,
      sigma_i0_1 = sigma_i0_1)
  )

}

# a <- generate_multimethod_data(n_method = 3, type = "continuous")
# b <- estimate_EM_continuous(a$generated_data, save_progress = TRUE)
