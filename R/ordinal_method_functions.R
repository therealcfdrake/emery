gen_multi_ord <-
  function(n_method = 3,
           n_level = 5,
           n_obs = 100,
           prev = 0.5,
           D = NULL,
           pmf_pos = matrix(rep(1:n_level - 1, n_method), nrow = n_method, byrow = TRUE),
           pmf_neg = matrix(rep(n_level:1 - 1, n_method), nrow = n_method, byrow = TRUE),
           method_names = NULL,
           level_names = NULL
  ){

    if(is.null(method_names)){method_names <- name_thing("method", n_method)}
    if(is.null(level_names)){level_names <- name_thing("level", n_level)}

    pmf_pos <- pmf_pos / rowSums(pmf_pos)
    pmf_neg <- pmf_neg / rowSums(pmf_neg)

    dimnames(pmf_pos) <- list(method_names, level_names)
    dimnames(pmf_neg) <- list(method_names, level_names)

    # dis <- define_disease_state(D, n_obs, prev)

    pos <- round(n_obs * prev, 0)
    neg <- n_obs - pos
    if (is.null(D)){D <- c(rep(1, pos), rep(0, neg))} else {pos = sum(D); neg = sum(1 - D)}

    ord_sample <-
      function(n, pmf){sample(1:length(pmf), n, replace = TRUE, prob = pmf)}

    generated_data <-
      rbind(
        lapply(setNames(1:n_method, method_names), function(x) ord_sample(pos, pmf_pos[x, ])) |> as.data.frame(),
        lapply(setNames(1:n_method, method_names), function(x) ord_sample(neg, pmf_neg[x, ])) |> as.data.frame()
      ) |> as.matrix()

    se_observed <- list()
    sp_observed <- list()

    for(x in 1:n_level){
      se_observed[[x]] <- colMeans(generated_data[D == 1, ] >= x)
      sp_observed[[x]] <- colMeans(generated_data[D == 0, ] < x)
    }

    return(
      list(generated_data = generated_data,
           n_method = n_method,
           n_level = n_level,
           n_obs = n_obs,
           D = D,
           pmf_pos = pmf_pos,
           pmf_neg = pmf_neg,
           method_names = method_names,
           level_names = level_names,
           se_observed = as.data.frame(se_observed, col.names = level_names),
           sp_observed = as.data.frame(sp_observed, col.names = level_names))
    )


  }


pollinate_EM_ordinal <- function(t_k, n_level, threshold_level = ceiling(n_level / 2), level_names = NULL){

    n_method <- ncol(t_k)
    n_obs <- nrow(t_k)

    method_names <- if(is.null(names(t_k))){name_thing("method", n_method)}else{names(t_k)}

    if(is.null(level_names)){level_names <- name_thing("level", n_level)}

    D_majority <- as.numeric(
      rowMeans(t_k, na.rm = TRUE) |>
        (\(x) x + runif(n_obs, -0.000001, 0.000001))() |> # break ties
        round() |>
        (\(x) x >= threshold_level)()
      )

    pi_1_1 <- mean(D_majority)

    t_k1 <- t_k[D_majority == 1, ]
    t_k0 <- t_k[D_majority == 0, ]

    phi_1ij_1 <-
      lapply(1:n_level, function(j) colMeans(t_k1 == j, na.rm = TRUE)) |>
      unlist() |>
      pmax(1e-10) |>
      pmin((1 - 1e-10)) |>
      matrix(nrow = n_level, ncol = n_method, byrow = TRUE, dimnames = list(level_names, method_names))

    phi_0ij_1 <-
      lapply(1:n_level, function(j) colMeans(t_k0 == j, na.rm = TRUE)) |>
      unlist() |>
      pmax(1e-10) |>
      pmin((1 - 1e-10)) |>
      matrix(nrow = n_level, ncol = n_method, byrow = TRUE, dimnames = list(level_names, method_names))

    return(
      list(
        pi_1_1 = pi_1_1,
        phi_1ij_1 = phi_1ij_1,
        phi_0ij_1 = phi_0ij_1,
        n_level = n_level)
    )

}


calc_l_cond_ordinal <- function(q_k1, q_k0, g_1, g_0){

  l_cond <-
    sum(c(
      q_k0 * log(g_0),
      q_k1 * log(g_1)
    ))

  return(l_cond)

}

calc_A_i <- function(phi_1ij, phi_0ij, n_level){

  outer_sum <- 0

  for(j in 1:(n_level - 1)){

    inner_sum <- colSums(phi_1ij[(j + 1):n_level, , drop = FALSE])

    outer_sum <- outer_sum + phi_0ij[j, ] * inner_sum

  }

  A_i <- outer_sum + 0.5 * colSums(phi_1ij * phi_0ij)
  return(A_i)

}


calc_y_k <- function(t_k, n_level){

  method_names <- if(is.null(names(t_k))){name_thing("method", ncol(t_k))}else{names(t_k)}
  level_names <- name_thing("level", n_level)

  n_method <- ncol(t_k) # -> i
# n_level                 -> j
  n_obs <- nrow(t_k)   #  -> k

  y_k <- list()

  for(k in 1:n_obs){
    tmp_y_k <- matrix(nrow = n_level, ncol = n_method, dimnames = list(level_names, method_names))
    for(j in 1:n_level){
      for(i in 1:n_method){

        tmp_y_k[j, i] <- as.numeric(t_k[k, i] == j)

      }

    }

    y_k[[k]] <- tmp_y_k

  }

  return(y_k)

}


calc_g_d <- function(phi_dij, y_k){

  g_d <- lapply(y_k, function(k) prod(phi_dij ^ k)) |> unlist() |> pmax(1e-300)

  return(g_d)

}


calc_q_kd <- function(g_1, g_0, prev, d){

  q_kd <-
    (prev * g_1 * d + (1 - prev) * g_0 * (1 - d)) /  # d terms added to toggle numerator
      ((1 - prev) * g_0 + prev * g_1)

  return(q_kd)

}


calc_next_pi_1 <- function(q_k1){

  mean(q_k1)

}


calc_next_phi_dij <- function(q_kd, y_k){

  lapply(1:length(q_kd), function(k){q_kd[k] * y_k[[k]]}) |> Reduce(f = "+", x = _) /
    sum(q_kd)

}


estimate_EM_ordinal <- function(data, init = list(pi_1_1 = NULL, phi_1ij_1 = NULL, phi_0ij_1 = NULL, n_level = NULL), n_iter = 1000){

  t_k <- data

  p_t <- pi_1_1
  phi_1ij_t <- phi_1ij_1
  phi_0ij_t <- phi_0ij_1
  y_k <- calc_y_k(t_k, n_level = n_level)

  list_p <- list()
  list_phi_1ij <- list()
  list_phi_0ij <- list()
  list_g_1 <- list()
  list_g_0 <- list()
  list_q_k1 <- list()
  list_q_k0 <- list()
  list_l_cond <- list()

  for(iter in 1:n_iter){


    g_1_t <- calc_g_d(phi_1ij_t, y_k = y_k)
    g_0_t <- calc_g_d(phi_0ij_t, y_k = y_k)

    q_k1_t <- calc_q_kd(g_1_t, g_0_t, p_t, d = 1)
    q_k0_t <- calc_q_kd(g_1_t, g_0_t, p_t, d = 0)

    l_cond_t <- calc_l_cond_ordinal(q_k1 = q_k1_t, q_k0 = q_k0_t, g_1 = g_1_t, g_0 = g_0_t)

    list_p[[iter]] <- p_t
    list_phi_1ij[[iter]] <- phi_1ij_t
    list_phi_0ij[[iter]] <- phi_0ij_t
    list_g_1[[iter]] <- g_1_t
    list_g_0[[iter]] <- g_0_t
    list_q_k1[[iter]] <- q_k1_t
    list_q_k0[[iter]] <- q_k0_t
    list_l_cond[[iter]] <- l_cond_t

    p_t <- calc_next_pi_1(q_k1_t)
    phi_1ij_t <- calc_next_phi_dij(q_kd = q_k1_t, y_k = y_k)
    phi_0ij_t <- calc_next_phi_dij(q_kd = q_k0_t, y_k = y_k)

    if(iter > 1){if(abs(list_l_cond[[iter]] - list_l_cond[[iter - 1]]) < 1e-20){break}}

  }

  list(list_p = list_p,
       list_phi_1ij = list_phi_1ij,
       list_phi_0ij = list_phi_0ij,
       list_g_1 = list_g_1,
       list_g_0 = list_g_0,
       list_q_k1 = list_q_k1,
       list_q_k0 = list_q_k0,
       list_l_cond = list_l_cond,
       y_k = y_k,
       iter = iter)

}
