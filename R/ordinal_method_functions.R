gen_multi_ord <-
  function(n_method = 3,
           n_obs = 100,
           prev = 0.5,
           D = NULL,
           n_level = 5,
           pmf_pos = matrix(rep(1:n_level - 1, n_method), nrow = n_method, byrow = TRUE),
           pmf_neg = matrix(rep(n_level:1 - 1, n_method), nrow = n_method, byrow = TRUE),
           method_names = NULL,
           level_names = NULL,
           obs_names = NULL
  ){

    if(is.null(method_names)){method_names <- name_thing("method", n_method)}
    if(is.null(level_names)){level_names <- name_thing("level", n_level)}
    if(is.null(obs_names)){obs_names <- name_thing("obs", n_obs)}

    pmf_pos <- pmf_pos / rowSums(pmf_pos)
    pmf_neg <- pmf_neg / rowSums(pmf_neg)

    dimnames(pmf_pos) <- list(method_names, level_names)
    dimnames(pmf_neg) <- list(method_names, level_names)

    dis <- define_disease_state(D, n_obs, prev)

    ord_sample <-
      function(n, pmf){sample(1:length(pmf), n, replace = TRUE, prob = pmf)}

    generated_data <-
      rbind(
        lapply(setNames(1:n_method, method_names), function(x) ord_sample(dis$pos, pmf_pos[x, ])) |> as.data.frame(),
        lapply(setNames(1:n_method, method_names), function(x) ord_sample(dis$neg, pmf_neg[x, ])) |> as.data.frame()
      ) |> as.matrix()

    dimnames(generated_data) <- list(obs_names, method_names)

    se_observed <- list()
    sp_observed <- list()

    for(x in 1:n_level){
      se_observed[[x]] <- colMeans(generated_data[dis$D == 1, ] > x)
      sp_observed[[x]] <- colMeans(generated_data[dis$D == 0, ] < x)
    }

    se_observed <- do.call(cbind, se_observed)
    sp_observed <- do.call(cbind, sp_observed)

    colnames(se_observed) <- level_names
    colnames(sp_observed) <- level_names

    params <-
      list(
        n_method = n_method,
        n_level = n_level,
        n_obs = dis$n_obs,
        prev = dis$prev,
        D = dis$D,
        pmf_pos = pmf_pos,
        pmf_neg = pmf_neg,
        se_observed = se_observed,
        sp_observed = sp_observed,
        method_names = method_names,
        level_names = level_names,
        obs_names = obs_names
      )

    return(
      list(
        generated_data = generated_data,
        params = params
      )
    )

  }

estimate_EM_ordinal <-
  function(data,
           init = list(pi_1_1 = NULL, phi_1ij_1 = NULL, phi_0ij_1 = NULL, n_level = NULL),
           level_names = NULL,
           max_iter = 1000,
           tol = 1e-7,
           save_progress = FALSE){

  calc_l_cond_ordinal <- function(){
    l_cond <-
      sum(c(
        q_k0_t * log(g_0_t),
        q_k1_t * log(g_1_t)
      ))
    return(l_cond)
  }
  calc_A_i <- function(phi_1ij, phi_0ij){
    outer_sum <- 0
    for(j in 1:(n_level - 1)){
      inner_sum <- colSums(phi_1ij[(j + 1):n_level, , drop = FALSE])
      outer_sum <- outer_sum + phi_0ij[j, ] * inner_sum
    }
    A_i <- outer_sum + 0.5 * colSums(phi_1ij * phi_0ij)
    return(A_i)
  }
  calc_y_k <- function(){
    # n_method <- ncol(t_k) -> i
    # n_level               -> j
    # n_obs <- nrow(t_k)    -> k
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
  calc_g_d <- function(phi_dij){
    g_d <- lapply(y_k, function(k) prod(phi_dij ^ k)) |> unlist() |> pmax(1e-300)
    return(g_d)
  }
  calc_q_kd <- function(d){
    q_kd <-
      (p_t * g_1_t * d + (1 - p_t) * g_0_t * (1 - d)) /  # d terms added to toggle numerator
      ((1 - p_t) * g_0_t + p_t * g_1_t)
    return(q_kd)
  }
  calc_next_prev <- function(q_k1){
    mean(q_k1)
  }
  calc_next_phi_dij <- function(q_kd){
    lapply(1:length(q_kd), function(k){q_kd[k] * y_k[[k]]}) |> Reduce(f = "+", x = _) /
      sum(q_kd)
  }

  t_k <- data
  n_method <- ncol(t_k)
  n_obs <- nrow(t_k)
  method_names <- if(is.null(colnames(t_k))){name_thing("method", n_method)}else{colnames(t_k)}
  obs_names <- if(is.null(rownames(t_k))){name_thing("obs", n_obs)}else{rownames(t_k)}

  if(is.null(init$n_level)){n_level <- length(unique(as.vector(t_k)))}
  if(is.null(level_names)){level_names <- name_thing("level", n_level)}
  if(any(sapply(init, is.null))){init <- pollinate_EM_ordinal(t_k, n_level = n_level, level_names = level_names)}

  p_t <- init$pi_1_1
  phi_1ij_t <- init$phi_1ij_1
  phi_0ij_t <- init$phi_0ij_1
  y_k <- calc_y_k()

  list_prev <- list()
  list_phi_1ij <- list()
  list_phi_0ij <- list()
  list_y_k <- list()
  list_g_1 <- list()
  list_g_0 <- list()
  list_q_k1 <- list()
  list_q_k0 <- list()
  list_l_cond <- list()

  for(iter in 1:max_iter){

    g_1_t <- calc_g_d(phi_1ij_t)
    g_0_t <- calc_g_d(phi_0ij_t)
    q_k1_t <- calc_q_kd(d = 1)
    q_k0_t <- calc_q_kd(d = 0)
    l_cond_t <- calc_l_cond_ordinal()
    list_prev <- c(list_prev, list(p_t))
    list_phi_1ij <- c(list_phi_1ij, list(phi_1ij_t))
    list_phi_0ij <- c(list_phi_0ij, list(phi_0ij_t))
    list_y_k <- c(list_y_k, list(y_k))
    list_g_1 <- c(list_g_1, list(g_1_t))
    list_g_0 <- c(list_g_0, list(g_0_t))
    list_q_k1 <- c(list_q_k1, list(q_k1_t))
    list_q_k0 <- c(list_q_k0, list(q_k0_t))
    list_l_cond <- c(list_l_cond, list(l_cond_t))

    if(iter > 1){if(abs(list_l_cond[[iter]] - list_l_cond[[iter - 1]]) < tol){break}}

    p_t <- calc_next_prev(q_k1_t)
    phi_1ij_t <- calc_next_phi_dij(q_kd = q_k1_t)
    phi_0ij_t <- calc_next_phi_dij(q_kd = q_k0_t)

  }

  output <-
    new("MultiMethodMLEstimate",
        results = list(
          prev_est = p_t,
          phi_1ij_est = phi_1ij_t,
          phi_0ij_est = phi_0ij_t,
          q_k1_est = q_k1_t),
        iter = iter,
        type = "ordinal")

  if(save_progress){
    output@prog <-
      list(
        prev = list_prev,
        phi_1ij = list_phi_1ij,
        phi_0ij = list_phi_0ij,
        y_k = list_y_k,
        g_1 = list_g_1,
        g_0 = list_g_0,
        q_k1 = list_q_k1,
        q_k0 = list_q_k0,
        l_cond = list_l_cond)
  }

  return(output)

}

pollinate_EM_ordinal <-
  function(data,
           n_level = NULL,
           threshold_level = ceiling(n_level / 2),
           level_names = NULL){

    t_k <- data
    n_method <- ncol(t_k)
    n_obs <- nrow(t_k)

    if(is.null(n_level)){n_level <- length(unique(as.vector(t_k)))}

    method_names <- if(is.null(colnames(t_k))){name_thing("method", ncol(t_k))}else{colnames(t_k)}
    obs_names <- if(is.null(rownames(t_k))){name_thing("obs", nrow(t_k))}else{rownames(t_k)}

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

# a <- generate_multimethod_data(n_method = 3, type = "ordinal")
# b <- estimate_EM_ordinal(a$generated_data, save_progress = TRUE)
