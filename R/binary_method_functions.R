## Binary

gen_multi_bin <-
  function(n_method = 3,
           n_obs = 100,
           prev = 0.5,
           D = NULL,
           se = rep(0.9, n_method),
           sp = rep(0.9, n_method),
           method_names = NULL,
           n_method_subset = n_method,
           first_reads_all = FALSE
  ){

    # "True" values have no suffix and are provided by user (population parameters)
    # "*_majority" - calculated using observed values with the majority class as the reference
    # "*_observed" - calculated using observed values with the "True" value as the reference

    if(is.null(method_names)){method_names <- thing_namer("method", n_method)}

    pos <- round(n_obs * prev, 0)
    neg <- n_obs - pos

    if(is.null(D)){D <- c(rep(1, pos), rep(0, neg))}

    # Create dataframe to (optionally) randomly censor observations such that only "n_method_subset" results are reported for each observation
    subset_matrix <-
      lapply(1:n_obs, function(i)
        if(!first_reads_all){sample(c(rep(1, n_method_subset), rep(NA, n_method - n_method_subset)), n_method, replace = FALSE)
        }else{
          c(1, sample(c(rep(1, n_method_subset - 1), rep(NA, n_method - n_method_subset)), n_method - 1, replace = FALSE))
        }
      ) |>
      do.call(what = rbind, args = _)

    # Build simulated data set based on input criteria
    generated_data <-
      lapply(1:n_method, function(i) rbinom(n_obs, 1, se[i] * D + (1 - sp[i]) * (1 - D))) |>
      do.call(what = cbind, args = _) * subset_matrix |>
      as.data.frame() |>
      setNames(method_names)

    # Calculate individual test se based on input criteria. This will differ slightly from "se" due to random sampling.
    se_observed <-
      cbind(generated_data, D = D) |>
      subset(D == 1, select = -D) |>
      colMeans(na.rm = TRUE)

    # Calculate individual test sp based on input criteria. This will differ slightly from "sp" due to random sampling.
    sp_observed <-
      cbind(generated_data, D = D) |>
      subset(D == 0, select = -D) |>
      colMeans(na.rm = TRUE) |>
      (\(x) 1 - x)()

    # Determine the majority classification. A small random number is added to each in order to break ties in cases with an even number of tests.
    D_majority <- generated_data |>
      rowMeans(na.rm = TRUE) |>
      (\(x) x + runif(n_obs, -0.000001, 0.000001))() |>
      round()

    # # Estimate pi wrt majority classification
    # prev_majority <- mean(D_majority, na.rm = TRUE)
    #
    # # Estimate individual test se wrt majority classification
    # se_majority <- generated_data[D_majority == 1, ] |> colMeans(na.rm = TRUE)
    #
    # # Estimate individual test sp wrt majority classification
    # sp_majority <- 1 - generated_data[D_majority == 0, ] |> colMeans(na.rm = TRUE)

    return(
      list(
        generated_data = generated_data,
        n_method = n_method,
        n_obs = n_obs,
        prev = prev,
        se = se,
        sp = sp,
        D = D,
        method_names = method_names,
        se_observed = se_observed,
        sp_observed = sp_observed
        )
    )
  }

pollinate_EM_binary <- function(data){

  method_names <- if(is.null(names(data))){thing_namer("method", ncol(data))}else{names(data)}

  n_obs <- nrow(data)

  D_majority <- data |>
    rowMeans(na.rm = TRUE) |>
    (\(x) x + runif(n_obs, -0.000001, 0.000001))() |>
    round()

  # Estimate prev_1 wrt majority classification
  prev_1 <- mean(D_majority, na.rm = TRUE) |> setNames("prev")

  # Estimate individual test se wrt majority classification
  se_1 <- data[D_majority == 1, ] |> colMeans(na.rm = TRUE) |> setNames(method_names)

  # Estimate individual test sp wrt majority classification
  sp_1 <- 1 - (data[D_majority == 0, ] |> colMeans(na.rm = TRUE) |> setNames(method_names))

  return(
    list(prev_1 = prev_1,
         se_1 = se_1,
         sp_1 = sp_1)
  )

}


calc_A2 <- function(data, se_m, prev_m){

  data |>
    apply(1, FUN = function(y) (se_m ^ y) * ((1 - se_m) ^ (1 - y))) |>
    apply(2, FUN = function(z) prod(z, na.rm = TRUE)) |>
    (\(x) x * prev_m)()

}

calc_B2 <- function(data, sp_m, prev_m){

  data |>
    apply(1, FUN = function(y) ((1 - sp_m) ^ y) * (sp_m ^ (1 - y))) |>
    apply(2, FUN = function(z) prod(z, na.rm = TRUE)) |>
    (\(x) x * (1 - prev_m))()

}


calc_qk <- function(A2, B2){

  A2 / (A2 + B2)

}


calc_next_se <- function(data, qk_m){

  data_mat <- as.matrix(data)
  data_mat[is.na(data_mat)] <- 0 # added to handle NA, used to calc appropriate denominators

  (qk_m %*% data_mat) / (qk_m %*% !is.na(as.matrix(data)))

  # (qk_m %*% as.matrix(data)) / sum(qk_m) # original function before changes to accept NA

}

calc_next_sp <- function(data, qk_m){

  data_mat <- 1 - as.matrix(data)
  data_mat[is.na(data_mat)] <- 0 # added to handle NA, used to calc appropriate denominators

  ((1 - qk_m) %*% data_mat) / ((1 - qk_m) %*% !is.na(as.matrix(data)))

  # ((1 - qk_m) %*% (1 - as.matrix(data))) / sum(1 - qk_m) # original function before changes to accept NA

}

calc_next_prev <- function(data, qk_m){

  mean(qk_m)

}

EM_estimator_binary <- function(data, init = list(prev_1 = NULL, se_1 = NULL, sp_1 = NULL), n_iter = 100, tol = 1e-7){

  if(any(sapply(init, is.null))){init <- pollinate_EM_binary(data)}

  method_names <- if(is.null(names(data))){thing_namer("method", ncol(data))}else{names(data)}

  # starting values
  se_m <- init$se_1
  sp_m <- init$sp_1
  prev_m <- init$prev_1

  # initialize lists
  list_se <- list(se_m)
  list_sp <- list(sp_m)
  list_prev <- list(prev_m)

  converged <- FALSE
  iter <- 1

  # iterate

  while(!converged & iter <= n_iter){

    A2_m <- calc_A2(data, se_m = se_m, prev_m)
    B2_m <- calc_B2(data, sp_m = sp_m, prev_m)
    qk_m <- calc_qk(A2_m, B2_m)

    `se_m+1` <- calc_next_se(data, qk_m)
    `sp_m+1` <- calc_next_sp(data, qk_m)
    `prev_m+1` <- calc_next_prev(data, qk_m)

    list_se <- c(list_se, list(`se_m+1`))
    list_sp <- c(list_sp, list(`sp_m+1`))
    list_prev <- c(list_prev, list(`prev_m+1`))

    if(max(abs(`se_m+1` - se_m), abs(`sp_m+1` - sp_m), abs(`prev_m+1` - prev_m)) < tol){converged <- TRUE}

    se_m <- `se_m+1`
    sp_m <- `sp_m+1`
    prev_m <- `prev_m+1`

    iter <- iter + 1

  }

  df_prev <- as.data.frame(do.call(rbind, list_prev))
  df_prev$iter <- row(df_prev)

  df_se <-
  reshape(
    as.data.frame(do.call(rbind, list_se)),
    direction = "long",
    varying = list(method_names),
    times = method_names,
    timevar = "method",
    v.names = "se",
    idvar = "iter")

  df_sp <-
  reshape(
    as.data.frame(do.call(rbind, list_sp)),
    direction = "long",
    varying = list(method_names),
    times = method_names,
    timevar = "method",
    v.names = "sp",
    idvar = "iter")

  estimates <-
    df_se |>
      merge(df_sp) |>
      merge(df_prev)

  return(estimates)

}






