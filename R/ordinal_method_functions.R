
#' @rdname generate_multimethod_data
#' @order 3
#' @export
#' @param n_level Used for ordinal methods. An integer representing the number of ordinal levels each method has
#' @param pmf_pos,pmf_neg Used for ordinal methods. A n_method by n_level matrix representing the probability mass functions for positive and negative results, respectively
#' @param level_names Used for ordinal methods. Optional vector of names used to identify each level
#'
generate_multimethod_ordinal <-
  function(
    n_method = 3,
    n_obs = 100,
    prev = 0.5,
    D = NULL,
    n_level = 5,
    pmf_pos = matrix(rep(1:n_level - 1, n_method), nrow = n_method, byrow = TRUE),
    pmf_neg = matrix(rep(n_level:1 - 1, n_method), nrow = n_method, byrow = TRUE),
    method_names = NULL,
    level_names = NULL,
    obs_names = NULL,
    n_method_subset = n_method,
    first_reads_all = FALSE
  ){

    # Internal functions
    ord_sample <- function(n, pmf){
      sample(1:length(pmf), n, replace = TRUE, prob = pmf)
      }

    if(is.null(method_names)){method_names <- name_thing("method", n_method)}
    if(is.null(level_names)){level_names <- name_thing("level", n_level)}
    if(is.null(obs_names)){obs_names <- name_thing("obs", n_obs)}

    pmf_pos <- pmf_pos / rowSums(pmf_pos)
    pmf_neg <- pmf_neg / rowSums(pmf_neg)

    dimnames(pmf_pos) <- list(method_names, level_names)
    dimnames(pmf_neg) <- list(method_names, level_names)

    dis <- define_disease_state(D, n_obs, prev)

    subset_matrix <-
      censor_data(
        n_obs = dis$n_obs,
        first_reads_all = first_reads_all,
        n_method_subset = n_method_subset,
        n_method = n_method)

    generated_data <-
      rbind(
        lapply(setNames(1:n_method, method_names), function(x) ord_sample(dis$pos, pmf_pos[x, ])) |> as.data.frame(),
        lapply(setNames(1:n_method, method_names), function(x) ord_sample(dis$neg, pmf_neg[x, ])) |> as.data.frame()
      ) |> as.matrix() * subset_matrix

    dimnames(generated_data) <- list(obs_names, method_names)

    se_observed <- list()
    sp_observed <- list()

    for(x in 1:n_level){
      se_observed[[x]] <- colMeans(generated_data[dis$D == 1, ] > x, na.rm = TRUE)
      sp_observed[[x]] <- colMeans(generated_data[dis$D == 0, ] < x, na.rm = TRUE)
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
        D = stats::setNames(dis$D, obs_names),
        pmf_pos = pmf_pos,
        pmf_neg = pmf_neg,
        se_observed = se_observed, # keep?
        sp_observed = sp_observed, # keep?
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

#' @rdname estimate_ML
#' @order 3
#' @export
#' @param level_names An optional, ordered, character vector of unique names corresponding to
#' the levels of the methods.
#'
estimate_ML_ordinal <-
  function(data,
           init = list(
             pi_1_1 = NULL,
             phi_1ij_1 = NULL,
             phi_0ij_1 = NULL,
             n_level = NULL),
           level_names = NULL,
           max_iter = 1000,
           tol = 1e-7,
           save_progress = TRUE){

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
          tmp_y_k[j, i] <-
            as.numeric(
              (!is.na(t_k[k, i])) & # addition to tolerate missing values
                t_k[k, i] == j # original code
              )
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
    denom <- as.vector(q_kd %*% !is.na(t_k)) # missing value correction. only observations which a method had a response are summed
    t(
      lapply(1:length(q_kd), function(k){q_kd[k] * y_k[[k]]}) |>
        Reduce(f = "+", x = _) |>
        t() / denom
    )
      # sum(q_kd) # original denominator
  }

  t_k <- as.matrix(data)
  n_method <- ncol(t_k)
  n_obs <- nrow(t_k)
  method_names <- if(is.null(colnames(t_k))){name_thing("method", n_method)}else{colnames(t_k)}
  obs_names <- if(is.null(rownames(t_k))){name_thing("obs", n_obs)}else{rownames(t_k)}

  dimnames(t_k) <- list(obs_names, method_names)

  if(is.null(init$n_level)){n_level <- sum(!is.na(unique(as.vector(t_k))))}else{n_level <- init$n_level}
  if(is.null(level_names)){level_names <- name_thing("level", n_level)}
  if(!all(c("pi_1_1", "phi_1ij_1", "phi_0ij_1", "n_level") %in% names(init)) |
     any(sapply(init, is.null))
     ){init <- pollinate_ML_ordinal(t_k, n_level = n_level, level_names = level_names)}

  p_t <- init$pi_1_1
  phi_1ij_t <- init$phi_1ij_1
  phi_0ij_t <- init$phi_0ij_1
  y_k <- calc_y_k()

  list_prev <- list()
  list_phi_1ij <- list()
  list_phi_0ij <- list()
  list_A_i <- list()
  list_y_k <- y_k # does not change
  list_g_1 <- list()
  list_g_0 <- list()
  list_q_k1 <- list()
  list_q_k0 <- list()
  list_l_cond <- list()

  for(iter in 1:max_iter){

    A_i <- calc_A_i(phi_1ij_t, phi_0ij_t)
    g_1_t <- calc_g_d(phi_1ij_t)
    g_0_t <- calc_g_d(phi_0ij_t)
    q_k1_t <- calc_q_kd(d = 1)
    q_k0_t <- calc_q_kd(d = 0)
    l_cond_t <- calc_l_cond_ordinal()
    list_prev <- c(list_prev, list(p_t))
    list_phi_1ij <- c(list_phi_1ij, list(phi_1ij_t))
    list_phi_0ij <- c(list_phi_0ij, list(phi_0ij_t))
    list_A_i <- c(list_A_i, list(A_i))
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
          A_i_est = A_i,
          phi_1ij_est = phi_1ij_t,
          phi_0ij_est = phi_0ij_t,
          q_k1_est = q_k1_t),
        names = list(
          method_names = method_names,
          obs_names = obs_names,
          level_names = level_names),
        data = t_k,
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

#' @rdname pollinate_ML
#' @order 3
#' @export
#' @param n_level Used for ordinal methods. Integer number of levels each method contains
#' @param threshold_level Used for ordinal methods. A value from 1 to `n_level` which
#' indicates the initial threshold used to define positive and negative disease states.
#' @param level_names Used for ordinal methods. Optional vector of length `n_level`
#' containing names for each level.
#' @importFrom stats runif

pollinate_ML_ordinal <-
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
        (\(x) x + stats::runif(n_obs, -0.000001, 0.000001))() |> # break ties
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


#' @rdname plot_ML
#' @order 3
#' @returns Ordinal:
#' \item{ROC}{The Receiver Operator Characteristic (ROC) curves estimated for
#' each method}
#' \item{q_k1}{A plot showing how the q values for each observation, k, change when d=1
#' over each iteration of the EM algorithm. Observations can be colored by True
#' state if it is passed (`params$D`).}
#' \item{q_k0}{A plot showing how the q values for each observation, k, change when d=0
#' over each iteration of the EM algorithm. Observations can be colored by True
#' state if it is passed by `params$D`.}
#' \item{q_k1_hist}{A histogram of q_1 values. Observations, k, can be colored by True
#' state if it is passed by `params$D`.}
#' \item{phi_d}{A stacked bar graph representing the estimated CMFs of each
#' method when `d=0` and `d=1`.}
#' @export
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @importFrom stats setNames
#' @importFrom purrr pluck
#' @import tibble
#'
plot_ML_ordinal <-
  function(
    ML_est,
    params = list(
      pi_1_1 = NULL,
      phi_1ij_1 = NULL,
      phi_0ij_1 = NULL,
      D = NULL)){

    # internal functions
    calc_pr_i <- function(b){
      phi_dij <- if(b){ML_est@results$phi_1ij_est}else{ML_est@results$phi_0ij_est}
      sapply(
        1:n_level,
        function(l) colSums(matrix(phi_dij[l:n_level, ], ncol = n_method)) |> pmax(0) |> pmin(1))
    }
    plot_ROC <- function(){

      AUC_data <-
        ML_est@results$A_i_est |>
        as.list() |>
        data.frame() |>
        tidyr::pivot_longer(everything(), names_to = "method", values_to = "value") |>
        dplyr::mutate(label = paste0(method, ": ", sprintf("%0.3f", value))) |>
        dplyr::pull(label) |>
        paste(collapse = "\n")

      AUC_label <- paste0("AUC\n", AUC_data)

      ROC_data |>
        dplyr::bind_rows(expand.grid(method = method_names, level = "", fpr = 0, tpr = 0)) |>
        ggplot2::ggplot(ggplot2::aes(x = fpr, y = tpr, group = method, color = method)) +
        ggplot2::geom_path() +
        ggplot2::geom_point() +
        ggplot2::geom_text(ggplot2::aes(label = level), vjust = "inward", hjust = "inward") + #ggrepel?
        ggplot2::geom_abline(slope = 1, lty = 2, color = "gray50") +
        ggplot2::annotate("text", x = 0.7, y = 0.1, label = AUC_label, hjust = 0, vjust = 0) +
        ggplot2::scale_y_continuous("TPR", limits = c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
        ggplot2::scale_x_continuous("FPR", limits = c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
        ggplot2::scale_color_brewer("", palette = "Set1", na.value = "gray30", drop = FALSE) +
        ggplot2::coord_fixed() +
        ggplot2::theme(panel.background = ggplot2::element_blank(),
                       panel.grid = ggplot2::element_line(color = "gray80"),
                       axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                       legend.position = "bottom") +
        ggplot2::ggtitle("ROC Curves")
    }
    reshape_phi_d <- function(phi_d, d){
      phi_d |>
        as.data.frame() |>
        tibble::rownames_to_column(var = "level") |>
        dplyr::mutate(d = d) |>
        dplyr::mutate(level = factor(level, levels = level_names, ordered = TRUE))
    }
    plot_q_kd <- function(q_kd){
      do.call(rbind, purrr::pluck(ML_est, "prog", q_kd)) |>
        as.data.frame() |>
        stats::setNames(obs_names) |>
        dplyr::mutate(iter = row_number()) |>
        tidyr::pivot_longer(!iter, names_to = "group", values_to = "value") |>
        dplyr::left_join(dis_data, by = "group") |>
        tidyr::replace_na(list(true_D = "Class unknown")) |>
        ggplot2::ggplot(ggplot2::aes(x = iter, y = value, group = group, color = true_D)) +
        ggplot2::geom_line() +
        ggplot2::scale_y_continuous(q_kd, limits = c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
        ggplot2::scale_x_continuous("Iteration", limits = c(0, ML_est@iter)) +
        ggplot2::scale_color_brewer("", palette = "Set1", na.value = "gray30", drop = FALSE) +
        ggplot2::theme(panel.background = ggplot2::element_blank(),
              panel.grid = ggplot2::element_line(color = "gray80"),
              axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
              legend.position = "bottom")
    }

    # dimensions and names
    method_names <- ML_est@names$method_names
    n_method <- length(method_names)
    obs_names <- ML_est@names$obs_names
    n_obs <- length(obs_names)
    level_names <- ML_est@names$level_names
    n_level <- length(level_names)
    if(is.null(params$D)){true_D <- NA}else{true_D <- params$D}

    ### ROC Plot

    fpr_i <- calc_pr_i(FALSE)
    tpr_i <- calc_pr_i(TRUE)

    dimnames(fpr_i) <- list(method_names, level_names)
    dimnames(tpr_i) <- list(method_names, level_names)

    ROC_data <-
      dplyr::left_join(
        fpr_i_long <- as.data.frame(as.table(fpr_i)),
        tpr_i_long <- as.data.frame(as.table(tpr_i)),
        by = c("Var1", "Var2")
      )

    colnames(ROC_data) <- c("method", "level", "fpr", "tpr")

    ROC_plot <- plot_ROC()

    # create progress plots
    # create data frame of disease state and observation name for coloring progress plots
    dis_data <- # disease status
      data.frame(
        true_D = as.character(as.numeric(true_D)),
        group = obs_names) |>
      tidyr::replace_na(list(true_D = "Unknown")) |>
      dplyr::mutate(true_D = paste("Class", true_D))

    q_k1_plot <- plot_q_kd("q_k1")
    q_k0_plot <- plot_q_kd("q_k0")

    # q_k1 histogram

    q_k1_hist <-
      data.frame(q_k1_est = ML_est@results$q_k1_est) |>
      dplyr::mutate(true_D = dis_data$true_D) |>
      ggplot2::ggplot(ggplot2::aes(x = q_k1_est, fill = true_D)) +
      ggplot2::geom_histogram(bins = 40) +
      ggplot2::scale_y_continuous("Observations", limits = c(0, NA), expand = c(0, 0.5)) +
      ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.05), expand = c(0, 0)) +
      ggplot2::scale_fill_brewer("", palette = "Set1", na.value = "gray30", drop = FALSE) +
      ggplot2::theme(panel.background = ggplot2::element_blank(),
            panel.grid = ggplot2::element_line(color = "gray80"),
            panel.spacing = unit(2, "lines"),
            strip.text.y.right = ggplot2::element_text(),
            strip.background = ggplot2::element_rect(fill = "gray80", color = "black"),
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "bottom")

    # phi plot

    phi_d_plot <-
      dplyr::bind_rows(
        reshape_phi_d(ML_est@results$phi_1ij_est, "d = 1"),
        reshape_phi_d(ML_est@results$phi_0ij_est, "d = 0")
      ) |>
        tidyr::pivot_longer(c(-level, -d)) |>
      ggplot2::ggplot(ggplot2::aes(x = name, y = value, fill = level)) +
      ggplot2::geom_col(width = 1, alpha = 0.75, color = "black") +
      ggplot2::scale_y_continuous("p", breaks = seq(0, 1, 0.10), expand = c(0, 0)) +
      ggplot2::scale_x_discrete("Method", expand = c(0, 0)) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::scale_fill_brewer(palette = "Blues") +
      ggplot2::facet_grid(d ~ .) +
      ggplot2::theme(panel.background = ggplot2::element_blank(),
              panel.grid = ggplot2::element_line(color = "gray80"),
              panel.grid.major.x = ggplot2::element_blank(),
              panel.spacing = unit(2, "lines"),
              strip.text.y.right = ggplot2::element_text(),
              strip.background = ggplot2::element_rect(fill = "gray80", color = "black"),
              axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ggplot2::ggtitle("Predicted CMF")

    return(
      list(
        ROC = ROC_plot,
        q_k1 = q_k1_plot,
        q_k0 = q_k0_plot,
        q_k1_hist = q_k1_hist,
        phi_d = phi_d_plot))

  }






