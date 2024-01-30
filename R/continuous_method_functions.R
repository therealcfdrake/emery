
# d disease status
# i methods
# j response
# k observations

### Build Sample Data Set

#' @rdname generate_multimethod_data
#' @order 4
#' @export
#' @param mu_i1,mu_i0 Used for continuous methods. Vectors of length n_method of the method mean values for positive (negative) observations
#' @param sigma_i1,sigma_i0 Used for continuous methods. Covariance matrices of method positive (negative) observations
#'
#' @importFrom mvtnorm rmvnorm
#'

generate_multimethod_continuous <-
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
    obs_names = NULL,
    n_method_subset = n_method,
    first_reads_all = FALSE
    ){

  if(is.null(method_names)){method_names <- name_thing("method", n_method)}
  if(is.null(obs_names)){obs_names <- name_thing("obs", n_obs)}

  dis <- define_disease_state(D, n_obs, prev)

  X <- mvtnorm::rmvnorm(n = dis$pos, mean = mu_i1, sigma = sigma_i1)
  Y <- mvtnorm::rmvnorm(n = dis$neg, mean = mu_i0, sigma = sigma_i0)

  subset_matrix <-
    censor_data(
      n_obs = dis$n_obs,
      first_reads_all = first_reads_all,
      n_method_subset = n_method_subset,
      n_method = n_method)

  generated_data <- rbind(X, Y) * subset_matrix

  dimnames(generated_data) <- list(obs_names, method_names)
  names(mu_i1) <- method_names
  dimnames(sigma_i1) <- list(method_names, method_names)
  names(mu_i0) <- method_names
  dimnames(sigma_i0) <- list(method_names, method_names)

  params <-
    list(
      n_method = n_method,
      n_obs = dis$n_obs,
      prev = dis$prev,
      D = stats::setNames(dis$D, obs_names),
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


#' @rdname estimate_ML
#' @order 3
#' @export
#' @importFrom stats setNames pnorm
#'
estimate_ML_continuous <-
  function(data,
           init = list(
             prev_1 = NULL,
             mu_i1_1 = NULL,
             sigma_i1_1 = NULL,
             mu_i0_1 = NULL,
             sigma_i0_1 = NULL),
           max_iter = 100,
           tol = 1e-7,
           save_progress = TRUE){

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
      sigma_d <- matrix(nrow = n_method, ncol = n_method, dimnames = list(method_names, method_names))
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
    stats::setNames(stats::pnorm(eta_j_m), method_names)
  }

  t_k <- as.matrix(data)
  n_method <- ncol(t_k)
  n_obs <- nrow(t_k)
  method_names <- if(is.null(colnames(t_k))){name_thing("method", n_method)}else{colnames(t_k)}
  obs_names <- if(is.null(rownames(t_k))){name_thing("obs", n_obs)}else{rownames(t_k)}

  dimnames(t_k) <- list(obs_names, method_names)

  if(!all(c("prev_1", "mu_i1_1", "sigma_i1_1", "mu_i0_1", "sigma_i0_1") %in% names(init)) |
     any(sapply(init, is.null))
  ){init <- pollinate_ML_continuous(t_k)}

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
        names = list(
          method_names = method_names,
          obs_names = obs_names),
        data = t_k,
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

#' @rdname pollinate_ML
#' @order 4
#' @export
#' @param prev A double between 0-1 representing the proportion of positives in the population
#' @param q_seeds Used for continuous methods. A vector of length 2 representing the quantiles at which the two groups are assumed to be centered
#' @param high_pos Used for continuous methods. A logical indicating whether larger values are considered "positive"
#' @importFrom stats quantile

pollinate_ML_continuous <-
  function(data,
           prev = 0.5,
           q_seeds = c((1 - prev) / 2, 1 - (prev / 2)),
           high_pos = TRUE){

  method_names <- if(is.null(colnames(data))){name_thing("method", ncol(data))}else{colnames(data)}
  obs_names <- if(is.null(rownames(data))){name_thing("obs", nrow(data))}else{rownames(data)}

  # adjust seeds depending on whether high (default) or low values are associated with "positive" diagnosis
  q_seeds <- sort(q_seeds, decreasing = high_pos)

  muD_mat <- apply(data, MARGIN = 2, stats::quantile, q_seeds)

  mu_i1_1 <- unlist(muD_mat[1, ])
  mu_i0_1 <- unlist(muD_mat[2, ])

  sigma_i1_1 <- diag(mu_i1_1 / 2)
  sigma_i0_1 <- diag(mu_i0_1 / 2)

  names(mu_i1_1) <- method_names
  dimnames(sigma_i1_1) <- list(method_names, method_names)
  names(mu_i0_1) <- method_names
  dimnames(sigma_i0_1) <- list(method_names, method_names)

  return(
    list(
      prev_1 = prev,
      mu_i1_1 = mu_i1_1,
      sigma_i1_1 = sigma_i1_1,
      mu_i0_1 = mu_i0_1,
      sigma_i0_1 = sigma_i0_1)
  )

}



#' @rdname plot_ML
#' @order 4
#' @returns Continuous:
#' \item{ROC}{The Receiver Operator Characteristic (ROC) curves estimated for
#' each method}
#' \item{z_k1}{A plot showing how the z_k1 values for each observation change
#' over each iteration of the EM algorithm. Observations can be colored by True
#' state if it is passed (`params$D`).}
#' \item{z_k0}{A plot showing how the z_k0 values for each observation change
#' over each iteration of the EM algorithm. Observations can be colored by True
#' state if it is passed (`params$D`).}
#' \item{z_k1_hist}{A histogram of z_k1 values. Observations can be colored by True
#' state if it is passed (`params$D`).}
#' @export
#' @import dplyr
#' @import ggplot2
#' @importFrom stats setNames
#' @importFrom purrr pluck
plot_ML_continuous <-
  function(
    ML_est,
    params = list(
      prev_1 = NULL,
      mu_i1_1 = NULL,
      sigma_i1_1 = NULL,
      mu_i0_1 = NULL,
      sigma_i0_1 = NULL,
      D = NULL)){

    # internal functions
    calc_pr_i <- function(b){ # (T/F) Positive Rate calculator
      z_kd <- if(b){ML_est@results$z_k1_est}else{ML_est@results$z_k0_est}
      apply(ML_est@data, 2, function(column){
        sapply(column, function(row, column){
          sum(z_kd[which(column >= row)])
        }, column = column)}) / sum(z_kd)
    }
    plot_ROC <- function(){

      AUC_data <-
        ML_est@results$A_j_est |>
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
        ggplot2::geom_abline(slope = 1, lty = 2, color = "gray50") +
        ggplot2::annotate("text", x = 0.6, y = 0.1, label = AUC_label, hjust = 0, vjust = 0) +
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
    plot_z_kd <- function(z_kd){
      do.call(rbind, purrr::pluck(ML_est, "prog", z_kd)) |>
        as.data.frame() |>
        stats::setNames(obs_names) |>
        dplyr::mutate(iter = row_number()) |>
        tidyr::pivot_longer(!iter, names_to = "group", values_to = "value") |>
        dplyr::left_join(dis_data, by = "group") |>
        tidyr::replace_na(list(true_D = "Class unknown")) |>
        ggplot2::ggplot(ggplot2::aes(x = iter, y = value, group = group, color = true_D)) +
        ggplot2::geom_line() +
        ggplot2::scale_y_continuous(z_kd, limits = c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
        ggplot2::scale_x_continuous("Iteration", limits = c(0, ML_est@iter)) +
        ggplot2::scale_color_brewer("", palette = "Set1", na.value = "gray30", drop = FALSE) +
        ggplot2::theme(panel.background = ggplot2::element_blank(),
                       panel.grid = ggplot2::element_line(color = "gray80"),
                       axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                       legend.position = "bottom")
    }

    method_names <- ML_est@names$method_names
    n_method <- length(method_names)
    obs_names <- ML_est@names$obs_names
    n_obs <- length(obs_names)
    if(is.null(params$D)){true_D <- NA}else{true_D <- params$D}

    fpr_i <- calc_pr_i(FALSE)
    tpr_i <- calc_pr_i(TRUE)

    ROC_data <-
      dplyr::left_join(
        fpr_i_long <- as.data.frame(as.table(fpr_i)),
        tpr_i_long <- as.data.frame(as.table(tpr_i)),
        by = c("Var1", "Var2")
      ) |> stats::setNames(c("obs", "method", "fpr", "tpr")) |>
      dplyr::arrange(method, desc(fpr), desc(tpr))

    ROC_plot <- plot_ROC()

    # create progress plots
    # create data frame of disease state and observation name for coloring progress plots
    dis_data <- # disease status
      data.frame(
        true_D = as.character(as.numeric(true_D)),
        group = obs_names) |>
      tidyr::replace_na(list(true_D = "Unknown")) |>
      dplyr::mutate(true_D = paste("Class", true_D))

    z_k1_plot <- plot_z_kd("z_k1")
    z_k0_plot <- plot_z_kd("z_k0")

    z_k1_hist <-
      data.frame(z_k1_est = ML_est@results$z_k1_est) |>
      dplyr::mutate(true_D = dis_data$true_D) |>
      ggplot2::ggplot(ggplot2::aes(x = z_k1_est, fill = true_D)) +
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

    return(
      list(
        ROC = ROC_plot,
        z_k1 = z_k1_plot,
        z_k0 = z_k0_plot,
        z_k1_hist = z_k1_hist))

    }
