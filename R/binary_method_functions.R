#' @rdname generate_multimethod_data
#' @order 2
#' @export
#' @param se,sp Used for binary methods. A vector of length n_method of
#' values between 0-1 representing the sensitivity and specificity of the methods.
#' @param n_method_subset Used for binary methods. An integer defining how
#' many methods to select at random to produce a result for each observation
#' @param first_reads_all Used for binary methods. A logical which forces method
#'  1 to have a result for every observation
#' @importFrom stats rbinom setNames

generate_multimethod_binary <-
  function(n_method = 3,
           n_obs = 100,
           prev = 0.5,
           D = NULL,
           se = rep(0.9, n_method),
           sp = rep(0.9, n_method),
           method_names = NULL,
           obs_names = NULL,
           n_method_subset = n_method,
           first_reads_all = FALSE
  ){

    # Create unique names for each method and observation if names are not provided
    if(is.null(method_names)){method_names <- name_thing("method", n_method)}
    if(is.null(obs_names)){obs_names <- name_thing("obs", n_obs)}

    # Define the True disease state based on input criteria
    dis <- define_disease_state(D, n_obs, prev)

    # Create dataframe to (optionally) randomly censor observations such that
    # only "n_method_subset" results are reported for each observation
      # lapply(1:dis$n_obs, function(i)
      #   if(!first_reads_all){sample(c(rep(1, n_method_subset), rep(NA, n_method - n_method_subset)), n_method, replace = FALSE)
      #   }else{
      #     c(1, sample(c(rep(1, n_method_subset - 1), rep(NA, n_method - n_method_subset)), n_method - 1, replace = FALSE))
      #   }
      # ) |>
      # do.call(what = rbind, args = _)

    subset_matrix <-
      censor_data(
        n_obs = dis$n_obs,
        first_reads_all = first_reads_all,
        n_method_subset = n_method_subset,
        n_method = n_method)

    # Build simulated data set based on input criteria
    generated_data <-
      lapply(1:n_method, function(i) stats::rbinom(dis$n_obs, 1, se[i] * dis$D + (1 - sp[i]) * (1 - dis$D))) |>
      do.call(what = cbind, args = _) * subset_matrix

    dimnames(generated_data) <- list(obs_names, method_names)

    # Calculate individual method se based on input criteria. This will differ
    # slightly from "se" due to random sampling.
    se_observed <-
      cbind(generated_data, D = dis$D) |>
      (\(x) x[which(x[, "D"] == 1), ])() |>
      subset(select = -D) |>
      colMeans(na.rm = TRUE)

    # Calculate individual method sp based on input criteria. This will differ
    # slightly from "sp" due to random sampling.
    sp_observed <-
      cbind(generated_data, D = dis$D) |>
      (\(x) x[which(x[, "D"] == 0), ])() |>
      subset(select = -D) |>
      colMeans(na.rm = TRUE) |>
      (\(y) 1 - y)()

    params <-
      list(
        n_method = n_method,
        n_obs = dis$n_obs,
        prev = dis$prev,
        se = se,
        sp = sp,
        D = stats::setNames(dis$D, obs_names),
        se_observed = se_observed,
        sp_observed = sp_observed,
        method_names = method_names,
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
#' @order 2
#' @export

estimate_ML_binary <-
  function(data,
           init = list(
             prev_1 = NULL,
             se_1 = NULL,
             sp_1 = NULL),
           max_iter = 100,
           tol = 1e-7,
           save_progress = TRUE){

  calc_A2 <- function(){
      data |>
        apply(1, FUN = function(y) (se_m ^ y) * ((1 - se_m) ^ (1 - y))) |>
        apply(2, FUN = function(z) prod(z, na.rm = TRUE)) |>
        (\(x) x * prev_m)()
    }
  calc_B2 <- function(){
      data |>
        apply(1, FUN = function(y) ((1 - sp_m) ^ y) * (sp_m ^ (1 - y))) |>
        apply(2, FUN = function(z) prod(z, na.rm = TRUE)) |>
        (\(x) x * (1 - prev_m))()

    }
  calc_qk <- function(A2, B2){
      A2 / (A2 + B2)
    }
  calc_next_se <- function(){
      data_mat <- as.matrix(data)
      data_mat[is.na(data_mat)] <- 0 # added to handle NA, used to calc appropriate denominators
      (qk_m %*% data_mat) / (qk_m %*% !is.na(as.matrix(data)))
      # (qk_m %*% as.matrix(data)) / sum(qk_m) # original function before changes to accept NA
    }
  calc_next_sp <- function(){
      data_mat <- 1 - as.matrix(data)
      data_mat[is.na(data_mat)] <- 0 # added to handle NA, used to calc appropriate denominators
      ((1 - qk_m) %*% data_mat) / ((1 - qk_m) %*% !is.na(as.matrix(data)))
      # ((1 - qk_m) %*% (1 - as.matrix(data))) / sum(1 - qk_m) # original function before changes to accept NA
    }
  calc_next_prev <- function(){
      mean(qk_m)
    }

  if(!all(c("prev_1", "se_1", "sp_1") %in% names(init)) | any(sapply(init, is.null))){init <- pollinate_ML_binary(data)}

  method_names <- if(is.null(colnames(data))){name_thing("method", ncol(data))}else{colnames(data)}
  obs_names <- if(is.null(rownames(data))){name_thing("obs", nrow(data))}else{rownames(data)}

  dimnames(data) <- list(obs_names, method_names)

  # starting values
  se_m <- init$se_1
  sp_m <- init$sp_1
  prev_m <- init$prev_1

  # initialize lists
  list_se <- list()
  list_sp <- list()
  list_prev <- list()
  list_A2 <- list()
  list_B2 <- list()
  list_qk <- list()

  # iterate

  for(iter in 1:(max_iter)){

    A2_m <- calc_A2()
    B2_m <- calc_B2()
    qk_m <- calc_qk(A2_m, B2_m)

    list_se <- c(list_se, list(se_m))
    list_sp <- c(list_sp, list(sp_m))
    list_prev <- c(list_prev, list(prev_m))
    list_A2 <- c(list_A2, list(A2_m))
    list_B2 <- c(list_B2, list(B2_m))
    list_qk <- c(list_qk, list(qk_m))

    if(iter > 1){
      if(
        max(abs(list_se[[iter]] - list_se[[iter - 1]]),
            abs(list_sp[[iter]] - list_sp[[iter - 1]]),
            abs(list_prev[[iter]] - list_prev[[iter - 1]])) < tol){
        break}
      }

    se_m <- calc_next_se()
    sp_m <- calc_next_sp()
    prev_m <- calc_next_prev()

  }

  prev_prog <- do.call(rbind, list_prev)
  se_prog <- do.call(rbind, list_se)
  sp_prog <- do.call(rbind, list_sp)
  A2_prog <- do.call(rbind, list_A2)
  B2_prog <- do.call(rbind, list_B2)
  qk_prog <- do.call(rbind, list_qk)

  output <-
    new("MultiMethodMLEstimate",
        results = list(
          prev_est = unlist(prev_m),
          se_est = unlist(se_m),
          sp_est = unlist(sp_m),
          qk_est = unlist(qk_m)),
        names = list(
          method_names = method_names,
          obs_names = obs_names),
        data = data,
        iter = iter,
        type = "binary")

  if(save_progress){
    output@prog <-
      list(
        prev = prev_prog,
        se = se_prog,
        sp = sp_prog,
        A2 = A2_prog,
        B2 = B2_prog,
        qk = qk_prog
        )
  }

  return(output)

}

#' @rdname pollinate_ML
#' @order 2
#' @export
#' @importFrom stats setNames runif

pollinate_ML_binary <-
  function(
    data,
    ...
    ){

  method_names <- if(is.null(colnames(data))){name_thing("method", ncol(data))}else{colnames(data)}

  n_obs <- nrow(data)

  D_majority <- data |>
    rowMeans(na.rm = TRUE) |>
    (\(x) x + stats::runif(n_obs, -0.000001, 0.000001))() |> # randomly break ties
    round()

  # Estimate initial prevalence wrt majority classification
  prev_1 <- mean(D_majority, na.rm = TRUE) |> stats::setNames("prev")

  # Estimate individual method se wrt majority classification
  se_1 <- data[D_majority == 1, ] |> colMeans(na.rm = TRUE) |> stats::setNames(method_names)

  # Estimate individual method sp wrt majority classification
  sp_1 <- 1 - (data[D_majority == 0, ] |> colMeans(na.rm = TRUE) |> stats::setNames(method_names))

  return(
    list(prev_1 = prev_1,
         se_1 = se_1,
         sp_1 = sp_1)
  )

}

#' @rdname plot_ML
#' @order 2
#' @returns Binary:
#' \item{prev}{A plot showing how the prevalence estimate changes with each
#' iteration of the EM algorithm}
#' \item{se}{A plot showing how the sensitivity estimates of each method change with each
#' iteration of the EM algorithm}
#' \item{sp}{A plot showing how the specificity estimates of each method change with each
#' iteration of the EM algorithm}
#' \item{qk}{A plot showing how the q values for each observation k change
#' over each iteration of the EM algorithm}
#' \item{qk_hist}{A histogram of q values. Observations, k, can be colored by True
#' state if it is passed by `params$D`.}
#' \item{se_sp}{A plot showing the path the sensitivity and specificity estimates
#' for each method follows during the EM algorithm. True sensitivity and specificity
#' values can be passed by `params$se` and `params$sp`, respectively. This is useful
#' for comparing algorithm results when applied to simulation data where True parameter
#' values are known.}
#' @export
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @importFrom purrr pluck
#' @importFrom stats setNames
#' @import tibble

plot_ML_binary <-
  function(
    ML_est,
    params = list(
      prev = NULL,
      se = NULL,
      sp = NULL,
      D = NULL)){

  prog_plots <- c("prev", "se", "sp", "qk") # A2 and B2 not plotted

  method_names <- ML_est@names$method_names
  n_method <- length(method_names)
  obs_names <- ML_est@names$obs_names
  n_obs <- length(obs_names)
  if(is.null(params$D)){true_D <- NA}else{true_D <- params$D}

  ### Se/Sp scatter plot

  # create long data frame of sequence of se/sp estimates
    se_sp_data <-
      dplyr::left_join(
        as.data.frame(ML_est@prog$se) |> dplyr::mutate(iter = row_number()) |> tidyr::pivot_longer(!iter, names_to = "method", values_to = "se"),
        as.data.frame(ML_est@prog$sp) |> dplyr::mutate(iter = row_number()) |> tidyr::pivot_longer(!iter, names_to = "method", values_to = "sp"),
        by = join_by(iter, method)
      )

  # create data frame containing final estimate and true value, if known. This is for assessing accuracy of estimates in simulations when true value is known.
    se_sp_result <-
      data.frame(
        method = rep(method_names, 2),
        se = c(ML_est@results$se_est, params$se) |> (\(x) rep(x, 2 * n_method / length(x)))(),
        sp = c(ML_est@results$sp_est, params$sp) |> (\(x) rep(x, 2 * n_method / length(x)))(),
        shape = c(
          rep("final", n_method),
          rep("truth", n_method)
        )
        )

  # create se/sp line plot showing change in estimates over time
    se_sp_plot <-
      se_sp_data |>
      ggplot2::ggplot(ggplot2::aes(x = sp, y = se, group = method, color = method)) +
      ggplot2::geom_path() +
      ggplot2::geom_point(data = se_sp_result, ggplot2::aes(shape = shape)) +
      ggplot2::scale_shape_manual(values = c(16, 1)) +
      ggplot2::geom_line(data = se_sp_result, lty = 2) +
      ggplot2::coord_fixed() +
      ggplot2::scale_y_continuous("Sensitivity", limits = c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
      ggplot2::scale_x_continuous("Specificity", limits = c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
      ggplot2::scale_color_brewer(palette = "Set1") +
      ggplot2::theme(panel.background = element_blank(),
              panel.grid = element_line(color = "gray80"),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



    ### ML progress plots

  # create data frame of disease state and observation name for coloring progress plots
    dis_data <- # disease status
      data.frame(
        true_D = as.character(true_D),
        group = obs_names) |>
      tidyr::replace_na(list(true_D = "Unknown")) |>
      dplyr::mutate(true_D = paste("Class", true_D))

  # create progress plots
    list_prog_plots <-
      lapply(prog_plots, function(j){
        df_plot <- as.data.frame(purrr::pluck(ML_est, "prog", j)) |>
          dplyr::mutate(iter = row_number())
        df_plot |>
          tidyr::pivot_longer(!iter, names_to = "group", values_to = "value") |>
          dplyr::left_join(dis_data, by = "group") |>
          dplyr::mutate(color_col = dplyr::case_when( # define which column to color lines by
            j %in% c("A2", "B2", "qk") ~ true_D,
            j %in% c("se", "sp") ~ group,
            j %in% c("prev") ~ "Estimate")) |>
          ggplot2::ggplot(ggplot2::aes(x = iter, y = value, group = group, color = color_col)) +
          ggplot2::geom_line() +
          ggplot2::scale_y_continuous(j, limits = c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
          ggplot2::scale_x_continuous("Iteration", limits = c(0, ML_est@iter)) +
          ggplot2::scale_color_brewer("", palette = "Set1", na.value = "gray30", drop = FALSE) +
          ggplot2::theme(panel.background = element_blank(),
                panel.grid = element_line(color = "gray80"),
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                legend.position = "bottom")
      })

    qk_hist <-
      data.frame(qk_est = ML_est@results$qk_est) |>
      dplyr::mutate(true_D = dis_data$true_D) |>
      ggplot2::ggplot(ggplot2::aes(x = qk_est, fill = true_D)) +
      ggplot2::geom_histogram(bins = 40, boundary = -0.025) +
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

    ### Output
    return(
      c(
        list_prog_plots,
        list(qk_hist),
        list(se_sp_plot)
        ) |> stats::setNames(c("prev", "se", "sp", "qk", "qk_hist", "se_sp"))
    )

  }

