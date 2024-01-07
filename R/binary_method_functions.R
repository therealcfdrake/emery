## Binary

gen_multi_bin <-
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

    # "True" values have no suffix and are provided by user (population parameters)
    # "*_observed" - calculated using observed values with the "True" value as the reference. These show how closely the random sample matches the desired input.

    # Create unique names for each method and observation if names are not provided
    if(is.null(method_names)){method_names <- name_thing("method", n_method)}
    if(is.null(obs_names)){obs_names <- name_thing("obs", n_obs)}

    # Define the True disease state based on input criteria
    dis <- define_disease_state(D, n_obs, prev)

    # Create dataframe to (optionally) randomly censor observations such that only "n_method_subset" results are reported for each observation
    subset_matrix <-
      lapply(1:dis$n_obs, function(i)
        if(!first_reads_all){sample(c(rep(1, n_method_subset), rep(NA, n_method - n_method_subset)), n_method, replace = FALSE)
        }else{
          c(1, sample(c(rep(1, n_method_subset - 1), rep(NA, n_method - n_method_subset)), n_method - 1, replace = FALSE))
        }
      ) |>
      do.call(what = rbind, args = _)

    # Build simulated data set based on input criteria
    generated_data <-
      lapply(1:n_method, function(i) rbinom(dis$n_obs, 1, se[i] * dis$D + (1 - sp[i]) * (1 - dis$D))) |>
      do.call(what = cbind, args = _) * subset_matrix

    dimnames(generated_data) <- list(obs_names, method_names)

    # Calculate individual method se based on input criteria. This will differ slightly from "se" due to random sampling.
    se_observed <-
      cbind(generated_data, D = dis$D) |>
      (\(x) x[which(x[, "D"] == 1), ])() |>
      subset(select = -D) |>
      colMeans(na.rm = TRUE)

    # Calculate individual method sp based on input criteria. This will differ slightly from "sp" due to random sampling.
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
        D = setNames(dis$D, obs_names),
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

estimate_ML_binary <-
  function(data,
           init = list(prev_1 = NULL, se_1 = NULL, sp_1 = NULL),
           max_iter = 100,
           tol = 1e-7,
           save_progress = FALSE){

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

  if(any(sapply(init, is.null))){init <- pollinate_ML_binary(data)}

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

  converged <- FALSE
  iter <- 1

  # iterate

  for(iter in 1:max_iter){

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

    iter <- iter + 1

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
          sp_est = unlist(sp_m)),
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

pollinate_ML_binary <-
  function(data){

  method_names <- if(is.null(names(data))){name_thing("method", ncol(data))}else{names(data)}

  n_obs <- nrow(data)

  D_majority <- data |>
    rowMeans(na.rm = TRUE) |>
    (\(x) x + runif(n_obs, -0.000001, 0.000001))() |> # randomly break ties
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

#' Title
#'
#' @param ML_est
#' @param plots
#'
#' @return A list of plots.
#' @export
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @importFrom purrr pluck
#' @import tibble
#'
#' @examples
#'

plot_ML_binary <-
  function(
    ML_est,
    params = list(prev = NULL, se = NULL, sp = NULL, D = NULL)){

  prog_plots <- c("prev", "se", "sp", "A2", "B2", "qk")

  n_method <- length(ML_est@results$se_est)
    se_sp_data <-
      left_join(
        as.data.frame(ML_est@prog$se) %>% mutate(iter = row_number()) %>% pivot_longer(!iter, names_to = "method", values_to = "se"),
        as.data.frame(ML_est@prog$sp) %>% mutate(iter = row_number()) %>% pivot_longer(!iter, names_to = "method", values_to = "sp"),
        by = join_by(iter, method)
      )

    ### Se/Sp scatter plot
    se_sp_result <-
      data.frame(
        method = rep(colnames(ML_est@results$se_est), 2),
        se = c(ML_est@results$se_est, params$se),
        sp = c(ML_est@results$sp_est, params$sp),
        shape = c(
          rep("final", n_method),
          rep("param", n_method)
        )
        )

    se_sp_plot <-
      se_sp_data %>%
        ggplot(aes(x = sp, y = se, group = method, color = method)) +
        geom_line() +
        geom_point(data = se_sp_result, aes(shape = shape)) +
        scale_shape_manual(values = c(16, 1)) +
        geom_line(data = se_sp_result, lty = 2) +
        coord_fixed() +
        scale_y_continuous("Sensitivity", limits = c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
        scale_x_continuous("Specificity", limits = c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
        scale_color_brewer(palette = "Set1") +
        theme(panel.background = element_blank(),
              panel.grid = element_line(color = "gray80"),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

    ### ML progress plots
    dis_data <- # disease status
      data.frame(
        D = as.character(as.numeric(params$D)),
        group = as.character(names(params$D)))

    list_prog_plots <-
      lapply(prog_plots, function(j){
        df_plot <- as.data.frame(pluck(ML_est, "prog", j)) %>%
          mutate(iter = row_number())
        color_col <- case_when( # define which column to color lines by
          j %in% c("A2", "B2", "qk") ~ "D",
          TRUE ~ "group")
        df_plot %>%
          pivot_longer(!iter, names_to = "group", values_to = "value") %>%
          left_join(dis_data, by = "group") %>%
          ggplot(aes(x = iter, y = value, group = group, color = .data[[color_col]])) +
          geom_line() +
          scale_y_continuous(j, limits = c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
          scale_x_continuous("Iteration", limits = c(0, ML_est@iter)) +
          scale_color_brewer(palette = "Set1") +
          theme(panel.background = element_blank(),
                panel.grid = element_line(color = "gray80"),
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                legend.position = "bottom")
      })

    ### Output
    return(
      c(
        list_prog_plots,
        list(se_sp_plot)
        )
    )

  }

plot_ML_binary_results()

a <- generate_multimethod_data(5, 100, 0.2, type = "binary", se = runif(5, 0.7, 0.9), sp = runif(5, 0.7, 0.9))
b <- estimate_ML_binary(a$generated_data, save_progress = TRUE)
# plot_ML_binary(b, plots = "all")
left_join(
  as.data.frame(b@prog$se) %>% mutate(iter = row_number()) %>% pivot_longer(!iter, names_to = "method", values_to = "se"),
  as.data.frame(b@prog$sp) %>% mutate(iter = row_number()) %>% pivot_longer(!iter, names_to = "method", values_to = "sp")
)
