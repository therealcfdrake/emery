# Binary

binary_EM_plotter <- function(data, EM_result, true_pi = NA_real_, true_Se = NA_real_, true_Sp = NA_real_, legend_on = TRUE){

  legend_pos <- if(legend_on){"left"}else{"none"}

  tests <- if(is.null(names(data))){paste0("test", str_pad(1:ncol(data), width = nchar(ncol(data)), side = "left", pad = "0"))}else{names(data)}

  true_pi <-
    ggplot(EM_result, aes(x = Iteration, y = pi)) +
    scale_x_continuous("Iteration", breaks = scales::pretty_breaks()) +
    geom_line() +
    geom_hline(yintercept = true_pi, lty = 2) +
    ggtitle("pi")

  Se_plot <-
    ggplot(EM_result, aes(x = Iteration, y = Se, color = test, group = test)) +
    geom_line() +
    geom_line(data = tibble(test = tests, Se = true_Se, Iteration = max(EM_result$Iteration) * 1.1) %>% bind_rows(filter(EM_result, Iteration == max(Iteration))), lty = 3) +
    geom_point(data = tibble(test = tests, Se = true_Se, Iteration = max(EM_result$Iteration) * 1.1), shape = 1) +
    scale_x_continuous("Iteration", breaks = scales::pretty_breaks()) +
    ggtitle("Se") +
    theme(legend.position = legend_pos)

  Sp_plot <-
    ggplot(EM_result, aes(x = Iteration, y = Sp, color = test, group = test)) +
    geom_line() +
    geom_line(data = tibble(test = tests, Sp = true_Sp, Iteration = max(EM_result$Iteration) * 1.1) %>% bind_rows(filter(EM_result, Iteration == max(Iteration))), lty = 3) +
    geom_point(data = tibble(test = tests, Sp = true_Sp, Iteration = max(EM_result$Iteration) * 1.1), shape = 1) +
    scale_x_continuous("Iteration", breaks = scales::pretty_breaks()) +
    ggtitle("Sp") +
    theme(legend.position = legend_pos)

  Se_span <- c(EM_result$Se, true_Se) %>% range(na.rm = TRUE) %>% diff()
  Se_min <- c(EM_result$Se, true_Se) %>% min(na.rm = TRUE)
  Se_max <- c(EM_result$Se, true_Se) %>% max(na.rm = TRUE)
  Se_center <- c(Se_min, Se_max) %>% mean()
  Sp_span <- c(EM_result$Sp, true_Sp) %>% range(na.rm = TRUE) %>% diff()
  Sp_min <- c(EM_result$Sp, true_Sp) %>% min(na.rm = TRUE)
  Sp_max <- c(EM_result$Sp, true_Sp) %>% max(na.rm = TRUE)
  Sp_center <- c(Sp_min, Sp_max) %>% mean()
  span <- c(Se_span, Sp_span) %>% max()

  limits_data <-
    data.frame(Se = pmin(c(1, 1 - span), c(Se_center + span / 2, Se_center - span / 2)),
               Sp = pmin(c(1, 1 - span), c(Sp_center + span / 2, Sp_center - span / 2)))


  Se_Sp_plot <-
    ggplot(EM_result, aes(x = Sp, y = Se, group = test, color = test)) +
    geom_path() +
    # geom_path(arrow = arrow(length = unit(0.03, "npc"), ends = "last")) +
    geom_point(data = filter(EM_result, Iteration == min(Iteration)), shape = 4, size = 2) +
    geom_point(data = filter(EM_result, Iteration == max(Iteration)), shape = 16, size = 2) +
    geom_point(data = tibble(test = tests, Se = true_Se, Sp = true_Sp), shape = 1, size = 2) +
    geom_line(data = tibble(test = tests, Se = true_Se, Sp = true_Sp) %>% bind_rows(filter(EM_result, Iteration == max(Iteration))), lty = 3) +
    geom_blank(data = limits_data, inherit.aes = FALSE, aes(x = Sp, y = Se)) +
    scale_x_continuous(breaks = seq(0, 1, 0.05)) +
    scale_y_continuous(breaks = seq(0, 1, 0.05)) +
    theme(aspect.ratio = 1) +
    theme(legend.position = legend_pos)

  list(true_pi = true_pi,
       Se_plot = Se_plot,
       Sp_plot = Sp_plot,
       Se_Sp_plot = Se_Sp_plot)

}

# Ordinal


# Continuous

sim_roc_plotter <- function(data){

  res <- 1000
  n_tests <- ncol(data$t_k)

  test_names <- if(is.null(names(data$t_k))){paste0("test_", 1:n_tests)}else{names(data$t_k)}

  test_data <-
    do.call(cbind, data) %>%
    data.frame %>%
    setNames(c(test_names, "D"))

  mins <- apply(data$t_k, MARGIN = 2, min) - 1
  maxs <- apply(data$t_k, MARGIN = 2, max) + 1
  seqs <- lapply(1:n_tests, function(x) seq(mins[x], maxs[x], length.out = res))

  sens <- lapply(1:n_tests, function(x){
    test_sens <- c()
    for(i in 1:res){test_sens[i] <- sum(data$t_k[, x] >= seqs[[x]][i] & data$D == 1) / sum(data$D == 1)}
    test_sens
  }
  )

  spcs <- lapply(1:n_tests, function(x){
    test_spcs <- c()
    for(i in 1:res){test_spcs[i] <- sum(data$t_k[, x] <= seqs[[x]][i] & data$D == 0) / sum(data$D == 0)}
    test_spcs
  }
  )

  plot_data <-
    data.frame(test = lapply(test_names, FUN = rep, res) %>% unlist,
               sen = unlist(sens),
               spc = unlist(spcs))

  plot_roc <- plot_data %>%
    ggplot(aes(x = 1 - spc, y = sen, group = test, color = test)) +
    geom_path() +
    geom_abline(slope = 1, lty = 2) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0))

  a_roc <- plot_data %>%
    group_by(test) %>%
    mutate(dA = (spc - lag(spc, 1)) * (sen + lag(sen, 1)) / 2) %>%
    summarize(A = sum(dA, na.rm = TRUE))

  list(plot = plot_roc, area = a_roc)

}
