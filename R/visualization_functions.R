

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
