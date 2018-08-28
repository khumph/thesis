maxPlots <- function(Q, dat_long, ex_ID, mon, n = 4, int, noise_pred, seed = 1) {
  dat <- filter(dat_long, month == mon)
  set.seed(seed)
  ids <- dat %>% filter(!is.na(M_next)) %>%
    distinct(ID) %>% sample_n(n) %>% flatten_dbl()
  dat <- dat %>% filter(ID %in% ids)
  
  df <- max_df(dat, model = Q$mod_list[[mon + 1]],
               truth = F, pred = F, nested = T) %>% ungroup() %>%
    mutate(ID = factor(ID))
  
  df_best <- dat %>% mutate(reward = Qhat) %>% 
    max_df(model = NULL, truth = T, pred = F, nested = T) %>%
    mutate(ID = factor(ID))
  
  ggplot(data = df) +
    geom_line(mapping = aes(x = dose, y = preds, color = ID)) +
    geom_line(
      data = df_best,
      aes(x = dose, y = preds, color = ID),
      linetype = 2
    ) +
    labs(caption = "dotted lines are true values")
}

plots_tab <- function(dat_test ) {
  dat_long_summ <- dat_test %>%
    # filter(!(group  %in% c(0.2, 0.4, 0.6, 0.8, 1))) %>% 
    group_by(rep, ID) %>%
    mutate(
      cumSurv = prod(1 - pdeath[1:6]),
      cured = sum(tumor_mass == 0, na.rm = T) > 0
    ) %>% group_by(group, month) %>% 
    summarise(
      median_dose = median(dose, na.rm = T),
      mean_tox = mean(toxicity, na.rm = T),
      mean_tumor = mean(tumor_mass, na.rm = T),
      median_reward = median(tot_reward, na.rm = T),
      mean_cumSurv = mean(cumSurv, na.rm = T),
      prop_cured = mean(cured, na.rm = T)
    ) %>% mutate(sum_means = mean_tox + mean_tumor)
  
  plot_dose <- ggplot(data = dat_long_summ) +
    geom_line(mapping = aes(
      x = month,
      y = median_dose,
      color = group,
      group = group
    ))
  
  plot_tox <- ggplot(data = dat_long_summ) +
    geom_line(mapping = aes(
      x = month,
      y = mean_tox,
      color = group,
      group = group
    ))
  
  plot_tumor <- ggplot(data = dat_long_summ) +
    geom_line(mapping = aes(
      x = month,
      y = mean_tumor,
      color = group,
      group = group
    ))
  
  plot_sum <- ggplot(data = dat_long_summ) +
    geom_line(mapping = aes(
      x = month,
      y = sum_means,
      color = group,
      group = group
    ))
  
  tab_deaths <- dat_long_summ %>%
    group_by(group) %>%
    summarise(mean_cumSurv = mean(mean_cumSurv)) %>%
    arrange(desc(mean_cumSurv))
  
  tab_cured <- dat_long_summ %>%
    group_by(group) %>%
    summarise(mean_prop_cured = mean(prop_cured)) %>%
    arrange(desc(mean_prop_cured))
  
  list(
    plot_dose = plot_dose,
    plot_tox = plot_tox,
    plot_tumor = plot_tumor,
    plot_sum = plot_sum,
    table_deaths = tab_deaths,
    table_cured = tab_cured
  )
}


reward_plot <- function(dat, Q_list) {
  plot_reward <- ggplot(data = dat) +
    geom_boxplot(mapping = aes(
      x = group,
      y = exp(tot_reward),
      color = group
    ), notch = T) +
    geom_point(mapping = aes(
      x = group,
      y = exp(mean(tot_reward))
    )) +
    labs(y = "Months of survival",
         x = "",
         title = "Survival times averaged over 25 test sets") + 
    theme(legend.position = "none") + coord_trans(y = "log")
  
  tab_reward <- dat %>%
    group_by(group) %>%
    summarise(median_tot_reward = median(tot_reward, na.rm = T),
              mean_reward = mean(tot_reward)) %>%
    arrange(desc(median_tot_reward))
  
  tab_MSE <- dat %>%
    filter(group %in% names(Q_list)) %>% group_by(group, month) %>% 
    summarise(MSE_dose = mean((dose - best_dose) ^ 2, na.rm = T),
              MAE_dose = mean(abs(dose - best_dose), na.rm = T))
  
  list(
    plot_reward = plot_reward,
    table_rewards = tab_reward,
    tab_MSE = tab_MSE
  )
}

results_plot <- function(dat) {
  ggplot(data = dat) +
    geom_boxplot(mapping = aes(
      x = group,
      y = exp(tot_reward) #,
      # color = group
    ), notch = T) +
    labs(y = "Months of survival",
         x = "Treatment regime") + #,
    # title = "Survival times averaged over 25 test sets") + 
    theme(legend.position = "none") + coord_trans(y = "log")
}

save_plot <- function(dat, int, noise, noise_pred) {
  type_nm <- paste_type_nm(int, noise, noise_pred)
  ggsave(paste0("results-plot-", type_nm, ".png"),
         plot = results_plot(dat),
         width = 7, height = 4)
}
