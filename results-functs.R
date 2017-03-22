# results functions -------------------------------------------------------

indPlot <- function(data, ex_ID) {
  tox_mass_plot <- ggplot(
    data = filter(data, ID == ex_ID)
  ) +
    geom_line(
      mapping = aes(x = month, y = toxicity, group = ID), color = "green"
    ) +
    geom_line(
      mapping = aes(x = month, y = tumor_mass, group = ID) 
    )
  
  dose_plot <- ggplot(
    data = filter(data, ID == ex_ID)
  ) +
    geom_line(
      mapping = aes(x = month, y = dose, group = ID), color = "red"
    ) + ylim(0, 1)
  
  list(grid.arrange(tox_mass_plot, dose_plot, nrow = 2, ncol = 1), ex_ID)
}

maxPlots <- function(Q, ex_ID, mon = 5) {
  dat <- filter(Q$data, month == mon)
  nested_df <-
    max_df(dat,
           Q$mod_list[[mon + 1]],
           Q$formula,
           mod_type = Q$mod_type,
           nested = T)
  
  onePlot <-
    ggplot(filter(nested_df, ID == ex_ID), aes(x = dose, y = preds)) +
    geom_point()
  
  ids <- sample(1:1000, 50)
  manyPlot <- ggplot(data = filter(nested_df, ID %in% ids)) +
    geom_line(mapping = aes(x = dose, y = preds, group = ID))
  
  list(onePlot, manyPlot, ex_ID, ids)
}

plots_tab <- function(dat_test_long) {
  dat_long_summ <- dat_test_long %>% group_by(ID) %>%
    mutate(
      cumSurv = prod(1 - pdeath[1:6]),
      cured = sum(tumor_mass == 0, na.rm = T) > 0
    ) %>% group_by(group, month) %>% 
    summarise(
      mean_tox = mean(toxicity, na.rm = T),
      mean_tumor = mean(tumor_mass, na.rm = T),
      mean_reward = mean(tot_reward, na.rm = T),
      mean_cumSurv = mean(cumSurv, na.rm = T),
      prop_cured = mean(cured, na.rm = T)
    ) %>% mutate(sum_means = mean_tox + mean_tumor)
  
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
  
  tab_reward <- dat_long_summ %>%
    group_by(group) %>%
    summarise(mean_tot_reward = mean(mean_reward, na.rm = T)) %>%
    arrange(desc(mean_tot_reward))
  
  plot_reward <- ggplot(data = tab_reward) +
    geom_bar(mapping = aes(
      x = group,
      y = exp(mean_tot_reward),
      fill = group
    ), stat = "identity") + 
    labs(y = "average months of survival")
  
  tab_MSE <- dat_test_long %>% group_by(month) %>% 
    filter(group == "optim") %>% summarise(MSE_dose = mean((dose - best_dose)^2, na.rm = T),
                                           MAE_dose = mean(abs(dose - best_dose), na.rm = T))
  
  list(
    plot_tox = plot_tox,
    plot_tumor = plot_tumor,
    plot_sum = plot_sum,
    plot_reward = plot_reward,
    table_deaths = tab_deaths,
    table_cured = tab_cured,
    table_rewards = tab_reward,
    tab_MSE = tab_MSE
  )
}