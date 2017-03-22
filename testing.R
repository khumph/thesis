set.seed(1)
meow <- sim_test(Q, int = int, noise = noise) %>% filter(group == "best")
meow  %>% View()

ids.1 <- meow %>% filter(group == "optim") %>% select(ID) %>% flatten_dbl()
ids1 <- meow %>% filter(group == "1on1off") %>% select(ID) %>% flatten_dbl()
idsbest <- meow %>% filter(group == "best") %>% select(ID) %>% flatten_dbl()
indPlot(meow, sample(idsbest, 1))


dat_long_summ <- meow %>% group_by(ID) %>%
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

tab_MSE <- meow %>% group_by(month) %>% 
  filter(group == "optim") %>% summarise(MSE_dose = mean((dose - best_dose)^2, na.rm = T),
                                         MAE_dose = mean(abs(dose - best_dose), na.rm = T))
