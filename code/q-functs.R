# Q-learning functions ----------------------------------------------------

makeForm <- function(formula, treatment, mod_type = "rcs") {
  form_char <- as.character(formula)
  response <- form_char[2]
  predictors <- form_char[3]
  predictor_names <- strsplit(predictors, " \\+ ")[[1]]
  if (mod_type == "rcs") {
    form_base_rcs <- paste(response,
                           "~",
                           paste0("rcs(", predictor_names, ")", collapse = " + "))
    ints <-
      paste0("rcs(", predictor_names, ")", " %ia% ", "rcs(", treatment, ")")
    trtbytrt <-
      paste0("rcs(", treatment, ")", " %ia% ", "rcs(", treatment, ")")
    ints <- ints[ints != trtbytrt]
    ints <- paste(ints, collapse = " + ")
    formula <- paste(c(form_base_rcs, ints), collapse = " + ")
  }
  covariates <- predictor_names[!(predictor_names %in% treatment)]
  list(
    formula = as.formula(formula),
    covariates = covariates,
    treatment = treatment
  )
}

max_df <- function(data, model, form, mod_type, 
                   x = seq(0, 1, by = 0.01), nested = F) {
    data_long <- data %>% filter(!dead) %>% 
      select(ID, one_of(form$covariates)) %>%
      mutate(dose = map(1:nrow(.), ~ x)) %>%
      unnest()
    data_preds <- data_long %>%
      mutate(preds = ifelse(is.na(tumor_mass), 0, predict(model, .))) %>%
      group_by(ID) %>%
      mutate(
        max = max(preds),
        best = ifelse(near(preds, max), dose, NA),
        best = ifelse(tumor_mass > 0,
                      max(best, na.rm = T),
                      min(best, na.rm = T))
      )
    if (nested) {
      data_preds
    } else {
      data_preds %>% filter(near(best, dose)) %>%
        bind_rows(filter(data, dead) %>% select(ID))
    }
  }

one_step_Q <- function(formula, treatment, data, mod_type, ...) {
  form <- makeForm(formula, treatment, mod_type)
  if (mod_type == "rcs") {
    require(rms)
    model <- ols(form$formula,
                 x = T,
                 y = T,
                 data,
                 ...)
  }
  if (mod_type == "caret") {
    require(caret)
    model <- train(formula,
                   data,
                   ...)
  }
  new_dat <- max_df(data, model, form, mod_type = mod_type)
  list(
    formula = form,
    model = model,
    max = new_dat$max,
    best = new_dat$best,
    data_new = new_dat
  )
}

Qlearn <- function(data, formula, treatment, mod_type, ...) {
  mod_list <- list()
  for (i in 4:0) {
    dat <- filter(data, month == i + 1)
    Q1 <- one_step_Q(formula, treatment, dat, mod_type, ...)
    mod_list[[i + 2]] <- Q1$model
    data[data$month == i, ]$Q_hat <- data[data$month == i, ]$reward + Q1$max
    data[data$month == (i + 1),]$best <- Q1$best
  }
  Q1 <- one_step_Q(formula, treatment, data = filter(data, month == 0),
                   mod_type, ...)
  mod_list[[1]] <- Q1$model
  data[data$month == i,]$best <- Q1$best
  list(
    data = data,
    mod_list = mod_list,
    formula = Q1$formula,
    mod_type = mod_type
  )
}

# results functions -------------------------------------------------------


indPlot <- function(data) {
  ex_ID <- sample(data$ID, 1)
  
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

maxPlots <- function(Q, mon = 5) {
  dat <- filter(Q$data, month == mon)
  nested_df <-
    max_df(dat,
           Q$mod_list[[mon + 1]],
           Q$formula,
           mod_type = Q$mod_type,
           nested = T)
  
  id <- sample(1:1000, 1)
  onePlot <-
    ggplot(filter(nested_df, ID == id), aes(x = dose, y = preds)) +
    geom_point()
  
  ids <- sample(1:1000, 50)
  manyPlot <- ggplot(data = filter(nested_df, ID %in% ids)) +
    geom_line(mapping = aes(x = dose, y = preds, group = ID))
  
  list(onePlot, manyPlot, id, ids)
}

plots_tab <- function(dat_test_long) {
  dat_long_summ <- dat_test_long %>% group_by(ID) %>%
    mutate(
      tot_reward = sum(reward, na.rm = T),
      cumSurv = prod(1 - pdeath[1:6])
    ) %>% group_by(group, month) %>% 
    summarise(
      mean_tox = mean(toxicity, na.rm = T),
      mean_tumor = mean(tumor_mass, na.rm = T),
      mean_reward = mean(reward, na.rm = T),
      mean_cumSurv = mean(cumSurv, na.rm = T)
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
  
  plot_reward <- ggplot(data = dat_long_summ) +
    geom_line(mapping = aes(
      x = month,
      y = mean_reward,
      color = group,
      group = group
    ))
  
  tab_deaths <- dat_long_summ %>%
    group_by(group) %>%
    summarise(mean_cumSurv = mean(mean_cumSurv)) %>%
    arrange(desc(mean_cumSurv))
  
  tab_reward <- dat_long_summ %>%
    group_by(group) %>%
    summarise(mean_tot_reward = sum(mean_reward, na.rm = T)) %>%
    arrange(desc(mean_tot_reward))
  
  list(
    plot_tox = plot_tox,
    plot_tumor = plot_tumor,
    plot_sum = plot_sum,
    plot_reward = plot_reward,
    table_deaths = tab_deaths,
    table_rewards = tab_reward
  )
}