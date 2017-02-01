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

max_df <-
  function(data,
           model,
           form,
           idvar = NULL,
           mod_type,
           x = seq(0, 1, by = 0.01),
           nested = F) {
    if (is.null(idvar)) {
      data <- data %>% mutate(ID = 1:nrow(.))
      idvar <- "ID"
    }
    data_long <- data %>%
      select(matches(idvar),
             one_of(form$covariates)) %>%
      mutate(dose = map(1:nrow(.), ~ x)) %>%
      unnest()
    data_preds <- data_long %>%
      mutate(preds = ifelse(is.na(tumor_mass), 0, predict(model, .))) %>%
      group_by_(idvar) %>%
      mutate(max = max(preds),
             best = ifelse(near(preds, max), dose, NA) %>% min(na.rm = T))
    if (nested == T) {
      data_preds
    } else {
      data_preds %>%
        select(-dose, -preds) %>% unique()
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
    Q1 <- one_step_Q(
      formula = formula,
      data = dat,
      treatment = treatment,
      mod_type = mod_type,
      ...
    )
    mod_list[[i + 2]] <- Q1$model
    data[data$month == i, ]$Q_hat <- dat$reward + Q1$max
    data[data$month == (i + 1),]$best <- Q1$best
    if (i == 0) {
      Q1 <- one_step_Q(
        formula = formula,
        data = filter(data, month == i),
        treatment = treatment,
        mod_type = mod_type,
        ...
      )
      mod_list[[i + 1]] <- Q1$model
      data[data$month == i,]$best <- Q1$best
    }
  }
  list(
    data = data,
    mod_list = mod_list,
    formula = Q1$formula,
    mod_type = mod_type
  )
}

# Getting results functions -----------------------------------------------

align_df <- function(Q) {
  Q$data %>%
    select(ID,
           month,
           dose,
           best,
           reward,
           Q_hat,
           tumor_mass,
           toxicity,
           died) %>%
    mutate(
      reward = lag(reward),
      Q_hat = lag(Q_hat),
      best = ifelse(died != 1 | is.na(died), best, NA)
    )
}

plots_tab <- function(dat_test_long) {
  dat_long_summ <- dat_test_long %>% group_by(group, month) %>%
    summarise(
      mean_tox = mean(toxicity, na.rm = T),
      mean_tumor = mean(tumor_mass, na.rm = T)
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
  
  tab <-
    dat_test_long %>% group_by(group) %>%
    summarise(num_died = sum(died, na.rm = T),
              prop_died = num_died / 200)
  
  list(
    plot_tox = plot_tox,
    plot_tumor = plot_tumor,
    plot_sum = plot_sum,
    table_deaths = tab
  )
}
