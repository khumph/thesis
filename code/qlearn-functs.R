 # Q-learning functions ----------------------------------------------------

makeRCS <- function(formula, treatment, mod_type = "rcs") {
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
  list(formula = as.formula(formula), covariates = covariates, treatment = treatment)
}

fit_rcs <- function(formula, data, ...) {
  ols(
    formula,
    x = T,
    y = T,
    data,
    ...
  )
}

fit_rpart <- function(formula, data, cpmod_type = "min", ...) {
  mod_rpart <- rpart(formula, 
                     data = data, mod_type = "class", ...)
  if (cpmod_type == "min_sd") {
    pruned <- mod_rpart$cptable %>%
      as_data_frame() %>%
      mutate(
        xerror_min = min(xerror),
        xerror_min_sd = ifelse(xerror == xerror_min, xstd, 0) %>% sum(na.rm = T)
      ) %>% filter(near(xerror, xerror_min + xerror_min_sd, tol = 0.01)) %>%
      filter(nsplit == min(nsplit)) %>% 
      select(CP) %>% flatten_dbl() %>%
      prune(mod_rpart, .)
  }
  if (cpmod_type == "min") {
    pruned <- mod_rpart$cptable %>%
      as_data_frame() %>%
      filter(xerror == min(xerror)) %>% 
      filter(nsplit == min(nsplit)) %>% 
      select(CP) %>% flatten_dbl() %>%
      prune(mod_rpart, .)
  }
  pruned
}


max_df <- function(data, model, form, idvar = NULL, mod_type, x = seq(0, 1, by = 0.01),  nested = F) {
  if (is.null(idvar)) {
    data <- data %>% mutate(ID = 1:nrow(.))
    idvar <- "ID"
  }
  data_long <- data %>%
    select(matches(idvar),
           one_of(form$covariates)) %>%
    mutate(dose = map(1:nrow(.), ~ x)) %>%
    unnest()
  if (mod_type == "rpart") {
    data_preds <- data_long %>%
      mutate(preds = predict(model, .) %>%
               apply(1, function(z) as.numeric(names(which.max(z))))) %>% 
      group_by_(idvar) %>%
      mutate(max = max(preds),
             best = ifelse(preds == max, dose, NA) %>% median(na.rm = T)
      )
  }
  if (mod_type == "ranger") {
    data_preds <- data_long %>%
      mutate(preds = ifelse(
        !is.na(tumor_mass),
        predict(model, data_long %>% na.omit()) %>% predictions(),
        0
      )) %>%
      group_by_(idvar) %>%
      mutate(max = max(preds),
             best = (which.max(preds) - 1) / 100) 
  } else {
    data_preds <- data_long %>%
      mutate(preds = ifelse(is.na(tumor_mass), 0, predict(model, .))) %>%
      group_by_(idvar) %>%
      mutate(max = max(preds),
             best = ifelse(near(preds, max), dose, NA) %>% min(na.rm = T))
  }
  if (nested == T) {
    data_preds
  } else {
    data_preds %>%
      select(-dose, -preds) %>% unique()  
  }
}

one_step_Q <- function(formula, treatment, data, mod_type, ...) {
  form <- makeRCS(formula, treatment, mod_type) 
  if (mod_type == "rcs") {
    model <- fit_rcs(form$formula, data, ...)
  }
  if (mod_type == "caret") {
    model <- train(formula,
                   data,
                   ...)
  }
  if (mod_type == "mars") {
    model <- earth(formula = formula,
                   data = na.omit(data),
                   ...)
  }
  if (mod_type == "rpart") {
    data <- data %>% mutate(Q_hat = factor(Q_hat))
    model <- fit_rpart(form$formula,
                       data,
                       ...)
  }
  if (mod_type == "ert") {
    model <- train(form$formula,
                   data,
                   method = 'extraTrees',
                   na.action = na.omit,
                   tuneGrid = expand.grid(mtry = 3,
                                          numRandomCuts = 1:3),
                   ntree = 50, # G in CRT paper
                   nodesize = 2, # nmin in CRT paper
                   ...)
  }
  if (mod_type == "ranger") {
    model <- ranger(formula,
                    na.omit(data),
                    always.split.variables = treatment,
                    ...)
  }
  new_dat <- max_df(data, model, form, mod_type = mod_type)
  list(formula = form, model = model, max = new_dat$max, best = new_dat$best, data_new = new_dat)
}

Qlearn <- function(data, formula, treatment, mod_type, ...) {
  mod_list <- list()
  for (i in 4:0) {
    dat <- filter(data, month == i + 1)
    Q1 <- one_step_Q(formula = formula,
                     data = dat,
                     treatment = treatment,
                     mod_type = mod_type,
                     ...)
    mod_list[[i + 2]] <- Q1$model
    data[data$month == i, ]$Q_hat <- dat$reward + Q1$max
    data[data$month == (i + 1), ]$best <- Q1$best
    if (i == 0) {
      subsetted_df <- filter(data, month == i)
      mod_list[[i + 1]] <- Q1$model
      Q1 <- one_step_Q(formula,
                       data = dat,
                       treatment = treatment, 
                       mod_type = mod_type)
      data[data$month == i, ]$best <- Q1$best
    }
  }
  list(data = data, mod_list = mod_list, formula = Q1$formula, mod_type = mod_type)
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
    mutate(reward = lag(reward),
           Q_hat = lag(Q_hat),
           best = ifelse(died != 1 | is.na(died), best, NA))
}

plots_tab <- function(dat_test_long) {
  dat_long_summ <- dat_test_long %>% group_by(group, month) %>% 
    summarise(
      mean_tox = mean(toxicity, na.rm = T),
      mean_tumor = mean(tumor_mass, na.rm = T)
    ) %>% mutate(
      sum_means = mean_tox + mean_tumor
    )
  
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
