# Q-learning functions ----------------------------------------------------

makeForm <- function(formula, treatment, mod_type = "rcs") {
  form_char <- as.character(formula)
  response <- form_char[2]
  predictors <- form_char[3]
  predictor_names <- strsplit(predictors, " \\+ ")[[1]]
  covariates <- predictor_names[!(predictor_names %in% treatment)]
  
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
  
  list(
    formula = as.formula(formula),
    covariates = covariates,
    treatment = treatment
  )
}

max_df <- function(data, model, form, mod_type, 
                   x = seq(0, 1, by = 0.01), nested = F) {
  dat <- data %>%
    filter(!is.na(tumor_mass)) %>% 
    mutate(dose = map(1:nrow(.), ~ x)) %>%
    unnest() %>%
    mutate(preds = predict(model, .)) %>%
    group_by(ID) %>%
    mutate(
      max = max(preds),
      best = ifelse(near(preds, max), dose, NA),
      best = ifelse(tumor_mass > 0,
                    quantile(best, probs = 0, na.rm = T, type = 3, names = F),
                    min(best, na.rm = T))
    )
  if (nested) {
    dat
  } else {
    dat %>% filter(near(best, dose)) %>%
      bind_rows(filter(data, is.na(tumor_mass)))
  }
}

one_step_Q <- function(form, data, mod_type, ...) {
  
  if (mod_type == "rcs") {
    model <- ols(form$formula,
                 x = T,
                 y = T,
                 data,
                 ...)
  }
  
  if (mod_type == "caret") {
    model <- train(form$formula,
                   data,
                   ...)
  }
  
  new_dat <- max_df(data, model, form, mod_type = mod_type)
  
  list(
    model = model,
    max = new_dat$max,
    best = new_dat$best
  )
}

Qlearn <- function(data, formula, treatment, mod_type, ...) {
  form <- makeForm(formula, treatment, mod_type)
  
  mod_list <- list()
  for (i in 4:0) {
    Q1 <- one_step_Q(form, filter(data, month == i + 1), mod_type, ...)
    mod_list[[i + 2]] <- Q1$model
    data[data$month == i, ]$Q_hat <- data[data$month == i, ]$reward + Q1$max
    data[data$month == (i + 1), ]$best <- Q1$best
  }
  
  Q1 <- one_step_Q(form, filter(data, month == 0), mod_type, ...)
  mod_list[[1]] <- Q1$model
  data[data$month == i,]$best <- Q1$best
  
  list(
    data = data,
    mod_list = mod_list,
    formula = form,
    mod_type = mod_type
  )
}
