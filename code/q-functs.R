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
    data[data$month == i, ]$Q_hat <- dat$reward + Q1$max
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

