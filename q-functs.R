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
                   x = seq(0, 1, by = 0.01), nested = F, pred = F) {
  data <- ungroup(data)
  if (pred) {
    dat <- data
  } else {
    dat <- data %>% filter(!is.na(M_next))
  }
  dat <- dat %>% 
    mutate(dose = map(1:nrow(.), ~ x)) %>%
    unnest() %>%
    mutate(preds = predict(model, .)) 
  if (!nested) {
    dat <- dat %>% group_by(ID) %>%
      mutate(
        max = max(preds),
        best = ifelse(near(preds, max), dose, NA),
        best = quantile(best, probs = 0, na.rm = T, type = 3, names = F)
        # best = ifelse(tumor_mass > 0,
        #               quantile(best, probs = 1, na.rm = T, type = 3, names = F),
        #               quantile(best, probs = 0, na.rm = T, type = 3, names = F))
      )
  }
  if (nested) {
    dat %>% ungroup()
  } else if (!pred) {
    dat %>% filter(near(best, dose)) %>%
      bind_rows(filter(data, is.na(M_next))) %>%
      arrange(ID) %>% ungroup()
  } else {
    dat %>% filter(near(best, dose)) %>%
      arrange(ID) %>% ungroup()
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
                   na.action = na.omit,
                   ...)
  }
  
  new_dat <- max_df(data, model, form, mod_type = mod_type)
  
  list(
    data = new_dat,
    model = model,
    max = new_dat$max,
    best = new_dat$best
  )
}

Qlearn <- function(data, formula, treatment, mod_type, boot = F, nstages = 6, ...) {
  form <- makeForm(formula, treatment, mod_type)
  
  mod_list <- list()
  for (i in (nstages - 1):1) {
    Q1 <- one_step_Q(form, filter(data, month == i), mod_type, ...)
    mod_list[[i + 1]] <- Q1$model
    data <- data %>% mutate(
      Qhat = ifelse(month == i - 1,
                    reward + ifelse(!is.na(Q1$max), Q1$max, 0),
                    Qhat),
      best = ifelse(month == i, Q1$best, best)
    )
  }
  Q1 <- one_step_Q(form, filter(data, month == 0), mod_type, ...)
  mod_list[[1]] <- Q1$model
  data <- data %>% mutate(
    best = ifelse(month == 0, Q1$best, best)
  )
  
  list(
    data = data,
    mod_list = mod_list,
    formula = form,
    mod_type = mod_type
  )
}
