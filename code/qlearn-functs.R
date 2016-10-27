# Q-learning functions ----------------------------------------------------

makeRCS <- function(formula, treatment) {
  form_char <- as.character(formula)
  response <- form_char[2]
  predictors <- form_char[3]
  predictor_names <- strsplit(predictors, " \\+ ")[[1]]
  form_base_rcs <- paste(response,
                         "~",
                         paste0("rcs(", predictor_names, ")", collapse = " + "))
  ints <-
    paste0("rcs(", predictor_names, ")", " %ia% ", "rcs(", treatment, ")")
  trtbytrt <-
    paste0("rcs(", treatment, ")", " %ia% ", "rcs(", treatment, ")")
  ints <- ints[ints != trtbytrt]
  ints <- paste(ints, collapse = " + ")
  formula_rcs <- paste(c(form_base_rcs, ints), collapse = " + ")
  list(formula = as.formula(formula_rcs), predictors = predictor_names, treatment = treatment)
}

subset_df <- function(data, occasion) {
  filter(data, month == occasion)
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


max_df <- function(data, model, form, idvar = NULL, x = seq(0, 1, by = 0.01)) {
  if (is.null(idvar)) {
    data <- data %>% mutate(ID = 1:nrow(.))
    idvar <- "ID"
  }
  data %>%
    select(matches(idvar),
           one_of(form$predictors),
           -matches(form$treatment)) %>%
    mutate(dose = map(1:nrow(.), ~ x)) %>%
    unnest() %>%
    mutate(preds = predict(model, .)) %>%
    group_by_(idvar) %>%
    mutate(max = max(preds),
           best = ifelse(is.na(max),
                               NA,
                               which.max(preds) - 1) / 100) %>%
    select(-dose, -preds) %>% unique()
}

one_step_Q <- function(formula, treatment, data) {
  form <- makeRCS(formula, treatment)
  model <- fit_rcs(form$formula, data)
  new_dat <- max_df(data, model, form)
  list(formula = form, model = model, max = new_dat$max, best = new_dat$best, data_new = new_dat)
}

# subset_df(occasion - 1)$reward + maxQ(occasion, mod)

Qlearn <- function(data, formula, treatment) {
  for (i in 4:0) {
    dat <- subset_df(data, i + 1)
    Q1 <- one_step_Q(formula,
                     data = dat,
                     treatment = treatment)
    mod_list[[i + 2]] <- Q1$model
    data[data$month == (i), ]$Q_hat <- dat$reward + Q1$max
    data[data$month == (i + 1), ]$best <- Q1$best
    if (i == 0) {
      subsetted_df <- subset_df(data, i)
      mod_list[[i + 1]] <- Q1$model
      Q1 <- one_step_Q(formula,
                       data = dat,
                       treatment = treatment)
      data[data$month == i, ]$best <- Q1$best
    }
  }
  list(data = data, mod_list = mod_list, formula = Q1$formula)
}


# testing -----------------------------------------------------------------

# d <- dat_long %>% filter(month == 5)
# form <- makeRCS(Q_hat ~ noise + tumor_mass + dose + toxicity, "dose")
# mod <- fit_rcs(form$formula, data = d)

# Q1 <- one_step_Q(
#   Q_hat ~ noise + tumor_mass + dose + toxicity,
#   data = mutate(dat_long, noise = runif(7000), Q_hat = reward) %>% filter(month == 5),
#   treatment = "dose"
# )
