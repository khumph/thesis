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

max_df <- function(data, model, form, idvar = "ID", mod_type, 
                   x = seq(0, 1, by = 0.01), nested = F) {
    if (is.null(idvar)) {
      data <- data %>% mutate(ID = 1:nrow(.))
      idvar <- "ID"
    }
    data_long <- data %>% filter(!dead) %>% 
      select(matches(idvar),
             one_of(form$covariates)) %>%
      mutate(dose = map(1:nrow(.), ~ x)) %>%
      unnest()
    data_preds <- data_long %>%
      mutate(preds = ifelse(is.na(tumor_mass), 0, predict(model, .))) %>%
      group_by_(idvar) %>%
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


# ex_plots <- function(Q, month, ids) {
#   dat <- align_df(Q)
#   
#   txr_plot <- ggplot(
#     filter(dat, month == month)) +
#       geom_point(aes(x = tumor_mass, y = reward)) +
#       geom_point(aes(x = tumor_mass, y = Q_hat), color = "blue") +
#       # labs(title = paste0("Tumor Mass by Reward month = ", as.character(j))
#       )
#   
#   dat <- filter(dat, month == month)
#   nested_df <-
#     max_df(
#       dat,
#       Q$mod_list[[month + 1]],
#       Q$formula,
#       idvar = "ID",
#       mod_type = "caret",
#       nested = T
#     )
#   
#   pred_plot <- ggplot(data = filter(nested_df, ID %in% ids)) +
#     geom_line(mapping = aes(x = dose, y = preds, group = ID))
#   list(txr_plot, pred_plot)
# }
# 
# ex_plots(Q, ex_plots(Q, 5, sample(1:1000, 30)))

