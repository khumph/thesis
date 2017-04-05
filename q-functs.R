# Q-learning functions ----------------------------------------------------

max_df <- function(data, model,
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

one_step_Q <- function(form, data, ...) {
  
  model <- train(form,
                 data,
                 na.action = na.omit,
                 ...)
  
  new_dat <- max_df(data, model)
  
  list(
    data = new_dat,
    model = model,
    max = new_dat$max,
    best = new_dat$best
  )
}

Qlearn <- function(data, form, boot = F, nstages = 6, ...) {
  
  mod_list <- list()
  for (i in (nstages - 1):0) {
    Q1 <- one_step_Q(form, filter(data, month == i), ...)
    mod_list[[i + 1]] <- Q1$model
    if (i > 0) {
      data <- data %>% mutate(
        Qhat = ifelse(month == i - 1,
                      reward + ifelse(!is.na(Q1$max), Q1$max, 0),
                      Qhat))
    }
    data <- data %>% mutate(best = ifelse(month == i, Q1$best, best))
  }
  
  list(
    data = data,
    mod_list = mod_list
  )
}
