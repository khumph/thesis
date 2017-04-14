max_df <- function(dat, model, truth, pred, nested = F) {
  dat <- ungroup(dat)
  if (!pred) {
    deads <- filter(dat, is.na(M_next)) %>% mutate(dose = NA)
    dat <- filter(dat, !is.na(M_next))
  }
  
  if (truth) {
    dat <- dat %>% 
      mutate(dose = map(1:nrow(.), ~ seq(0, 1, by = 0.01))) %>%
      unnest()
    dat <- Mnext(dat, truth = truth)
    dat <- Wnext(dat, truth = truth)
    dat <- dat %>% mutate(
      expect_surv_time = 1 / lambda(M_next, W_next, Z = noise_chng),
      preds = log(expect_surv_time)
    )
  } else {
    dat <- dat %>% 
      mutate(dose = map(1:nrow(.), ~ seq(0, 1, by = 0.01))) %>%
      unnest() %>%
      mutate(preds = predict(model, .))
  }
  
  if (!nested) {
    dat <- dat %>% group_by(ID) %>%
      filter(preds == max(preds)) %>% 
      slice(1)
    if (!pred) {
      dat <- dat %>% bind_rows(deads) %>% arrange(ID) %>% ungroup() 
    }
  }
  dat %>% ungroup()
}

Qlearn <- function(form, dat_long, nstages = 6, ...) {
  mod_list <- list()
  for (i in (nstages - 1):0) {
    dat <- filter(dat_long, month == i)
    mod <- train(form, dat, na.action = na.omit, ...)
    new_dat <- max_df(dat, mod, truth = F, pred = F)
    mod_list[[i + 1]] <- mod
    bestR <- new_dat$preds
    bestD <- new_dat$dose
    if (i > 0) {
      dat_long <- dat_long %>% mutate(
        Qhat = ifelse(month == i - 1,
                      reward + ifelse(!is.na(bestR), bestR, 0),
                      Qhat))
    }
    dat_long <- dat_long %>%
      mutate(best = ifelse(month == i, bestD, best))
  }
  list(
    mod_list = mod_list
  )
}
