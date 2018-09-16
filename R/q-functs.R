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
    dat$W_next <- updateW(dat$tumor_mass, dat$toxicity, dat$dose, c = dat$cW)
    dat$M_next <- updateM(dat$tumor_mass, dat$toxicity, dat$dose, c = dat$cM)
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
      filter(preds == max(preds)) %>% slice(1)
    if (!pred) {
      dat <- dat %>% bind_rows(deads) %>% arrange(ID) %>% ungroup()
    }
  }
  dat %>% ungroup()
}

Qlearn <- function(formula, data, n_stages = 3, method, ...) {
  mod_list <- list()
  for (i in (n_stages - 1):0) {
    dat <- filter(data, month == i)
    if (method == "rpart") {
      mod <- rpart(formula, dat, ...)
      mod <- prune(mod, cp = mod$cptable[which.min(mod$cptable[, "xerror"]), "CP"])
    } else {
      mod <- train(form, dat, method = method, na.action = na.omit, ...)
    }
    new_dat <- max_df(dat, mod, truth = F, pred = F)
    mod_list[[i + 1]] <- mod
    bestR <- new_dat$preds
    if (i > 0) {
      data <- data %>% mutate(
        Qhat = ifelse(month == i - 1,
                      reward + ifelse(!is.na(bestR), bestR, 0),
                      Qhat))
    }
  }
  list(
    mod_list = mod_list
  )
}
