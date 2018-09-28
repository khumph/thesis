makeFormula <- function(data) {
  predictor_names <- c("tumor_mass", "toxicity",
                       grep("\\X|\\Z|\\V", names(data), value = T))
  treatment_name <- "dose"
  formula <- paste("Qhat ~ ", paste(c(treatment_name, predictor_names),
                                    collapse = " + "))
  return(list(
    formula = as.formula(formula),
    predictor_names = predictor_names
  ))
}

Qlearn <- function(formula, data, n_stages = 3, method, ...) {

  mod_list <- list()
  for (i in (n_stages - 1):0) {

    dat <- data[data$month == i, ]
    dat <- data.table::as.data.table(dat)

    # remove people who are already dead at this stage, can't optimize dose if dead
    dat <- dat[!is.na(dat$M_next), ]

    if (method == "rpart") {
      mod <- rpart::rpart(formula, dat, ...)
      mod <- rpart::prune(mod, cp = mod$cptable[which.min(mod$cptable[, "xerror"]), "CP"])
    } else {
      mod <- caret::train(formula, dat, method = method, na.action = na.omit, ...)
    }
    mod_list[[paste0('month', i)]] <- mod

    n <- nrow(dat)
    doses <- seq(from = 0, to = 1, by = 0.01)
    dat_expand <- dat[rep(seq_len(n), each = length(doses)), ]
    dat_expand[ , dose := rep(doses, n)]
    dat_expand[ , pred := predict(mod, dat_expand)]
    max_preds <- dat_expand[ , max(pred), by = ID]

    if (i > 0) {
      ind <- (data$month == i - 1) & (data$ID %in% max_preds$ID)
      data[ind, ]$Qhat <- data[ind, ]$reward + max_preds$V1
    }
  }
  return(mod_list)
}
