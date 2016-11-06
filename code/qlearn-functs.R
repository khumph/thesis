# Q-learning functions ----------------------------------------------------

makeRCS <- function(formula, treatment, method = "rcs") {
  form_char <- as.character(formula)
  response <- form_char[2]
  predictors <- form_char[3]
  predictor_names <- strsplit(predictors, " \\+ ")[[1]]
  if (method == "rcs") {
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

fit_rpart <- function(formula, data, cpmethod = "min", ...) {
  mod_rpart <- rpart(formula, 
                     data = data, ...)
  if (cpmethod == "min_sd") {
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
  if (cpmethod == "min") {
    pruned <- mod_rpart$cptable %>%
      as_data_frame() %>%
      filter(xerror == min(xerror)) %>% 
      filter(nsplit == min(nsplit)) %>% 
      select(CP) %>% flatten_dbl() %>%
      prune(mod_rpart, .)
  }
  pruned
}


max_df <- function(data, model, form, idvar = NULL, x = seq(0, 1, by = 0.01), method = "rcs") {
  if (is.null(idvar)) {
    data <- data %>% mutate(ID = 1:nrow(.))
    idvar <- "ID"
  }
  data_long <- data %>%
    select(matches(idvar),
           one_of(form$covariates)) %>%
    mutate(dose = map(1:nrow(.), ~ x)) %>%
    unnest() 
  if (method == "rcs") {
    data_preds <- data_long %>%
      mutate(preds = predict(model, .))
  }
  if (method == "rpart") {
    data_preds <- data_long %>%
      mutate(preds = predict(model, .) %>%
               apply(1, function(z) as.numeric(names(which.max(z)))))
  }
  data_preds %>%
    group_by_(idvar) %>%
    mutate(max = max(preds),
           best = ifelse(is.na(max),
                               NA,
                               which.max(preds) - 1) / 100) %>%
    select(-dose, -preds) %>% unique()
}

one_step_Q <- function(formula, treatment, data, method = "rcs", ...) {
  form <- makeRCS(formula, treatment, method = method) 
  if (method == "rcs") {
    model <- fit_rcs(form$formula, data, ...)
  }
  if (method == "rpart") {
    data <- data %>% mutate(Q_hat = factor(Q_hat))
    model <- fit_rpart(form$formula,
                       data,
                       ...)
  }
  new_dat <- max_df(data, model, form, method = method)
  list(formula = form, model = model, max = new_dat$max, best = new_dat$best, data_new = new_dat)
}

# subset_df(occasion - 1)$reward + maxQ(occasion, mod)

Qlearn <- function(data, formula, treatment, method = "rcs", ...) {
  mod_list <- list()
  for (i in 4:0) {
    dat <- subset_df(data, i + 1)
    Q1 <- one_step_Q(formula,
                     data = dat,
                     treatment = treatment,
                     method = method,
                     ...)
    mod_list[[i + 2]] <- Q1$model
    if (method == "rcs") {
      data[data$month == (i), ]$Q_hat <- dat$reward + Q1$max
    }
    if (method == "rpart") {
      data[data$month == (i), ]$Q_hat <- dat$reward + Q1$max
    }
    data[data$month == (i + 1), ]$best <- Q1$best
    if (i == 0) {
      subsetted_df <- subset_df(data, i)
      mod_list[[i + 1]] <- Q1$model
      Q1 <- one_step_Q(formula,
                       data = dat,
                       treatment = treatment, 
                       method = method)
      data[data$month == i, ]$best <- Q1$best
    }
  }
  list(data = data, mod_list = mod_list, formula = Q1$formula, method = method)
}


# testing -----------------------------------------------------------------

# data <- dat_long %>% filter(month == 5) %>% mutate(Q_hat = factor(Q_hat))
# formula <- Q_hat ~ tumor_mass + dose + toxicity
# method <- "rcs"
# treatment <- "dose"
# model <- fit_rpart(form$formula, data = data)
# form <- makeRCS(formula = formula, treatment = "dose", method = "rpart")
# idvar <- "ID"
# method = "rpart"
# x <- seq(0,1,0.01)
# max_df(data = data, model = model, form, method = "rpart")
# 
# if (is.null(idvar)) {
#   data <- data %>% mutate(ID = 1:nrow(.))
#   idvar <- "ID"
# }
# data_long <- data %>%
#   select(matches(idvar),
#          one_of(form$covariates)) %>%
#   mutate(dose = map(1:nrow(.), ~ x)) %>%
#   unnest() 
# if (method == "rcs") {
#   data_preds <- data_long %>%
#     mutate(preds = predict(model, .))
# }
# if (method == "rpart") {
#   data_preds <- data_long %>%
#     mutate(preds = predict(model, .) %>%
#              apply(1, function(x) as.numeric(names(which.max(x)))))
# }
# 
# predict(model, data_long)
# 
# 
# data_preds %>%
#   group_by_(idvar) %>%
#   mutate(max = max(preds),
#          best = ifelse(is.na(max),
#                        NA,
#                        which.max(preds) - 1) / 100) %>%
#   select(-dose, -preds) %>% unique()


# 
# # Q1 <- one_step_Q(
# #   Q_hat ~ noise + tumor_mass + dose + toxicity,
# #   data = mutate(dat_long, noise = runif(7000), Q_hat = reward) %>% filter(month == 5),
# #   treatment = "dose"
# # )
# 
#   form_char <- as.character(formula)
#   response <- form_char[2]
#   predictors <- form_char[3]
#   predictor_names <- strsplit(predictors, " \\+ ")[[1]]
#   if (method == "rcs") {
#     form_base_rcs <- paste(response,
#                            "~",
#                            paste0("rcs(", predictor_names, ")", collapse = " + "))
#     ints <-
#       paste0("rcs(", predictor_names, ")", " %ia% ", "rcs(", treatment, ")")
#     trtbytrt <-
#       paste0("rcs(", treatment, ")", " %ia% ", "rcs(", treatment, ")")
#     ints <- ints[ints != trtbytrt]
#     ints <- paste(ints, collapse = " + ")
#     formula <- paste(c(form_base_rcs, ints), collapse = " + ")
#   }
#   covariates <- predictor_names[!(predictor_names %in% treatment)]
#   list(formula = as.formula(formula), covariates = covariates, treatment = treatment)
# 
# data %>%
#   select(matches(idvar),
#          one_of(form$covariates)) %>%
#   mutate(dose = map(1:nrow(.), ~ x)) %>%
#   unnest() %>%
#   mutate(preds = predict(model, .)) %>%
#   group_by_(idvar) %>%
#   mutate(max = max(preds),
#          best = ifelse(is.na(max),
#                        NA,
#                        which.max(preds) - 1) / 100) %>%
#   select(-dose, -preds) %>% unique()


# sim testing function ----------------------------------------------------

sim_test <- function(Q) {
  if (str_detect(Q$formula$covariates, "noise")) {
    noise_vars <- T
  } else {
    noise_vars <- F
  }
  # 200 patients per each of 11 treatments
  N <- 200 * 11
  # 6 months treatment
  Ttot <- 6 
  
  # The initial values of W0 and M0 for the patients were randomly chosen from the 
  # same uniform distribution used in the training data.
  set.seed(5)
  W0 <- runif(N, min = 0, max = 2)
  M0 <- runif(N, min = 0, max = 2)
  if (noise_vars == T) {
    noise <- replicate(Ttot, runif(N, min = 0, max = 2))
  }
  
  # treatments consisting of the estimated optimal treatment regime and each of the 10 possible fixed dose levels ranging from 0.1 to 1.0 with increments of size 0.1.
  D1 <- map(seq(from = 0.1, to = 1, by = 0.1), ~ rep(., 200)) %>% flatten_dbl()
  
  # estimate optimal treatment regime for 200
  if (noise_vars == T) {
    dat <-
      tibble(
        ID = tail(1:N, 200),
        toxicity = tail(W0, 200),
        tumor_mass = tail(M0, 200),
        noise = tail(noise[, 1], 200)
      )
  } else {
    dat <-
      tibble(
        ID = tail(1:N, 200),
        toxicity = tail(W0, 200),
        tumor_mass = tail(M0, 200)
      )
  }
  
  D0 <- max_df(data = dat, model = Q$mod_list[[1]], form = Q$formula, idvar = "ID", method = Q$method)$best
  
  D1 <- c(D1, D0)
  
  D <- replicate(6, D1)
  
  # function for how toxicity changes
  Wdot <- function(M, D, a1 = 0.1, b1 = 1.2, d1 = 0.5) { 
    a1 * M + b1 * (D - d1)
  }
  
  # function for how tumor mass changes
  Mdot <- function(M, W, D, a2 = 0.15, b2 = 1.2, d2 = 0.5) { 
    (a2 * W - b2 * (D - d2)) * ifelse(M > 0, 1, 0)
  }
  
  # reward functions
  R2 <- function(W1, W0) {
    ifelse(W1 - W0 <= -0.5, 5,
           ifelse(W1 - W0 >= 0.5, -5, 0))
  }
  R3 <- function(M1, M0) {
    ifelse(M1 == 0 & M0 == 0, 0,
           ifelse(M1 == 0, 15,
                  ifelse(M1 - M0 <= -0.5, 5,
                         ifelse(M1 - M0 >= 0.5, -5, 0)))
    )
  }
  
  lambda <- function(W, M, mu0 = -7, mu1 = 1, mu2 = 1) {
    exp(mu0 + mu1 * W + mu2 * M)
  }
  
  r <- matrix(c(numeric(N * Ttot)), ncol = Ttot)
  died <- matrix(c(numeric(N * Ttot)), ncol = Ttot)
  W <- matrix(c(W0, numeric(N * Ttot)), ncol = Ttot + 1)
  M <- matrix(c(M0, numeric(N * Ttot)), ncol = Ttot + 1)
  
  for (i in 1:(Ttot)) {
    for (j in 1:N) {
      if (sum(died[j, 1:i], na.rm = T) == 0) {
        M_next <- ifelse(M[j, i] == 0, # if patient cured (tumor mass == 0)
                         0, # mass stays 0
                         Mdot(M[j, i], W[j, i], D[j, i]) + M[j, i]) # otherwise find next mass
        W_next <- Wdot(M[j, i], D[j, i]) + W[j, i]
        # mass and toxicity bounded at 0
        W[j, i + 1] <- ifelse(W_next <= 0,
                              0,
                              W_next)
        M[j, i + 1] <- ifelse(M_next <= 0,
                              0,
                              M_next)
        # determine if patient died
        lam <- lambda(W[j, i + 1], M[j, i + 1])
        deltaF <- exp(-lam)
        p <- 1 - deltaF
        died[j, i] <- rbinom(1, 1, p)
        # add up rewards
        r[j, i] <-
          R2(W[j, i + 1], W[j, i]) +
          R3(M[j, i + 1], M[j, i]) +
          ifelse(died[j, i] == 1, -60, 0)
      } else {
        # if patient already died, can't die again or get rewards, mass and tox to NA
        r[j, i] <- 0
        if (i != Ttot) {
          died[j, i + 1] <- 0
        } else {died[j, i] <- 0}
        M[j, i + 1] <- NA
        W[j, i + 1] <- NA
      }
    }
    if (i < Ttot) {
      if (noise_vars == T) {
        df <- tibble(
          ID = tail(1:N, 200),
          toxicity = tail(W[, i + 1], 200),
          tumor_mass = tail(M[, i + 1], 200),
          noise = tail(noise[, i + 1], 200)
        )
      } else {
        df <- tibble(
          ID = tail(1:N, 200),
          toxicity = tail(W[, i + 1], 200),
          tumor_mass = tail(M[, i + 1], 200)
        )
      }
      D[df$ID, i + 1] <- max_df(data = df,
                                model = Q$mod_list[[i + 1]],
                                form = Q$formula,
                                method = Q$method)$best
    }
  }
  
  colnames(D) <- 0:5
  colnames(M) <- 0:6
  colnames(W) <- 0:6
  colnames(r) <- 0:5
  colnames(died) <- 1:6
  if (noise_vars == T) {
    colnames(noise) <- 0:5
    dat_test <-
      data.frame(
        D = D,
        M = M,
        W = W,
        r = r,
        d = died,
        n = noise 
      ) %>% tbl_df()
  } else {
    dat_test <-
      data.frame(
        D = D,
        M = M,
        W = W,
        r = r,
        d = died
      ) %>% tbl_df()
  }
  
  dat_test <- dat_test %>%
    rownames_to_column(var = "ID") %>%
    mutate(ID = as.numeric(ID),
           group = ifelse(ID < 2001, D1, "optim"),
           group = factor(group))
  
  dat_test <- dat_test %>%
    gather(key, value, -ID, -group) %>%
    extract(col = key,
            into = c("var", "month"),
            regex = "(.)\\.(.)") %>%
    spread(var, value)
  
  if (noise_vars == T) {
    dat_test %>%
      select(
        ID,
        month,
        group,
        noise = n,
        dose = D,
        tumor_mass = M,
        toxicity = W,
        reward = r,
        died = d
      )
  } else {
    dat_test %>%
      select(
        ID,
        month,
        group,
        dose = D,
        tumor_mass = M,
        toxicity = W,
        reward = r,
        died = d
      )
  }
}



