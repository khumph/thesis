# simulation functions ----------------------------------------------------

# function for how toxicity changes
Wdot <- function(M, D, X, a1 = 0.1, b1 = 1.2, d1 = 0.5) { 
  a1 * M + b1 * (D - d1)
}

# function for how tumor mass changes
Mdot <- function(M, W, D, X, a2 = 0.15, b2 = 1.2, d2 = 0.5) {
  (a2 * W - b2 * (D - d2) # + ifelse(X > 0.5 & D > 0.5, -1, 0)
  ) * ifelse(M > 0, 1, 0)
}

lambda <- function(W, M, mu0 = -7, mu1 = 1, mu2 = 1) {
  exp(mu0 + mu1 * W + mu2 * M)
}

# reward functions
R2 <- function(W1, W0) {
  -10 * (W1 - W0)
}

R3 <- function(M1, M0) {
  ifelse(M1 == 0 & M0 == 0, 0,
         ifelse(M1 == 0, 30,
                -10 * (M1 - M0)))
}

sim <- function(N = 1000, Ttot = 6) {
  # 1000 patients
  # N <- 1000
  # 6 months treatment
  # Ttot <- 6 
  
  # choose random initial toxicities and tumor masses
  W0 <- runif(N, min = 0, max = 2)
  M0 <- runif(N, min = 0, max = 2)
  
  # choose random baseline characteristic X that doubles the effect on tumor size if greater than 0.5
  # X <- runif(N, min = 0, max = 1)
  
  # choose random starting dose > 0.5
  D0 <- runif(N, min = 0.5, max = 1)
  # choose random subsequent doses between 0 and 1 (max tolerable)
  D <- replicate(5, runif(N, min = 0, max = 1))
  D <- cbind(D0, D)
  
  r <- matrix(c(numeric(6000)), ncol = 6)
  died <- matrix(c(numeric(6000)), ncol = 6)
  W <- matrix(c(W0, numeric(6000)), ncol = 7)
  M <- matrix(c(M0, numeric(6000)), ncol = 7)
  
  for (i in 1:(Ttot)) {
    for (j in 1:N) {
      if (sum(died[j, 1:i]) == 0) {
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
  }
  
  colnames(D) <- 0:5
  colnames(M) <- 0:6
  colnames(W) <- 0:6
  colnames(r) <- 0:5
  colnames(died) <- 1:6
  # colnames(noise) <- 0:6
  dat <- data.frame(D = D, M = M, W = W, r = r, d = died) %>% tbl_df() %>% 
    rownames_to_column(var = "ID") %>% mutate(ID = as.numeric(ID))
  
  dat %>%
    gather(key, value, -ID) %>%
    extract(col = key,
            into = c("var", "month"),
            regex = "(.)\\.(.)") %>%
    spread(var, value) %>%
    select(
      ID,
      month,
      dose = D,
      tumor_mass = M,
      toxicity = W,
      reward = r,
      died = d
    )
}


# sim testing function ----------------------------------------------------

sim_test <- function(Q) {
  if (sum(str_detect(Q$formula$covariates, "noise")) > 0) {
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
