# simulation functions ----------------------------------------------------

# function for how toxicity changes
updateW <- function(M, W, D, a1 = 0.1, b1 = 1.2, d1 = 0.5) { 
    W_next <- a1 * M + b1 * (D - d1) + W
    ifelse(W_next > 0, W_next, 0)
}

# function for updating tumor mass
updateM <- function(M, W, D, a2 = 0.15, b2 = 1.2, d2 = 0.5) {
  # takes current value of tumor mass, outputs next value
    if (M > 0) {
      Mdot <- (a2 * W - b2 * (D - d2))
      ifelse(Mdot + M > 0, Mdot + M, 0)
    } else {
      0
    }
}

# reward functions
R2 <- function(W_next, W) {
  -10 * (W_next - W)
}

R3 <- function(M_next, M) {
  ifelse(M_next == 0 & M == 0, 0,
         ifelse(M_next == 0, 30,
                -10 * (M_next - M)))
}

reward <- function(M_next, M, W_next, W, died) {
    R2(W_next, W) +
      R3(M_next, M) +
      ifelse(died == 1, -60, 0)
}

lambda <- function(W, M, mu0 = -7, mu1 = 1, mu2 = 1) {
  exp(mu0 + mu1 * W + mu2 * M)
}

determineDeath <- function(M_next, W_next) {
    lam <- lambda(W_next, M_next)
    deltaF <- exp(-lam)
    p <- 1 - deltaF
    return(rbinom(1, 1, p))
}

sim2 <- function(N = 1000, Ttot = 6) {
  # 1000 patients, 6 months treatment
  
  # choose random initial toxicities and tumor masses
  W0 <- runif(N, min = 0, max = 2)
  M0 <- runif(N, min = 0, max = 2)
  
  # generate 10 noise variables
  X <- replicate(10, runif(N, min = 0, max = 1))
  
  # choose random starting dose > 0.5
  D0 <- runif(N, min = 0.5, max = 1)
  # choose random subsequent doses between 0 and 1 (max tolerable)
  D <- replicate(Ttot - 1, runif(N, min = 0, max = 1))
  D <- cbind(D0, D)
  
  r <- matrix(c(numeric(N * Ttot)), ncol = 6)
  died <- matrix(c(numeric(N * Ttot)), ncol = 6)
  dead <- logical(N)
  W <- matrix(c(W0, numeric(N * Ttot)), ncol = 7)
  M <- matrix(c(M0, numeric(N * Ttot)), ncol = 7)
  
  for (j in 1:Ttot) {
    for (i in 1:N) {
      if (dead[i]) {
        died[i, j] <- 0 
        M[i, j + 1] <- NA
        W[i, j + 1] <- NA
        r[i, j] <- 0
      } else {
        M[i, j + 1] <- updateM(M[i, j], W[i, j], D[i, j])
        W[i, j + 1] <- updateW(M[i, j], W[i, j], D[i, j])
        died[i, j] <- determineDeath(M[i, j + 1], W[i, j + 1])
        r[i, j] <- reward(M[i, j + 1],
                          M[i, j],
                          W[i, j + 1],
                          W[i, j],
                          died[i , j])
      }
      if (died[i, j] == 1) {dead[i] <- T}
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

# set.seed(1)
# dat <- sim()
# set.seed(1)
# dat2 <- sim2()
# 
# all_equal(dat, dat2)


# sim testing function ----------------------------------------------------

sim_test2 <- function(Q) {
  # 200 patients per each of 11 treatments
  N <- 200 * 11
  # 6 months treatment
  Ttot <- 6 
  
  # The initial values of W0 and M0 for the patients were randomly chosen from the 
  # same uniform distribution used in the training data.
  set.seed(5)
  W0 <- runif(N, min = 0, max = 2)
  M0 <- runif(N, min = 0, max = 2)
  
  # treatments consisting of the estimated optimal treatment regime and each of the 10 possible fixed dose levels ranging from 0.1 to 1.0 with increments of size 0.1.
  D1 <- map(seq(from = 0.1, to = 1, by = 0.1), ~ rep(., 200)) %>% flatten_dbl()
  
  # estimate optimal treatment regime for 200
    dat <-
      tibble(
        ID = tail(1:N, 200),
        toxicity = tail(W0, 200),
        tumor_mass = tail(M0, 200)
      )
  
  D0 <-
    max_df(
      data = dat,
      model = Q$mod_list[[1]],
      form = Q$formula,
      idvar = "ID",
      mod_type = Q$mod_type
    )$best
  
  D1 <- c(D1, D0)
  
  D <- replicate(6, D1)
  
  r <- matrix(c(numeric(N * Ttot)), ncol = Ttot)
  died <- matrix(c(numeric(N * Ttot)), ncol = Ttot)
  dead <- logical(N)
  W <- matrix(c(W0, numeric(N * Ttot)), ncol = Ttot + 1)
  M <- matrix(c(M0, numeric(N * Ttot)), ncol = Ttot + 1)
  
  for (j in 1:(Ttot)) {
    for (i in 1:N) {
      if (dead[i]) {
        died[i, j] <- 0
        M[i, j + 1] <- NA
        W[i, j + 1] <- NA
        r[i, j] <- 0
      } else {
        M[i, j + 1] <- updateM(M[i, j], W[i, j], D[i, j])
        W[i, j + 1] <- updateW(M[i, j], W[i, j], D[i, j])
        died[i, j] <- determineDeath(M[i, j + 1], W[i, j + 1])
        r[i, j] <- reward(M[i, j + 1],
                          M[i, j],
                          W[i, j + 1],
                          W[i, j],
                          died[i , j])
      }
      if (died[i, j] == 1) {
        dead[i] <- T
      }
    }
    if (j < Ttot) {
      if (noise_vars == T) {
        df <- tibble(
          ID = tail(1:N, 200),
          toxicity = tail(W[, j + 1], 200),
          tumor_mass = tail(M[, j + 1], 200),
          noise = tail(noise[, j + 1], 200)
        )
      } else {
        df <- tibble(
          ID = tail(1:N, 200),
          toxicity = tail(W[, j + 1], 200),
          tumor_mass = tail(M[, j + 1], 200)
        )
      }
      D[df$ID, j + 1] <- max_df(
        data = df,
        model = Q$mod_list[[j + 1]],
        form = Q$formula,
        mod_type = Q$mod_type
      )$best
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

# set.seed(1)
# t1 <- sim_test(Q)
# set.seed(1)
# t2 <- sim_test2(Q)
# all_equal(t1, t2)
