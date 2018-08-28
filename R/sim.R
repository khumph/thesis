# simulation functions ----------------------------------------------------

updateW <- function(M, W, D, a1 = 0.1, b1 = 1.2, c1, d1 = 0.5, truth) {
  # truth <- rep(truth, length(M))
  W_next <-
    a1 * M + b1 * (c1 * D  - d1) + W #+
    # ifelse(!truth, rnorm(length(M), 0, 0.05), 0)
  ifelse(W_next > 0, W_next, 0)
}

Wnext <- function(dat, truth = F) {
  if (truth) {
    dat %>% mutate(W_next = updateW(tumor_mass, toxicity, dose,
                                    c1 = c1,
                                    truth = truth))
  } else {
    dat %>%
      mutate(W_next = ifelse(!dead,
                             updateW(tumor_mass, toxicity, dose,
                                     c1 = c1,
                                     truth = truth),
                             NA))
  }
}

updateM <- function(M, W, D, a2 = 0.15, b2 = 1.2, c2, d2 = 0.5,
                    X = 0, Z = 0, a3 = 1e-3, truth) {
  # truth <- rep(truth, length(M))
  M_next <- ifelse(M > 0,
                   (a2 * W - b2 * (c2 * D  - d2)) +
                     M + sum(a3 * Z), #+
                     # ifelse(!truth, rnorm(length(M), 0, 0.05), 0),
                   0)
  ifelse(M_next > 0, M_next, 0)
}

Mnext <- function(dat, truth = F) {
  if (truth) {
    dat %>% mutate(M_next = updateM(tumor_mass, toxicity, dose, 
                                    c2 = c2, truth = truth))
  } else {
    dat %>% mutate(M_next = ifelse(!dead,
                                   updateM(tumor_mass, toxicity, dose, 
                                           c2 = c2, truth = truth),
                                   NA))
  }
}

genIntNoise <- function(dat, int, noise, noise_pred) {
  if (int) {
    dat <- dat %>%
      mutate(
        X1 = runif(nrow(.), min = 0, max = 1),
        X2 = runif(nrow(.), min = 0, max = 1),
        c1 = ifelse(X1 > 0.5, 1.5, 1),
        c2 = ifelse(X2 > 0.5, 1.5, 1)
      )
  } else {
    dat <- dat %>% mutate(c1 = 1, c2 = 1)
  }
  if (noise) {
    dat <- dat %>%
      mutate(
        Z1  = rnorm(nrow(.), mean = 1),
        Z2  = rnorm(nrow(.), mean = 1),
        Z3  = rnorm(nrow(.), mean = 1),
        Z4  = rnorm(nrow(.), mean = 1),
        Z5  = rnorm(nrow(.), mean = 1),
        Z6  = rnorm(nrow(.), mean = -1),
        Z7  = rnorm(nrow(.), mean = -1),
        Z8  = rnorm(nrow(.), mean = -1),
        Z9  = rnorm(nrow(.), mean = -1),
        Z10 = rnorm(nrow(.), mean = -1)
      ) %>% 
      bind_cols(
        replicate(90, rnorm(nrow(.))) %>% as.data.frame())
  } 
  if (noise_pred) {
    dat <- dat %>% mutate(
      noise_chng = Z1 + Z2 + Z3 + Z4 + Z5 + Z6 + Z7 + Z8 + Z9 + Z10)
  } else {
    dat <- dat %>% mutate(noise_chng = 0)
  }
  dat
}

# defined as in NSCLC paper
lambda <- function(M, W, Z, mu0 = -5.5, mu1 = 1, mu2 = 1.2, mu3 = 0.75,
                   a3 = 0.05) {
  exp(mu0 + mu1 * W + mu2 * M + mu3 * W * M + a3 * Z)
}


# simulation --------------------------------------------------------------

simMonth <- function(dat) {
  dat <- Mnext(dat)
  dat <- Wnext(dat)
  dat %>% mutate(
    beta = 1 / lambda(M_next, W_next, Z = noise_chng),
    d_next = runif(nrow(.), min = 0, max = 1),
    surv_time = rexp(nrow(.), 1 / beta),
    dead = ifelse(dead, dead, surv_time < 1)
  )
}

sim <- function(N = 1000, Ttot = 6, int = F, noise = F, noise_pred = F, seed = 1) {
  set.seed(seed)
  dat <- tibble(
    ID = 1:N,
    month = rep(0, N),
    tumor_mass = runif(N, min = 0, max = 2),
    toxicity = runif(N, min = 0, max = 2),
    dose = runif(N, min = 0, max = 1),
    dead = rep(F, N)
  )
  
  dat <- genIntNoise(dat, int, noise, noise_pred)
  
  d <- simMonth(dat) %>%
    mutate(reward = ifelse(dead, log(surv_time), 0))
  out <- d
  for (i in 1:(Ttot - 1)) {
    d <- d %>% mutate(month = i,
                      tumor_mass = M_next,
                      toxicity = W_next,
                      dose = d_next) %>%
      simMonth() %>%
      mutate(reward = ifelse(dead, log(i + 1 + surv_time), 0))
    out <- bind_rows(out, d)
  }
  out %>% mutate(
    reward = ifelse(month == Ttot - 1 & !dead,
                    log(surv_time + Ttot - 1), reward),
    Qhat = reward
  )
}
