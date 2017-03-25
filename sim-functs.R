# simulation functions ----------------------------------------------------

updateW <- function(M, W, D, a1 = 0.1, b1, d1 = 0.5,
                    Z = 0, a3 = 1e-3, truth) {
  truth <- rep(truth, length(M))
  W_next <-
    a1 * M + b1 * (D - d1) + W + sum(a3 * Z) + ifelse(!truth, rnorm(length(M), 0, 0.05), 0)
  ifelse(W_next > 0, W_next, 0)
}

Wnext <- function(dat, int, noise_pred, b1 = 1.2, truth = F) {
  if (int) {
    dat %>%
      mutate(W_next = ifelse(!dead,
                             updateW(tumor_mass, toxicity, dose,
                                     b1 = ifelse(X1 > 0.5, 1.5, 1), truth = truth),
                             NA))
  } else if (noise_pred) {
    dat %>%
      mutate(W_next = ifelse(!dead,
                             updateW(tumor_mass, toxicity, dose, b1 = b1,
                                     Z = c(Z7, Z8, Z9, Z10, Z11, Z12), truth = truth),
                             NA))
  } else {
    dat %>%
      mutate(W_next = ifelse(!dead, updateW(tumor_mass, toxicity, dose,
                                            b1 = b1, truth = truth), NA))
  }
}

updateM <- function(M, W, D, a2 = 0.15, b2 = 1.2, d2 = 0.5,
                    X = 0, Z = 0, a3 = 1e-3, truth) {
  truth <- rep(truth, length(M))
  M_next <- ifelse(M > 0,
                   (a2 * W - b2 * (D  - d2)) +
                     M + sum(a3 * Z) + ifelse(!truth, rnorm(length(M), 0, 0.05), 0),
                   0)
  ifelse(M_next > 0, M_next, 0)
}

# (5 + dat$month[1])/10
Mnext <- function(dat, int, noise_pred, b2 = 1.2, truth = F) {
  if (int) {
    dat %>%
      mutate(M_next = ifelse(!dead,
                             updateM(tumor_mass, toxicity, dose, 
                                     b2 = ifelse(X2 > 0.5, 1.5, 1), truth = truth),
                             NA))
  } else if (noise_pred) {
    dat %>%
      mutate(M_next = ifelse(!dead,
                             updateM(tumor_mass, toxicity, dose, b2 = b2,
                                     Z = c(Z1, Z2, Z3, Z4, Z5, Z6), truth = truth),
                             NA))
  } else {
    dat %>%
      mutate(M_next = ifelse(!dead, updateM(tumor_mass, toxicity, dose, b2 = b2,
                                            truth = truth), NA))
  }
}

genIntNoise <- function(dat, int, noise) {
  if (int) {
    dat %>%
      mutate(
        X1 = runif(nrow(.), min = 0, max = 1),
        X2 = runif(nrow(.), min = 0, max = 1)
      )
  } else if (noise) {
    dat %>% mutate(
      Z1  = rnorm(nrow(.), mean = 0.05, sd = 0.25),
      Z2  = rnorm(nrow(.), mean = 0.05, sd = 0.25),
      Z3  = rnorm(nrow(.), mean = 0.05, sd = 0.25),
      Z4  = rnorm(nrow(.), mean = -0.05, sd = 0.25),
      Z5  = rnorm(nrow(.), mean = -0.05, sd = 0.25),
      Z6  = rnorm(nrow(.), mean = -0.05, sd = 0.25),
      Z7  = rnorm(nrow(.), mean = 0.05, sd = 0.25),
      Z8  = rnorm(nrow(.), mean = 0.05, sd = 0.25),
      Z9  = rnorm(nrow(.), mean = 0.05, sd = 0.25),
      Z10 = rnorm(nrow(.), mean = -0.05, sd = 0.25),
      Z11 = rnorm(nrow(.), mean = -0.05, sd = 0.25),
      Z12 = rnorm(nrow(.), mean = -0.05, sd = 0.25)
    ) %>% 
      bind_cols(
        replicate(90, rnorm(nrow(.))) %>% as.data.frame())
  } else {
    dat
  }
}

# defined as in NSCLC paper
lambda <- function(M, W, mu0 = -10, mu1 = 1.75, mu2 = 2, mu3 = 2) {
  exp(mu0 + mu1 * W + mu2 * M + mu3 * W * M)
}

# our new definition
# lambda <- function(M, W, mu0 = -6.5, mu1 = 1, mu2 = 1) {
#   g <- (1)*(M - W)^2
#   exp(mu0 + mu1 * W + mu2 * M + g)
# }

# original definition
# lambda <- function(M, W, mu0 = -6.5, mu1 = 0.9, mu2 = 1) {
#   exp(mu0 + mu1 * W + mu2 * M)
# }
