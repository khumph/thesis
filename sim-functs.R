# simulation functions ----------------------------------------------------

updateW <- function(M, W, D, a1 = 0.1, b1 = 1.2, c1, d1 = 0.5, truth) {
  truth <- rep(truth, length(M))
  W_next <-
    a1 * M + b1 * (c1 * D  - d1) + W #+ ifelse(!truth, rnorm(length(M), 0, 0.05), 0)
  ifelse(W_next > 0, W_next, 0)
}

Wnext <- function(dat, int, noise_pred, c1 = 1, truth = F) {
  if (int & truth) {
    dat %>% mutate(W_next = updateW(tumor_mass, toxicity, dose,
                                    c1 = ifelse(X1 > 0.5, 1.5, 1),
                                    truth = truth))
  } else if (int) {
    dat %>%
      mutate(W_next = ifelse(!dead,
                             updateW(tumor_mass, toxicity, dose,
                                     c1 = ifelse(X1 > 0.5, 1.5, 1),
                                     truth = truth),
                             NA))
  } else if (truth & !int) {
    dat %>% mutate(W_next = updateW(tumor_mass, toxicity, dose,
                              c1 = c1, truth = truth))
  } else {
    dat %>% mutate(W_next = ifelse(!dead, updateW(tumor_mass, toxicity, dose,
                                            c1 = c1, truth = truth), NA))
  }
}

updateM <- function(M, W, D, a2 = 0.15, b2 = 1.2, c2, d2 = 0.5,
                    X = 0, Z = 0, a3 = 1e-3, truth) {
  truth <- rep(truth, length(M))
  M_next <- ifelse(M > 0,
                   (a2 * W - b2 * (c2 * D  - d2)) +
                     M + sum(a3 * Z), #+ ifelse(!truth, rnorm(length(M), 0, 0.05), 0),
                   0)
  ifelse(M_next > 0, M_next, 0)
}

Mnext <- function(dat, int, noise_pred, c2 = 1, truth = F) {
  if (int & truth) {
    dat %>% mutate(M_next = updateM(tumor_mass, toxicity, dose, 
                              c2 = ifelse(X2 > 0.5, 1.5, 1), truth = truth))
  } else if (int) {
    dat %>%
      mutate(M_next = ifelse(!dead,
                             updateM(tumor_mass, toxicity, dose, 
                                     c2 = ifelse(X2 > 0.5, 1.5, 1), truth = truth),
                             NA))
  } else if (truth & !int) {
    dat %>% mutate(M_next = updateM(tumor_mass, toxicity, dose,
                                    c2 = c2, truth = truth))
  } else {
    dat %>% mutate(M_next = ifelse(!dead, updateM(tumor_mass, toxicity, dose, c2 = c2,
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
    dat %>%
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
  } else {
    dat
  }
}


lamNext <- function(dat, int, noise_pred) {
  if (noise_pred) {
    dat %>% mutate(lam = lambda(M_next, W_next,
                                Z = Z1 + Z2 + Z3 + Z4 + Z5 + Z6 + Z7 + Z8 + Z9 + Z10))
  } else {
    dat %>% mutate(lam = lambda(M_next, W_next))
  }
}
# defined as in NSCLC paper
lambda <- function(M, W,  mu0 = -7, mu1 = 1, mu2 = 1.2, mu3 = 1, a3 = 0.05, Z = 0) {
  exp(mu0 + mu1 * W + mu2 * M + mu3 * W * M + a3 * Z)
}

# our new definition
# lambda <- function(M_next, W_next, M, W, dose, mu0 = -6.5, mu1 = 1, mu2 = 1) {
#   dw <- (-0.1 * M - W) / 1.2 + 0.5
#   dm <- (-0.15 * W + M) / 1.2 + 0.5
#   mu0 + 3 * mu1 * exp(-(dose - dw)^2) +
#     3 * mu2 * exp(-(dose - dm)^2) - M_next - W_next
#   g <- (1)*(M - W)^2
#   exp(mu0 + mu1 * W + mu2 * M + g)
# }

# original definition
# lambda <- function(M, W, mu0 = -6.5, mu1 = 0.9, mu2 = 1) {
#   exp(mu0 + mu1 * W + mu2 * M)
# }
