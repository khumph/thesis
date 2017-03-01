# simulation functions ----------------------------------------------------

# function for how toxicity changes
updateW <- function(M, W, D, a1 = 0.1, b1 = 1.2, d1 = 0.5,
                    Z = 0, a3 = 1e-3) {
  W_next <- a1 * M + b1 * (D - d1) + W + sum(a3 * Z)
  ifelse(W_next > 0, W_next, 0)
}

Wnext <- function(dat, int = F, noise = F) {
  if (noise) {
    dat %>%
      mutate(W_next = ifelse(!dead,
                             updateW(tumor_mass, toxicity, dose,
                                     Z = c(Z7, Z8, Z9, Z10, Z11, Z12)),
                             NA))
  } else {
    dat %>%
      mutate(W_next = ifelse(!dead, updateW(tumor_mass, toxicity, dose), NA))
  }
}

# function for updating tumor mass
updateM <- function(M, W, D, a2 = 0.15, b2 = 1.2, d2 = 0.5,
                    X = 0, Z = 0, a3 = 1e-3) {
  M_next <- ifelse(M > 0,
                   (a2 * W - b2 * (D * ifelse(X > 0.5, 2, 1) - d2)) +
                     M + sum(a3 * Z),
                   0)
  ifelse(M_next > 0, M_next, 0)
}

Mnext <- function(dat, int = F, noise = F) {
  if (int) {
    dat %>%
      mutate(M_next = ifelse(!dead,
                             updateM(tumor_mass, toxicity, dose, X = X),
                             NA))
  } else if (noise) {
    dat %>%
      mutate(M_next = ifelse(!dead,
                             updateM(tumor_mass, toxicity, dose,
                                     Z = c(Z1, Z2, Z3, Z4, Z5, Z6)),
                             NA))
  } else {
    dat %>%
      mutate(M_next = ifelse(!dead, updateM(tumor_mass, toxicity, dose), NA))
  }
}

genIntNoise <- function(dat, int, noise) {
  if (int) {
    dat %>%
      mutate(X = runif(N, min = 0, max = 1))
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

lambda <- function(M_next, W_next, M, W, dose, mu0 = -7, mu1 = 1, mu2 = 1) {
  diff <- (M - W)/(M + W)
  g <- ifelse(M > W,
             (dose - diff) ^ 2,
             ifelse(M == 0 & W == 0,
                    dose ^ 2,
                    (dose + diff) ^ 2))
  g <- 0
  exp(mu0 + mu1 * W_next + mu2 * M_next + g)
}
