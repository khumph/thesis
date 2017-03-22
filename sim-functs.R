# simulation functions ----------------------------------------------------

# function for how toxicity changes
updateW <- function(M, W, D, a1 = 0.1, b1 = 1.2, d1 = 0.5,
                    Z = 0, a3 = 1e-3) {
  W_next <- a1 * M + b1 * (D - d1) + W + sum(a3 * Z)
  ifelse(W_next > 0, W_next, 0)
}

Wnext <- function(dat, int, noise) {
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
updateM <- function(M, W, D, a2 = 0.15, b2 = 1.2, d2,
                    X = 0, Z = 0, a3 = 1e-3) {
  M_next <- ifelse(M > 0,
                   (a2 * W - b2 * (D  - d2)) +
                     M + sum(a3 * Z),
                   0)
  ifelse(M_next > 0, M_next, 0)
}

# , d2 = (5 + dat$month[1])/10
Mnext <- function(dat, int, noise, d2 = 0.5) {
  if (int) {
    dat %>%
      mutate(M_next = ifelse(!dead,
                             updateM(tumor_mass, toxicity, dose, d2 = X),
                             NA))
  } else if (noise) {
    dat %>%
      mutate(M_next = ifelse(!dead,
                             updateM(tumor_mass, toxicity, dose, d2 = d2,
                                     Z = c(Z1, Z2, Z3, Z4, Z5, Z6)),
                             NA))
  } else {
    dat %>%
      mutate(M_next = ifelse(!dead, updateM(tumor_mass, toxicity, dose, d2 = d2), NA))
  }
}

genIntNoise <- function(dat, int, noise) {
  if (int) {
    dat %>%
      mutate(X = runif(nrow(.), min = 0, max = 1),
             X = ifelse(X > 0.5, 0.8, 0.5))
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
lambda <- function(M, W, mu0 = -5, mu1 = 0.9, mu2 = 1) {
  exp(mu0 + mu1 * W + mu2 * M + 0.5 * W * M)
}

# our new definition
# lambda <- function(M, W, mu0 = -6.5, mu1 = 1, mu2 = 1) {
#   g <- (1)*(M - W)^2
#   exp(mu0 + mu1 * W + mu2 * M + g)
# }

# original definition
# lambda <- function(M, W, mu0 = -6.5, mu1 = 1, mu2 = 1) {
#   exp(mu0 + mu1 * W + mu2 * M)
# }
