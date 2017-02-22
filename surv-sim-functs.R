# simulation functions ----------------------------------------------------

# function for how toxicity changes
updateW <- function(M, W, D, a1 = 0.1, b1 = 1.2, d1 = 0.5) {
  W_next <- a1 * M + b1 * (D - d1) + W
  ifelse(W_next > 0, W_next, 0)
}

# function for updating tumor mass
updateM <- function(M, W, D, a2 = 0.15, b2 = 1.2, d2 = 0.5,
                    X = 0, V = 0, a3 = 1e-4) {
  M_next <- ifelse(M > 0,
                   (a2 * W - b2 * (D * ifelse(X > 0.5, 2, 1) - d2)) +
                     M + sum(a3 * V),
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
                                     V = c(V1, V2, V3, V4, V5,
                                           V6, V7, V8, V9, V10)),
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
    dat %>%
      bind_cols(replicate(100, rnorm(nrow(.), mean = -0.05, sd = 1)) %>% as.data.frame())
  } else {
    dat
  }
}




reward <- function(M_next, M, W_next, W, dead) {
  ifelse(!is.na(M_next), ifelse(dead, -2, 0), NA)
}

reward_est <- function(M_next, M, W_next, W, pdied) {
    1 - pdied
}

lambda <- function(W, M, mu0 = -7, mu1 = 1, mu2 = 1) {
  exp(mu0 + mu1 * W + mu2 * M)
}
