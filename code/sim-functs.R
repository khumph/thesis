# simulation functions ----------------------------------------------------

# function for how toxicity changes
updateW <- function(M, W, D, a1 = 0.1, b1 = 1.1, d1 = 0.5) {
  W_next <- a1 * M + b1 * (D - d1) + W
  ifelse(W_next > 0, W_next, 0)
}

# function for updating tumor mass
updateM <- function(M, W, D, a2 = 0.15, b2 = 1.2, d2 = 0.5) {
  M_next <- ifelse(M > 0,
                   (a2 * W - b2 * (D - d2)) + M,
                   0)
  ifelse(M_next > 0, M_next, 0)
}

updateMint <- function(M, W, D, X, a2 = 0.15, b2 = 1.2, d2 = 0.5) {
  M_next <- ifelse(M > 0,
                   (a2 * W - b2 * (D * ifelse(X > 0.5, 2, 1) - d2)) + M,
                   0)
  ifelse(M_next > 0, M_next, 0)
}

# reward functions
R2 <- function(W_next, W) {
  -10 * (W_next - W)
}

R3 <- function(M_next, M) {
  ifelse(M_next == 0 & M == 0,
         0,
         ifelse(M_next == 0,
                30,
                -10 * (M_next - M)))
}

reward <- function(M_next, M, W_next, W, dead) {
  R2(W_next, W) +
    R3(M_next, M) +
    ifelse(dead, -60, 0)
}

reward_est <- function(M_next, M, W_next, W, pdied) {
  R2(W_next, W) +
    R3(M_next, M) +
    -60 * pdied
}

lambda <- function(W, M, mu0 = -7, mu1 = 1, mu2 = 1) {
  exp(mu0 + mu1 * W + mu2 * M)
}

pDeath <- function(M_next, W_next) {
  lam <- lambda(W_next, M_next)
  deltaF <- exp(-lam)
  1 - deltaF
}