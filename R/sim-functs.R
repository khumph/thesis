updateW <- function(M, W, D, a = 0.1, b = 1.2, c, d = 0.5) {
  W_next <- a * M + b * (c * D  - d) + W
  return(replace(W_next, W_next < 0, 0))
}

Wnext <- function(dat) {
  dat$W_next <- updateW(dat$tumor_mass, dat$toxicity, dat$dose, c = dat$cW)
  dat$W_next <- replace(dat$W_next, dat$dead, NA_real_)
  return(dat)
}

updateM <- function(M, W, D, a = 0.15, b = 1.2, c, d = 0.5, Z = 0, g = 1e-3) {
  M_next <- a * W - b * (c * D - d) + M + sum(g * Z)
  return(replace(M_next, M <= 0 | M_next < 0, 0))
}

Mnext <- function(dat) {
  dat$M_next <- updateM(dat$tumor_mass, dat$toxicity, dat$dose, c = dat$cM)
  dat$M_next <- replace(dat$M_next, dat$dead, NA_real_)
  return(dat)
}

lambda <- function(M, W, Z, mu0 = -5.5, mu1 = 1, mu2 = 1.2, mu3 = 0.75, a3 = 0.05) {
  exp(mu0 + mu1 * W + mu2 * M + mu3 * W * M + a3 * Z)
}

