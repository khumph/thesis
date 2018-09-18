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

genIntNoise <- function(dat, int, noise, noise_pred) {

  n <- nrow(dat)

  dat$cW <- rep(1, n)
  dat$cM <- rep(1, n)

  if (int) {
    dat$X1 <- runif(n, min = 0, max = 1)
    dat$X2 <- runif(n, min = 0, max = 1)
    dat$cW <- replace(dat$cW, dat$X1 > 0.5, 1.5)
    dat$cM <- replace(dat$cM, dat$X2 > 0.5, 1.5)
  }

  if (noise) {
    Z1 <- lapply(1:5, function(x) rnorm(n, mean = 1))
    Z2 <- lapply(1:5, function(x) rnorm(n, mean = -1))
    V <- lapply(1:90, function(x) rnorm(n))

    names(Z1) <- paste0("Z", 1:5)
    names(Z2) <- paste0("Z", 6:10)
    names(V) <- paste0("V", 1:90)

    dat <- bind_cols(dat, Z1, Z2, V)
  }

  if (noise_pred) {
    dat$noise_chng <- dat$Z1 + dat$Z2 + dat$Z3 + dat$Z4 + dat$Z5 +
      dat$Z6 + dat$Z7 + dat$Z8 + dat$Z9 + dat$Z10
  } else {
    dat$noise_chng <- rep(0, n)
  }
  return(dat)
}


lambda <- function(M, W, Z, mu0 = -5.5, mu1 = 1, mu2 = 1.2, mu3 = 0.75, a3 = 0.05) {
  exp(mu0 + mu1 * W + mu2 * M + mu3 * W * M + a3 * Z)
}


simMonth <- function(dat) {
  dat <- Mnext(dat)
  dat <- Wnext(dat)
  dat$lambda <- lambda(dat$M_next, dat$W_next, Z = dat$noise_chng)
  dat$d_next <- runif(nrow(dat), min = 0, max = 1)
  dat$surv_time <- rexp(nrow(dat), dat$lambda)
  dat$dead <- replace(dat$dead, dat$surv_time < 1, T)
  return(dat)
}

sim <- function(n_subjects = 1000, n_stages = 3, int = F, noise = F,
                noise_pred = F, seed = 1) {

    set.seed(seed)

    dat <- data.frame(
      ID = 1:n_subjects,
      month = rep(0, n_subjects),
      tumor_mass = runif(n_subjects, min = 0, max = 2),
      toxicity = runif(n_subjects, min = 0, max = 2),
      dose = runif(n_subjects, min = 0, max = 1),
      dead = rep(F, n_subjects)
    )

    dat <- genIntNoise(dat, int, noise, noise_pred)

    dat <- simMonth(dat)
    dat$reward <- replace(log(dat$surv_time), !(dat$dead), 0)

    out <- dat
    for (i in 1:(n_stages - 1)) {
      dat$month <- rep(i, nrow(dat))
      dat$tumor_mass <- dat$M_next
      dat$toxicity <- dat$W_next
      dat$dose <- dat$d_next
      dat <- simMonth(dat)
      dat$reward <- replace(log(i + 1 + dat$surv_time), !dat$dead, 0)
      out <- bind_rows(out, dat)
    }

    out$reward <- ifelse(out$month == n_stages - 1 & !out$dead,
                         log(out$surv_time + n_stages - 1), out$reward)
    out$Qhat <- out$reward

    return(out)
}
