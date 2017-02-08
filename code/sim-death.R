# simulation functions ----------------------------------------------------

# function for how toxicity changes
updateW <- function(M, W, D, a1 = 0.1, b1 = 1.2, d1 = 0.5) {
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
  # takes current value of tumor mass, outputs next value
  M_next <- ifelse(M > 0,
                   (a2 * W - b2 * (D * ifelse(X > 0.5, 2, 1) - d2)) + M,
                   0)
  ifelse(M_next > 0, M_next, 0)
}

lambda <- function(W, M, mu0 = -6, mu1 = 1, mu2 = 1) {
  exp(mu0 + mu1 * W + mu2 * M)
}

pDeath <- function(M_next, W_next) {
  lam <- lambda(W_next, M_next)
  deltaF <- exp(-lam)
  1 - deltaF
}

# sim ---------------------------------------------------------------------

simMonth <- function(dat, int) {
  if (int) {
    dat <- dat %>%
      mutate(M_next = ifelse(!dead, updateMint(tumor_mass, toxicity, dose, X), NA))
  } else {
    dat <- dat %>%
      mutate(M_next = ifelse(!dead, updateM(tumor_mass, toxicity, dose), NA))
  }
  dat %>% mutate(
    W_next = ifelse(!dead, updateW(tumor_mass, toxicity, dose), NA),
    d_next = runif(nrow(.), min = 0, max = 1),
    died = ifelse(!dead, rbinom(nrow(.), 1, pDeath(M_next, W_next)) == 1, F),
    dead = ifelse(dead, T, died),
    reward = ifelse(died, -60, 0) + ifelse(M_next == 0 & tumor_mass == 0, 0,
                                             ifelse(M_next == 0, 15, 0))
  )
}

sim <- function(N = 1000, Ttot = 6, int = F, noise = F) {
  dat <- tibble(
    ID = 1:N,
    month = rep(0, N),
    tumor_mass = runif(N, min = 0, max = 2),
    toxicity = runif(N, min = 0, max = 2),
    dose = runif(N, min = 0, max = 1),
    dead = rep(F, N)
  )
  if (int) {
    dat <- dat %>%
      mutate(X = runif(N, min = 0, max = 1))
  } else if (noise) {
    dat <- dat %>%
      bind_cols(replicate(10, runif(N, min = 0, max = 1)) %>% as.data.frame())
  }
  d <- simMonth(dat, int = int)
  out <- d
  for (i in 1:(Ttot - 1)) {
    d <- d %>% mutate(month = i,
                      tumor_mass = M_next,
                      toxicity = W_next,
                      dose = d_next) %>% simMonth(int = int)
    out <- bind_rows(out, d)
  }
  d <- d %>% mutate(month = Ttot,
                    tumor_mass = M_next,
                    toxicity = W_next,
                    dose = NA, dead = NA, reward = NA)
  bind_rows(out, d) %>% mutate(
    Q_hat = ifelse(month == 6, NA,
                   ifelse(month == 5, reward, 0)),
    best = NA
  )
}

# best possible -----------------------------------------------------------

maxMonth <- function(dat, len = 101, int) {
  dat <- dat %>% ungroup() %>% 
    mutate(dose = map(1:nrow(.), ~ seq(0, 1, length.out = len))) %>%
    unnest()
  if (int) {
    dat <- dat %>%
      mutate(M_next = updateMint(tumor_mass, toxicity, dose, X))
  } else {
    dat <- dat %>%
      mutate(M_next = updateM(tumor_mass, toxicity, dose))
  }
  dat %>% mutate(
    W_next = updateW(tumor_mass, toxicity, dose),
    pdied = pDeath(M_next, W_next),
    reward = -60 * pdied + ifelse(M_next == 0 & tumor_mass == 0, 0,
                                  ifelse(M_next > 0, 0, 15))
  ) %>% group_by(ID) %>% mutate(
    max = max(reward),
    best = ifelse(near(reward, max), dose, NA),
    best = ifelse(
      M_next > 0,
      max(best, na.rm = T),
      min(best, na.rm = T)
    )
  ) %>% filter(near(dose, best))
}

# sim test -----------------------------------------------------

simMonthT <- function(dat, Q, int) {
  optimD <- max_df(
    data = filter(dat, group == "optim"),
    model = Q$mod_list[[dat$month[1] + 1]],
    form = Q$formula,
    mod_type = Q$mod_type
  )$best
  
  bestD <- maxMonth(filter(dat, group == "best"), int = int)$dose
  
  if (int) {
    dat <- dat %>%
      mutate(
        M_next = ifelse(!dead, updateMint(tumor_mass, toxicity, dose, X), NA)
      )
  } else {
    dat <- dat %>%
      mutate(M_next = ifelse(!dead, updateM(tumor_mass, toxicity, dose), NA))
  }
  
  dat %>%
    mutate(
      W_next = updateW(tumor_mass, toxicity, dose),
      dose = ifelse(group == "optim", optimD,
                      ifelse(group == "best", bestD, dose)), 
      pdeath = pDeath(M_next, W_next),
      reward = -60 * pdeath + ifelse(M_next == 0 & tumor_mass == 0, 0,
                                       ifelse(M_next > 0, 0, 15))
    )
}

sim_test <- function(Q, npergroup = 200, ngroups = 12, Ttot = 6,
                     int = F, noise = F) {
  M0 <- runif(npergroup, min = 0, max = 2)
  W0 <- runif(npergroup, min = 0, max = 2)
  
  dat <- tibble(
    ID = 1:npergroup,
    month = rep(0, npergroup),
    tumor_mass = M0,
    toxicity = W0,
    dead = rep(F, npergroup)
  )
  
  if (int) {
    X <- runif(npergroup, min = 0, max = 1)
    dat <- dat %>%
      mutate(X = X)
  } else if (noise) {
    V <- replicate(10, runif(npergroup, min = 0, max = 1))
    dat <- dat %>% bind_cols(V %>% as.data.frame())
  }
  
  D1 <- rep(seq(from = 0.1, to = 1, by = 0.1), each = npergroup)
  
  D0 <-
    max_df(
      data = dat,
      model = Q$mod_list[[1]],
      form = Q$formula,
      mod_type = Q$mod_type
    )$best

  Dbest <- maxMonth(dat, int = int)$dose
  
  groups <- c(
    seq(from = 0.1, to = 1, by = 0.1) %>% as.character(),
    "best",
    "optim")
  
  dat <- dat[rep(seq_len(nrow(dat)), ngroups), ] %>% 
    mutate(
      ID = rep(1:(npergroup * ngroups)),
      group = rep(groups, each = npergroup),
      dose = c(D1, Dbest, D0)
    )
  
  d <- simMonthT(dat, Q, int = int)
  out <- d
  for (i in 1:(Ttot - 1)) {
    d <- d %>% mutate(month = i,
                      tumor_mass = M_next,
                      toxicity = W_next) %>% simMonthT(Q, int = int)
    out <- bind_rows(out, d)
  }
  d <- d %>% mutate(month = Ttot,
                    tumor_mass = M_next,
                    toxicity = W_next,
                    dose = NA,
                    reward = NA)
  bind_rows(out, d)
}