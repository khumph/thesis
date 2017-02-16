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
    dead = ifelse(dead,
                  T,
                  ifelse(rbinom(
                    nrow(.), 1, pDeath(M_next, W_next)
                  ) == 1, T, F)),
    reward = reward(
      M_next = M_next,
      M = tumor_mass,
      W = toxicity,
      W_next = W_next,
      dead
    )
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
