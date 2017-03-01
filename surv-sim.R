# sim ---------------------------------------------------------------------

simMonth <- function(dat, int, noise) {
  dat <- Mnext(dat, int, noise)
  dat <- Wnext(dat, int, noise)
  dat %>% mutate(
    d_next = runif(nrow(.), min = 0, max = 1),
    surv_time = rexp(nrow(.), lambda(M_next, W_next, tumor_mass, toxicity, dose)),
    dead = ifelse(dead,
                  T,
                  surv_time < 1)
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
  
  dat <- genIntNoise(dat, int, noise)
  
  d <- simMonth(dat, int = int, noise = noise) %>% mutate(
    reward = ifelse(dead, log(surv_time), log(1))
  )
  out <- d
  for (i in 1:(Ttot - 1)) {
    d <- d %>% mutate(month = i,
                      tumor_mass = M_next,
                      toxicity = W_next,
                      dose = d_next) %>%
      simMonth(int = int, noise = noise) %>%
      mutate(reward = ifelse(!dead, log(i + 1), log(i + 1 + surv_time)))
    out <- bind_rows(out, d)
  }
  out %>% mutate(
    reward = ifelse(month == Ttot - 1 & !dead, log(surv_time + Ttot - 1), reward),
    Qhat = reward,
    best = NA
  )
}
