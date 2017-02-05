# simulation functions ----------------------------------------------------

# function for how toxicity changes
updateW <- function(M, W, D, a1 = 0.1, b1 = 1.2, d1 = 0.5) {
  W_next <- a1 * M + b1 * (D - d1) + W
  ifelse(W_next > 0, W_next, 0)
}

# function for updating tumor mass
updateM <- function(M, W, D, X, a2 = 0.15, b2 = 1.2, d2 = 0.5) {
  # takes current value of tumor mass, outputs next value
  M_next <- ifelse(M > 0,
                   (a2 * W - b2 * (D * ifelse(X > 0.5, 2, 1) - d2)) + M,
                   0)
  ifelse(M_next > 0, M_next, 0)
  # (a2 * W - b2 * (D - d2)) + M
}

# reward functions
R2 <- function(W_next, W) {
  -10 * (W_next - W)
}

R3 <- function(M_next, M) {
  ifelse(M_next == 0 & M == 0,
         0,
         ifelse(M_next == 0,
                15,
                -10 * (M_next - M)))
}

reward <- function(M_next, M, W_next, W, died) {
  R2(W_next, W) +
    R3(M_next, M) +
    ifelse(died == 1, -60, 0)
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

# sim ---------------------------------------------------------------------

simMonth <- function(dat) {
  dat %>%
    mutate(
      M_next = ifelse(!dead, updateM(tumor_mass, toxicity, dose, X), NA),
      W_next = ifelse(!dead, updateW(tumor_mass, toxicity, dose), NA),
      d_next = runif(nrow(.), min = 0, max = 1),
      died = rbinom(nrow(.), 1, pDeath(M_next, W_next)),
      dead = ifelse(dead,
                    T,
                    ifelse(died == 1, T, F)),
      reward = reward(
        M_next = M_next,
        M = tumor_mass,
        W = toxicity,
        W_next = W_next,
        died
      )
    )
}

sim <- function(N = 1000, Ttot = 6) {
  dat <- tibble(
    ID = 1:N,
    month = rep(0, N),
    tumor_mass = runif(N, min = 0, max = 2),
    toxicity = runif(N, min = 0, max = 2),
    X = runif(N, min = 0, max = 1),
    dose = runif(N, min = 0.5, max = 1),
    dead = rep(F, N)
  )
  d <- simMonth(dat)
  out <- d
  for (i in 1:(Ttot - 1)) {
    d <- d %>% mutate(month = i,
                      tumor_mass = M_next,
                      toxicity = W_next,
                      dose = d_next) %>% simMonth()
    out <- bind_rows(out, d)
  }
  d <- d %>% mutate(month = Ttot,
                    tumor_mass = M_next,
                    toxicity = W_next,
                    dose = NA, died = NA, dead = NA,
                    reward = NA)
  bind_rows(out, d) %>%
    select(ID, month, tumor_mass,
           toxicity, X, dose, died, dead, reward)
}

# best possible -----------------------------------------------------------

maxMonth <- function(dat, len = 101) {
  dat %>% ungroup() %>% 
    mutate(dose = map(1:nrow(.), ~ seq(0, 1, length.out = len))) %>%
    unnest() %>%
    mutate(
      M_next = updateM(tumor_mass, toxicity, dose, X),
      W_next = updateW(tumor_mass, toxicity, dose),
      pdied = pDeath(M_next, W_next),
      reward = reward_est(
        M_next = M_next,
        M = tumor_mass,
        W = toxicity,
        W_next = W_next,
        pdied
      )
    ) %>% group_by(ID) %>% mutate(
      max = max(reward),
      best = ifelse(near(reward, max), dose, NA),
      best = ifelse(
        tumor_mass > 0,
        max(best, na.rm = T),
        min(best, na.rm = T)
      )
    ) %>% filter(near(dose, best))
}

# sim test refactored -----------------------------------------------------

simMonthT <- function(dat) {
  optimD <- max_df(
    data = filter(dat, group == "optim"),
    model = Q$mod_list[[dat$month[1] + 1]],
    form = Q$formula,
    mod_type = Q$mod_type
  )$best
  
  bestD <- maxMonth(filter(dat, group == "best"))$dose
  
  dat %>%
    mutate(
      M_next = updateM(tumor_mass, toxicity, dose, X),
      W_next = updateW(tumor_mass, toxicity, dose),
      d_next = ifelse(group == "optim", optimD,
                      ifelse(group == "best", bestD, dose)), 
      pdeath = pDeath(M_next, W_next),
      reward = reward_est(M_next, tumor_mass, W_next, toxicity, pdeath)
    )
}

sim_test_df <- function(Q, npergroup = 200, ngroups = 12, Ttot = 6) {
  M0 <- runif(npergroup, min = 0, max = 2)
  W0 <- runif(npergroup, min = 0, max = 2)
  X <- runif(npergroup, min = 0, max = 1)

  dat <- tibble(
    ID = 1:npergroup,
    month = rep(0, npergroup),
    tumor_mass = M0,
    toxicity = W0,
    X = X,
    dead = rep(F, npergroup)
  )
  
  D1 <- rep(seq(from = 0.1, to = 1, by = 0.1), each = npergroup)
  
  D0 <-
    max_df(
      data = dat,
      model = Q$mod_list[[1]],
      form = Q$formula,
      mod_type = Q$mod_type
    )$best
  
  Dbest <- maxMonth(dat)$dose
  
  groups <- c(
    seq(from = 0.1, to = 1, by = 0.1) %>% as.character(),
    "best",
    "optim")
  
  dat <- tibble(
    ID = rep(1:(npergroup * ngroups)),
    group = rep(groups, each = npergroup),
    month = rep(0, npergroup * ngroups),
    tumor_mass = rep(M0, ngroups),
    toxicity = rep(W0, ngroups),
    X = rep(X, ngroups),
    dose = c(D1, Dbest, D0),
    dead = rep(F, npergroup * ngroups)
  )
  
  d <- simMonthT(dat)
  out <- d
  for (i in 1:(Ttot - 1)) {
    d <- d %>% mutate(month = i,
                      tumor_mass = M_next,
                      toxicity = W_next,
                      dose = d_next) %>% simMonthT()
    out <- bind_rows(out, d)
  }
  d <- d %>% mutate(month = Ttot,
                    tumor_mass = M_next,
                    toxicity = W_next,
                    dose = NA,
                    reward = NA)
  out <- bind_rows(out, d)
  out %>%
    select(ID, group, month, tumor_mass, toxicity, dose, pdeath, reward) %>%
    group_by(ID) %>%
    mutate(tot_reward = sum(reward, na.rm = T))
}


# results functions -------------------------------------------------------

align_df <- function(Q) {
  Q$data %>%
    select(ID,
           month,
           dose,
           best,
           reward,
           Q_hat,
           tumor_mass,
           toxicity,
           X,
           died, dead) %>%
    mutate(
      reward = lag(reward),
      Q_hat = lag(Q_hat),
      best = ifelse(died != 1 | is.na(died), best, NA)
    )
}

plots_tab <- function(dat_test_long) {
  dat_long_summ <- dat_test_long %>% group_by(group, month) %>%
    summarise(
      mean_tox = mean(toxicity, na.rm = T),
      mean_tumor = mean(tumor_mass, na.rm = T),
      mean_reward = mean(reward, na.rm = T)
    ) %>% mutate(sum_means = mean_tox + mean_tumor)
  
  plot_tox <- ggplot(data = dat_long_summ) +
    geom_line(mapping = aes(
      x = month,
      y = mean_tox,
      color = group,
      group = group
    ))
  
  plot_tumor <- ggplot(data = dat_long_summ) +
    geom_line(mapping = aes(
      x = month,
      y = mean_tumor,
      color = group,
      group = group
    ))
  
  plot_sum <- ggplot(data = dat_long_summ) +
    geom_line(mapping = aes(
      x = month,
      y = sum_means,
      color = group,
      group = group
    ))
  
  plot_reward <- ggplot(data = dat_long_summ) +
    geom_line(mapping = aes(
      x = month,
      y = mean_reward,
      color = group,
      group = group
    ))
  
  tab_deaths <-
    dat_test_long %>% group_by(group) %>%
    summarise(pdeath = sum(pdeath, na.rm = T) / n()) %>%
    arrange(pdeath)
  
  tab_reward <-
    dat_test_long %>% select(ID, group, tot_reward) %>% unique() %>%
    group_by(group) %>%
    summarise(avg_tot_reward = mean(tot_reward)) %>%
    arrange(desc(avg_tot_reward))
  
  list(
    plot_tox = plot_tox,
    plot_tumor = plot_tumor,
    plot_sum = plot_sum,
    plot_reward = plot_reward,
    table_deaths = tab_deaths,
    table_rewards = tab_reward
  )
}
