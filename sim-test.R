# best possible -----------------------------------------------------------

maxMonth <- function(dat, len = 101, int) {
  dat <- dat %>% ungroup() %>% 
    mutate(dose = map(1:nrow(.), ~ seq(0, 1, length.out = len))) %>%
    unnest()
  if (int) {
    dat <- dat %>%
      mutate(M_next = updateM(tumor_mass, toxicity, dose, X = X))
  } else {
    dat <- dat %>%
      mutate(M_next = updateM(tumor_mass, toxicity, dose))
  }
  dat %>% mutate(
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
      # (nnet::which.is.max(max) - 1) / 100,
      quantile(best, probs = 0, na.rm = T, type = 3, names = F),
      min(best, na.rm = T)
    )
  ) %>% filter(near(dose, best))
}

# test simulation -----------------------------------------------------

simMonthT <- function(dat, Q, int) {
  optimD <- max_df(
    data = filter(dat, group == "optim"),
    model = Q$mod_list[[dat$month[1] + 1]],
    form = Q$formula,
    mod_type = Q$mod_type
  )$best
  
  bestD <- maxMonth(filter(dat, group == "best"), int = int, len = length(optimD))$dose
  bestDopt <- maxMonth(filter(dat, group == "optim"), int = int, len = length(optimD))$dose
  
  if (int) {
    dat <- dat %>%
      mutate(
        M_next = ifelse(!dead, updateM(tumor_mass, toxicity, dose, X = X), NA)
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
      best_dose = ifelse(group == "optim", bestDopt, dose),
      pdeath = pDeath(M_next, W_next),
      reward = reward_est(M_next, tumor_mass, W_next, toxicity, pdeath)
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
    V <- replicate(100, runif(npergroup, min = -0.5, max = 1))
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
      dose = c(D1, Dbest, D0),
      best_dose = ifelse(group == "optim", Dbest, NA)
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
                    best_dose = NA,
                    reward = NA)
  bind_rows(out, d)
}