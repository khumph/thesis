# best possible -----------------------------------------------------------

maxMonth <- function(dat, len = 101, int, noise) {
  dat <- dat %>% ungroup() %>% 
    mutate(dose = map(1:nrow(.), ~ seq(0, 1, by = 0.05 ))) %>%
    unnest()
  dat <- Mnext(dat, int, noise)
  dat <- dat %>% mutate(
    W_next = updateW(tumor_mass, toxicity, dose),
    lam = lambda(W_next, M_next),
    pdeath = pexp(lam),
    surv_time = 1 / lam,
    reward = log(surv_time)
  ) %>% group_by(ID) %>% mutate(
    bestR = max(reward),
    bestD = ifelse(near(reward, bestR), dose, NA),
    best = ifelse(
      tumor_mass > 0,
      quantile(bestD, probs = 1, na.rm = T, type = 3, names = F),
      min(bestD, na.rm = T)
    )
  ) 
  dat <- dat %>% filter(near(dose, best))
}

# test simulation -----------------------------------------------------

simMonthT <- function(dat, Q, int, noise) {
  optimD <- max_df(
    data = filter(dat, group == "optim"),
    model = Q$mod_list[[dat$month[1] + 1]],
    form = Q$formula,
    mod_type = Q$mod_type
  )$best
  
  bestD <- maxMonth(filter(dat, group == "best"),
                    int = int, noise = noise)$dose
  bestDopt <- maxMonth(filter(dat, group == "optim"),
                       int = int, noise = noise)$dose
  
  dat <- Mnext(dat, int, noise)
  
  dat %>%
    mutate(
      W_next = updateW(tumor_mass, toxicity, dose),
      dose = ifelse(group == "optim", optimD,
                    ifelse(group == "best", bestD, dose)), 
      best_dose = ifelse(group == "optim", bestDopt, dose),
      lam = lambda(W_next, M_next),
      pdeath = pexp(lam),
      surv_time = 1 / lam,
      reward = log(surv_time)
    )
}

sim_test <- function(Q, int, noise, npergroup = 200, ngroups = 12, Ttot = 6) {
  M0 <- runif(npergroup, min = 0, max = 2)
  W0 <- runif(npergroup, min = 0, max = 2)
  
  dat <- tibble(
    ID = 1:npergroup,
    month = rep(0, npergroup),
    tumor_mass = M0,
    toxicity = W0,
    dead = rep(F, npergroup)
  )
  
  dat <- genIntNoise(dat, int, noise)
  
  D1 <- rep(seq(from = 0.1, to = 1, by = 0.1), each = npergroup)
  
  D0 <-
    max_df(
      data = dat,
      model = Q$mod_list[[1]],
      form = Q$formula,
      mod_type = Q$mod_type
    )$best
  
  Dbest <- maxMonth(dat, int = int, noise = noise)$dose
  
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
  
  d <- simMonthT(dat, Q, int = int, noise = noise)
  out <- d
  for (i in 1:(Ttot - 1)) {
    d <- d %>% mutate(month = i,
                      tumor_mass = M_next,
                      toxicity = W_next) %>%
      simMonthT(Q, int = int, noise = noise)
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


# following estimated optimal regime --------------------------------------

mon1 <- function(dat, int, noise) {
  dat <- Mnext(dat, int, noise)
  dat %>% mutate(
    W_next = ifelse(!dead, updateW(tumor_mass, toxicity, dose), NA)
  )
}

optDat <- function(dat, int, noise, Ttot = 6) {
  m1 <- dat %>% group_by(ID) %>% filter(month == 0) %>% mon1(int, noise)
  out <- m1
  for (i in 2:Ttot) {
    m1 <- m1 %>% mutate(tumor_mass = M_next,
                        toxicity = W_next,
                        dose = dat$best[i],
                        month = i - 1) %>% mon1(int, noise)
    out <- bind_rows(out, m1)
  }
  out
}
