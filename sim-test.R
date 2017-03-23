# best possible -----------------------------------------------------------

maxMonth <- function(dat, int, noise, nested = F) {
  deads <- filter(dat, dead) %>% mutate(dose = NA)
  dat <- filter(dat, !dead) %>% ungroup() %>% 
    mutate(dose = map(1:nrow(.), ~ seq(0, 1, by = 0.005))) %>%
    unnest()
  dat <- Mnext(dat, int, noise)
  dat <- Wnext(dat, int, noise)
  dat <- dat %>% mutate(
    lam = lambda(M_next, W_next),
    pdeath = pexp(lam),
    surv_time = 1 / lam,
    reward = log(surv_time)
  ) %>% group_by(ID) %>% mutate(
    bestR = max(reward),
    bestD = ifelse(near(reward, bestR), dose, NA),
    best = quantile(bestD,
                    probs = 0,
                    na.rm = T, type = 3, names = F)
    # best = ifelse(
    #   tumor_mass > 0,
    #   quantile(bestD, probs = 0, na.rm = T, type = 3, names = F),
    #   min(bestD, na.rm = T)
    # )
  ) 
  if (nested) {
    dat
  } else {
    dat <- dat %>% filter(near(dose, best))
    bind_rows(dat, deads) %>% arrange(ID)
  }
}

# test simulation -----------------------------------------------------

simMonthT <- function(dat, Q, int, noise) {
  optimD <- max_df(
    data = filter(dat, !dead, group == "optim"),
    model = Q$mod_list[[dat$month[1] + 1]],
    form = Q$formula,
    mod_type = Q$mod_type,
    pred = T
  ) %>% select(ID, dose)
  
  deads <- filter(dat, dead, group == "optim") %>% mutate(dose = NA)
  optimD <- bind_rows(optimD, deads) %>% arrange(ID)
  optimD <- optimD$dose
  
  bestD <- maxMonth(filter(dat, group == "best"),
                    int = int, noise = noise)$dose
  bestDopt <- maxMonth(filter(dat, group == "optim"),
                       int = int, noise = noise)$dose
  D1on1off <- ifelse(dat$month[1] %% 2 == 1, 0, 1)
  
  dat <- dat %>%
    mutate(
      dose = ifelse(group == "optim",
                    optimD,
                    ifelse(
                      group == "best",
                      bestD,
                      ifelse(group == "1on1off",
                             D1on1off, dose))),
      best_dose = ifelse(group == "optim", bestDopt, NA)
    )
  dat <- Mnext(dat, int, noise)
  dat <- Wnext(dat, int, noise)
  dat %>%
    mutate(
      lam = lambda(M_next, W_next),
      pdeath = pexp(lam),
      surv_time = 1 / lam,
      dead = ifelse(dead,
                    T,
                    surv_time < 1)
    )
}

sim_test <- function(Q, int, noise, npergroup = 200, ngroups = 13, Ttot = 6) {
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
      mod_type = Q$mod_type,
      pred = T
    )$best
  
  Dbest <- maxMonth(dat, int = int, noise = noise)$dose
  
  D1on1off <- rep(1, npergroup)
  
  groups <- c(
    seq(from = 0.1, to = 1, by = 0.1) %>% as.character(),
    "best",
    "optim", "1on1off")
  
  dat <- dat[rep(seq_len(nrow(dat)), ngroups), ] %>% 
    mutate(
      ID = rep(1:(npergroup * ngroups)),
      group = rep(groups, each = npergroup),
      dose = c(D1, Dbest, D0, D1on1off),
      best_dose = ifelse(group == "optim", Dbest, NA)
    )
  
  d <- simMonthT(dat, Q, int = int, noise = noise) %>% mutate(
    reward = ifelse(dead, log(surv_time), 0)
  )
  out <- d
  for (i in 1:(Ttot - 1)) {
    d <- d %>% mutate(month = i,
                      tumor_mass = M_next,
                      toxicity = W_next) %>%
      simMonthT(Q, int = int, noise = noise) %>%
      mutate(reward = ifelse(!dead, 0, log(i + 1 + surv_time)))
    out <- bind_rows(out, d)
  }
  out %>% group_by(ID) %>% mutate(
    reward = ifelse(month == Ttot - 1 & !dead, log(surv_time + Ttot - 1), reward),
    pdeath = ifelse(is.na(pdeath), 1, pdeath),
    tot_reward = sum(reward, na.rm = T),
    Qhat = reward,
    best = NA
  ) %>% arrange(ID)
}


# following estimated optimal regime --------------------------------------

# mon1 <- function(dat, int, noise) {
#   dat <- Mnext(dat, int, noise)
#   dat <- Wnext(dat, int, noise)
#   dat
# }
# 
# optDat <- function(dat, int, noise, Ttot = 6) {
#   m1 <- dat %>% group_by(ID) %>% filter(month == 0) %>% mon1(int, noise)
#   out <- m1
#   for (i in 2:Ttot) {
#     m1 <- m1 %>% mutate(tumor_mass = M_next,
#                         toxicity = W_next,
#                         dose = dat$best[i],
#                         month = i - 1) %>% mon1(int, noise)
#     out <- bind_rows(out, m1)
#   }
#   out
# }
