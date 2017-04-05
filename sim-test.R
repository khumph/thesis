# best possible -----------------------------------------------------------

maxMonth <- function(dat, int, noise_pred, nested = F, pred = F) {
  dat <- ungroup(dat)
  if (!pred) {
    deads <- filter(dat, is.na(M_next)) %>% mutate(dose = NA)
    dat <- filter(dat, !is.na(M_next))
  }
  dat <- dat %>% 
    mutate(dose = map(1:nrow(.), ~ seq(0, 1, by = 0.005))) %>%
    unnest()
  dat <- Mnext(dat, int, noise_pred, truth = T)
  dat <- Wnext(dat, int, noise_pred, truth = T)
  dat <- lamNext(dat, int, noise_pred)
  dat <- dat %>% mutate(
    expect_surv_time = 1 / lam,
    reward = log(expect_surv_time)
  )
  if (!nested) {
    dat <- dat %>% group_by(ID) %>% mutate(
      bestR = max(reward),
      bestD = ifelse(near(reward, bestR), dose, NA),
      best = quantile(bestD, probs = 0, na.rm = T, type = 3, names = F)
    )
  }
  if (nested) {
    dat %>% ungroup()
  } else if (pred) {
    dat %>% filter(near(dose, best)) %>% arrange(ID) %>% ungroup()
  } else {
    dat %>% filter(near(dose, best)) %>% bind_rows(deads) %>%
      arrange(ID) %>% ungroup()
  }
}


# test simulation -----------------------------------------------------

getPreds <- function(Q, name, dat, pred, int, noise_pred) {
  regex_name <- paste0("^", name, "$")
  if (!pred) {
    dat_opt <- filter(dat, !is.na(M_next), str_detect(group, regex_name))
  } else {
    dat_opt <- filter(dat, str_detect(group, regex_name)) 
  }
  
  optimD <- max_df(
    data = dat_opt,
    model = Q$mod_list[[dat$month[1] + 1]],
    form = Q$formula,
    mod_type = Q$mod_type,
    pred = pred
  )
  bestDopt <- maxMonth(optimD, int = int, noise_pred = noise_pred, pred = pred)$dose
  optimD <- optimD %>% mutate(best_dose = bestDopt)
  
  if (!pred) {
    deads <- filter(dat, is.na(M_next), str_detect(group, regex_name)) %>%
      mutate(dose = NA)
    bind_rows(optimD, deads) %>% arrange(ID)
  } else {
    optimD
  }
}

simMonthT <- function(dat, Q, int, noise_pred, pred = F) {
  out <- dat %>% filter(group != "best",
                        !str_detect(group, paste(names(Q), collapse = "|"))) %>%
    mutate(dose = ifelse(group == "1on1off",
                         ifelse(dat$month[1] %% 2 == 1, 0, 1),
                         dose))
  
  bestD <- maxMonth(filter(dat, group == "best"),
                    int = int, noise_pred = noise_pred, pred = pred)
  
  optimD <- map2_df(Q, names(Q),
                    ~ getPreds(.x, .y, dat = dat, pred = pred,
                               int = int, noise_pred = noise_pred))
  
  dat <- bind_rows(out, bestD, optimD) %>% arrange(ID) 
  dat <- Mnext(dat, int, noise_pred)
  dat <- Wnext(dat, int, noise_pred)
  dat <- lamNext(dat, int, noise_pred)
  dat %>%
    mutate(
      pdeath = pexp(lam),
      surv_time = rexp(nrow(.), lam),
      expect_surv_time = 1 / lam,
      dead = ifelse(dead, T, surv_time < 1)
    )
}

sim_test <- function(Q, int, noise, noise_pred, npergroup = 200, Ttot = 6) {
  ngroups <- 12 + length(Q)
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
  
  groups <- c(seq(from = 0.1, to = 1, by = 0.1) %>% as.character(),
    "best",
    names(Q), "1on1off")
  
  dat <- dat[rep(seq_len(nrow(dat)), ngroups), ] %>% 
    mutate(
      ID = rep(1:(npergroup * ngroups)),
      group = rep(groups, each = npergroup),
      dose = c(D1, rep(NA, npergroup * (ngroups - 10))),
      best_dose = NA
    )
  
  d <- simMonthT(dat, Q, int = int, noise_pred = noise_pred, pred = T) %>%
    mutate(reward = ifelse(dead, log(surv_time), 0))
    # mutate(reward = -M_next)
    # mutate(reward = log(expect_surv_time))
  out <- d
  for (i in 1:(Ttot - 1)) {
    d <- d %>% mutate(month = i,
                      tumor_mass = M_next,
                      toxicity = W_next) %>%
      simMonthT(Q, int = int, noise_pred = noise_pred) %>%
      mutate(reward = ifelse(!dead, 0, log(i + 1 + surv_time)))
      # mutate(reward = -M_next)
      # mutate(reward = log(expect_surv_time))
    out <- bind_rows(out, d)
  }
  out %>% group_by(ID) %>% mutate(
    reward = ifelse(month == Ttot - 1 & !dead, log(surv_time + Ttot - 1), reward),
    # reward = -M_next,
    # reward = log(expect_surv_time),
    pdeath = ifelse(is.na(pdeath), 1, pdeath),
    tot_reward = sum(reward[1:6], na.rm = T),
    Qhat = reward,
    best = NA
  ) %>% arrange(ID) %>% ungroup()
}
