# test simulation -----------------------------------------------------

getPreds <- function(Q, name, dat, pred) {
  optimD <- max_df(
    dat = filter(dat, str_detect(group, paste0("^", name, "$"))),
    model = Q$mod_list[[dat$month[1] + 1]],
    truth = F, pred = pred)
  
  bestDopt <- max_df(dat = optimD, model = NULL,
                     truth = T, pred = pred)$dose
  optimD %>% mutate(best_dose = bestDopt)
}

simMonthT <- function(dat, Q, pred) {
  out <- dat %>% filter(group != "best",
                        !str_detect(group, paste(names(Q), collapse = "|"))) %>%
    mutate(dose = ifelse(group == "1on1off",
                         ifelse(dat$month[1] %% 2 == 1, 0, 1),
                         dose))
  
  bestD <- max_df(dat = filter(dat, group == "best"),
                  model = NULL, truth = T, pred = pred)
  
  optimD <- map2_df(Q, names(Q),
                    ~ getPreds(.x, .y, dat = dat, pred = pred))
  
  dat <- bind_rows(out, bestD, optimD) %>% arrange(ID) 
  dat <- Mnext(dat)
  dat <- Wnext(dat)
  dat %>%
    mutate(
      lam = lambda(M_next, W_next, Z = noise_chng),
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
  
  dat <- genIntNoise(dat, int, noise, noise_pred)
  
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
  
  d <- simMonthT(dat, Q, pred = T) %>%
    mutate(reward = ifelse(dead, log(surv_time), 0))
  out <- d
  for (i in 1:(Ttot - 1)) {
    d <- d %>% mutate(month = i,
                      tumor_mass = M_next,
                      toxicity = W_next) %>%
      simMonthT(Q, pred = F) %>%
      mutate(reward = ifelse(!dead, 0, log(i + 1 + surv_time)))
    out <- bind_rows(out, d)
  }
  out %>% group_by(ID) %>% mutate(
    reward = ifelse(month == Ttot - 1 & !dead, log(surv_time + Ttot - 1), reward),
    pdeath = ifelse(is.na(pdeath), 1, pdeath),
    tot_reward = sum(reward[1:6], na.rm = T),
    Qhat = reward,
    best = NA
  ) %>% arrange(ID) %>% ungroup()
}
