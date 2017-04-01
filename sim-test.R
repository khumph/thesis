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
    pdeath = pexp(lam),
    surv_time = 1 / lam,
    reward = log(surv_time) #
    # reward = -M_next
  )
  if (!nested) {
    dat <- dat %>% group_by(ID) %>% mutate(
      bestR = max(reward),
      bestD = ifelse(near(reward, bestR), dose, NA),
      best = quantile(bestD, probs = 0, na.rm = T, type = 3, names = F)
      # best = ifelse(tumor_mass > 0,
      #               quantile(bestD, probs = 1, na.rm = T, type = 3, names = F),
      #               quantile(bestD, probs = 0, na.rm = T, type = 3, names = F))
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

# bp

max_best <- function(data, int, noise_pred) {
  data <- data %>% mutate(
    pdeath = NA
  )
  out <- NULL
  for (i in 5:1) {
    Q1 <- maxMonth(filter(data, month == i), int = int, noise_pred = noise_pred)
    d0 <- data %>% filter(month == i - 1) %>% ungroup() %>% mutate(
      reward = reward + ifelse(!is.na(Q1$dose), Q1$reward, 0)
    )
    out <- bind_rows(out, Q1)
  }
  i <- 0
  Q1 <- maxMonth(filter(data, month == i), int = int, noise_pred = noise_pred)
  bind_rows(out, Q1) %>% arrange(ID, month)
}


# test simulation -----------------------------------------------------

simMonthT <- function(dat, Q, int, noise_pred, pred = F) {
  out <- dat %>% filter(group != "best",
                        !str_detect(group, paste(names(Q), collapse = "|"))) %>%
    mutate(dose = ifelse(group == "1on1off",
                         ifelse(dat$month[1] %% 2 == 1, 0, 1),
                         dose))
  
  bestD <- maxMonth(filter(dat, group == "best"),
                    int = int, noise_pred = noise_pred, pred = pred)
  
  out <- bind_rows(out, bestD)

  for (i in seq_along(Q)) {
    if (!pred) {
      dat_opt <- filter(dat, !is.na(M_next), str_detect(group, names(Q)[i])) 
    } else {
      dat_opt <- filter(dat, str_detect(group, names(Q)[i])) 
    }
    
    optimD <- max_df(
      data = dat_opt,
      model = Q[[i]]$mod_list[[dat$month[1] + 1]],
      form = Q[[i]]$formula,
      mod_type = Q[[i]]$mod_type,
      pred = pred
    )
    bestDopt <- maxMonth(optimD, int = int, noise_pred = noise_pred, pred = pred)$dose
    optimD <- optimD %>% mutate(best_dose = bestDopt)
    
    if (!pred) {
      deads <- filter(dat, is.na(M_next), str_detect(group, names(Q)[i])) %>%
        mutate(dose = NA)
      optimD <- bind_rows(optimD, deads) %>% arrange(ID)
    }
    
    out <- out %>% bind_rows(optimD) %>% arrange(ID)
  }
  dat <- out 
  dat <- Mnext(dat, int, noise_pred)
  dat <- Wnext(dat, int, noise_pred)
  dat <- lamNext(dat, int, noise_pred)
  dat %>%
    mutate(
      pdeath = pexp(lam),
      surv_time = 1 / lam,
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
  out <- d
  for (i in 1:(Ttot - 1)) {
    d <- d %>% mutate(month = i,
                      tumor_mass = M_next,
                      toxicity = W_next) %>%
      simMonthT(Q, int = int, noise_pred = noise_pred) %>%
      mutate(reward = ifelse(!dead, 0, log(i + 1 + surv_time)))
      # mutate(reward = -M_next)
    out <- bind_rows(out, d)
  }
  out %>% group_by(ID) %>% mutate(
    reward = ifelse(month == Ttot - 1 & !dead, log(surv_time + Ttot - 1), reward),
    # reward = -M_next,
    pdeath = ifelse(is.na(pdeath), 1, pdeath),
    tot_reward = sum(reward[1:6], na.rm = T),
    Qhat = reward,
    best = NA
  ) %>% arrange(ID) %>% ungroup()
}



# recursive best dose -----------------------------------------------------

# best_df <- max_best(filter(out, group == "best"), int = int, noise_pred = noise_pred) %>%
#   mutate(group = "bestR") %>% arrange(ID)
# 
# bestDR <- best_df$dose
# 
# d <- simMonthT(dat, Q, int = int, noise_pred = noise_pred) %>% mutate(
#   reward = ifelse(dead, log(surv_time), 0)
# )
# out <- d
# for (i in 1:(Ttot - 1)) {
#   d <- d %>% mutate(month = i,
#                     tumor_mass = M_next,
#                     toxicity = W_next) %>%
#     simMonthT(Q, int = int, noise_pred = noise_pred) %>%
#     mutate(reward = ifelse(!dead, 0, log(i + 1 + surv_time)))
#   out <- bind_rows(out, d)
# }
# out <- 
#   out %>% group_by(ID) %>% mutate(
#     reward = ifelse(month == Ttot - 1 & !dead, log(surv_time + Ttot - 1), reward),
#     pdeath = ifelse(is.na(pdeath), 1, pdeath),
#     tot_reward = sum(reward, na.rm = T),
#     Qhat = reward,
#     best = NA
#   ) %>% arrange(ID)

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
