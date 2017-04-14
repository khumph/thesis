
getPreds <- function(Q, name, dat, pred) {
  optimD <- max_df(
    dat = filter(dat, str_detect(group, paste0("^", name, "$"))),
    model = Q$mod_list[[dat$month[1] + 1]],
    truth = F, pred = pred)
  
  # bestDopt <- max_df(dat = optimD, model = NULL,
  #                    truth = T, pred = pred)$dose
  optimD 
  # %>% mutate(best_dose = bestDopt)
}

simMonthT <- function(dat, Q, bestDoses, pred) {
  out <- dat %>% filter(group != "greedy",
                        !str_detect(group, paste(names(Q), collapse = "|"))) %>%
    mutate(dose = ifelse(group == "1on1off",
                         ifelse(dat$month[1] %% 2 == 1, 0, 1),
                         dose))
  
  optimD <- map2_df(Q, names(Q),
                    ~ getPreds(.x, .y, dat = dat, pred = pred))
  
  if (!is.null(bestDoses)) {
    bDoses <- bestDoses %>% select_(paste0("dose", dat$month[1])) %>% flatten_dbl()
    bestD <- dat %>% filter(group == "greedy") %>% mutate(dose = bDoses)
    dat <- bind_rows(out, bestD, optimD) %>% arrange(ID) 
  } else {
    dat <- bind_rows(out, optimD) %>% arrange(ID) 
  }
  
  dat <- Mnext(dat)
  dat <- Wnext(dat)
  dat %>%
    mutate(
      lam = lambda(M_next, W_next, Z = noise_chng),
      pdeath = pexp(lam),
      # surv_time = rexp(nrow(.), lam),
      expect_surv_time = 1 / lam,
      dead = ifelse(dead, T, expect_surv_time < 1)
    )
}

filterBest <- function(bestD) {
  bestD %>% group_by(ID) %>%
    filter(expect_surv_time == max(expect_surv_time)) %>%
    select(ID, starts_with("dose")) %>% ungroup() %>%
    mutate(totDose = rowSums(.[starts_with("dose", vars = names(.))])) %>%
    group_by(ID) %>% filter(totDose == min(totDose)) %>% ungroup()
}

getBestDoses <- function(dat, Ttot) {
  if (Ttot > 3) {
    stop("Too many stages")
  }
  for (i in 0:(Ttot - 1)) {
    if (i == 0) {
      bestD <- dat %>%
        select(ID, tumor_mass, toxicity,  noise_chng, starts_with("c")) %>% 
        max_df(model = NULL, truth = T, pred = T, nested = T) %>% 
        select(ID, tumor_mass, toxicity,  noise_chng, starts_with("c"),
               M_next, W_next, expect_surv_time, starts_with("dose"))
    } else {
      bestD <- max_df(bestDalive, model = NULL, truth = T, pred = T, nested = T)
    }
    if (i < (Ttot - 1)) {
      bestD <- bestD %>% mutate(tumor_mass = M_next,
                                toxicity = W_next)
    }
    bestD <- bestD %>% 
      mutate_(.dots = setNames("dose", paste0("dose", i))) %>% select(-dose)
    bestDalive <- bestD %>% filter(expect_surv_time > 1)
    doseDead1 <- bestD %>% filter(!(ID %in% bestDalive$ID)) %>% filterBest()
    if (i == 0) {
      doseDead <- doseDead1
    } else {
      doseDead <- bind_rows(doseDead, doseDead1)
    }
    if (nrow(bestDalive) == 0) {
      break
    }
  }
  
  bind_rows(filterBest(bestDalive), doseDead) %>% arrange(ID)
}

sim_test <- function(Q, int, noise, noise_pred, npergroup = 100, Ttot = 3, seed, best = F) {
  set.seed(seed)
  ngroups <- 11 + length(Q) + ifelse(best, 1, 0)
  M0 <- runif(npergroup, min = 0, max = 2)
  W0 <- runif(npergroup, min = 0, max = 2)
  
  dat <- tibble(
    ID = 1:npergroup,
    month = 0,
    tumor_mass = M0,
    toxicity = W0,
    dead = F
  )
  
  dat <- genIntNoise(dat, int, noise, noise_pred)
  
  D1 <- rep(seq(from = 0.1, to = 1, by = 0.1), each = npergroup)
  
  groups <- c(seq(from = 0.1, to = 1, by = 0.1) %>% as.character(), names(Q), "1on1off")
  
  if (best) {
    groups <- c(groups, "greedy")
  }
  
  dat <- dat[rep(seq_len(nrow(dat)), ngroups), ] %>% 
    mutate(
      ID = rep(1:(npergroup * ngroups)),
      group = rep(groups, each = npergroup),
      dose = c(D1, rep(NA, npergroup * (ngroups - 10))),
      best_dose = NA
    )

  if (best) {
    doseDat <- filter(dat, group == "greedy") %>%
      mutate(rep = 1:n() %/% (20.000001))
    reps <- doseDat %>% group_by(rep) %>% count() %>% select(rep) %>% flatten_dbl()
    bestDoses <- map_df(reps, ~ getBestDoses(dat = filter(doseDat, rep == .x),
                                            Ttot = Ttot))
  } else {
    bestDoses <- NULL
  }
  
  d <- simMonthT(dat, Q, bestDoses = bestDoses, pred = T) %>%
    mutate(reward = ifelse(dead, log(expect_surv_time), 0))
  out <- d
  for (i in 1:(Ttot - 1)) {
    d <- d %>% mutate(month = i,
                      tumor_mass = M_next,
                      toxicity = W_next) %>%
      simMonthT(Q, bestDoses = bestDoses, pred = F) %>%
      mutate(reward = ifelse(!dead, 0, log(i + expect_surv_time)))
    out <- bind_rows(out, d)
  }
  out %>% group_by(ID) %>% mutate(
    reward = ifelse(month == Ttot - 1 & !dead, log(expect_surv_time + Ttot - 1), reward),
    pdeath = ifelse(is.na(pdeath), 1, pdeath),
    tot_reward = sum(reward[1:Ttot], na.rm = T),
    Qhat = reward,
    best = NA
  ) %>% arrange(ID) %>% ungroup()
}
