# simulation functions ----------------------------------------------------

# function for how toxicity changes
updateW <- function(M, W, D, a1 = 0.1, b1 = 1.2, d1 = 0.5) {
  W_next <- a1 * M + b1 * (D - d1) + W
  ifelse(W_next > 0, W_next, 0)
}

# function for updating tumor mass
updateM <- function(M, W, D, a2 = 0.15, b2 = 1.2, d2 = 0.5) {
  # takes current value of tumor mass, outputs next value
  M_next <- ifelse(M > 0,
                   (a2 * W - b2 * (D - d2)) + M,
                   0)
  ifelse(M_next > 0, M_next, 0)
  # (a2 * W - b2 * (D - d2)) + M
}

# # reward functions
# R2 <- function(W_next, W) {
#   dW <- W_next - W
#   ifelse(
#     dW <= -0.5,
#     5,
#     ifelse(
#       dW > 0.5,
#       -5,
#       0
#     )
#   )
# }
# 
# R3 <- function(M_next, M) {
#   ifelse(M_next == 0 & M == 0,
#          0,
#          ifelse(M_next == 0,
#                 15,
#                 ifelse(
#                   M_next - M < -0.5,
#                   5,
#                   ifelse(M_next - M > 0.5,
#                          -5,
#                          0)
#                 )))
# }

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
      M_next = ifelse(!dead, updateM(tumor_mass, toxicity, dose), NA),
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
           toxicity, dose, died, dead, reward)
}

# sim testing function ----------------------------------------------------

sim_test2 <- function(Q) {
  # 200 patients per each of 11 treatments
  N <- 200 
  ngroups <- 12
  # 6 months treatment
  Ttot <- 6 
  
  # The initial values of W0 and M0 for the patients were randomly chosen from the 
  # same uniform distribution used in the training data.
  W0 <- runif(N, min = 0, max = 2)
  M0 <- runif(N, min = 0, max = 2)
  
  # treatments consisting of the estimated optimal treatment regime and each of the 10 possible fixed dose levels ranging from 0.1 to 1.0 with increments of size 0.1.
  D1 <- map(seq(from = 0.1, to = 1, by = 0.1), ~ rep(., N)) %>% flatten_dbl()
  
  # estimate optimal treatment regime for 200
  dat <-
    tibble(
      ID = 1:N,
      toxicity = W0,
      tumor_mass = M0
    )
  
  D0 <-
    max_df(
      data = dat,
      model = Q$mod_list[[1]],
      form = Q$formula,
      mod_type = Q$mod_type
    )$best
  
  Dbest <- maxMonth(dat)$dose
  
  D1 <- c(D1, Dbest, D0)
  
  D <- matrix(rep(D1, Ttot), ncol = Ttot)
  
  num_cells <- N * ngroups * Ttot 
  r <- matrix(c(numeric(num_cells)), ncol = Ttot)
  pdied <- matrix(c(numeric(num_cells)), ncol = Ttot)
  dead <- logical(N * ngroups)
  W <- matrix(c(rep(W0, ngroups), numeric(num_cells)), ncol = Ttot + 1)
  M <- matrix(c(rep(M0, ngroups), numeric(num_cells)), ncol = Ttot + 1)
  
  for (j in 1:(Ttot)) {
    for (i in 1:(N * ngroups)) {
      if (dead[i]) {
        pdied[i, j] <- 0
        M[i, j + 1] <- NA
        W[i, j + 1] <- NA
        r[i, j] <- 0
      } else {
        M[i, j + 1] <- updateM(M[i, j], W[i, j], D[i, j])
        W[i, j + 1] <- updateW(M[i, j], W[i, j], D[i, j])
        pdied[i, j] <- pDeath(M[i, j + 1], W[i, j + 1])
        r[i, j] <- reward_est(
          M[i, j + 1],
          M[i, j],
          W[i, j + 1],
          W[i, j],
          pdied[i , j]
        )
      }
      if (pdied[i, j] == 1) {
        dead[i] <- T
      }
    }
    if (j < Ttot) {
      df <- tibble(
        ID = tail(1:(N * ngroups), 200),
        toxicity = tail(W[, j + 1], 200),
        tumor_mass = tail(M[, j + 1], 200)
      )
      D[df$ID, j + 1] <- max_df(
        data = df,
        model = Q$mod_list[[j + 1]],
        form = Q$formula,
        mod_type = Q$mod_type
      )$best
      df <- tibble(
        ID = tail(1:(N * ngroups), 400) %>% head(200),
        toxicity = tail(W[, j + 1], 400) %>% head(200),
        tumor_mass = tail(M[, j + 1], 400) %>% head(200)
      )
      D[df$ID, j + 1] <- maxMonth(df)$dose
    }
  }
  
  colnames(D) <- 0:5
  colnames(M) <- 0:6
  colnames(W) <- 0:6
  colnames(r) <- 0:5
  colnames(pdied) <- 1:6
  dat_test <-
    data.frame(
      D = D,
      M = M,
      W = W,
      r = r,
      d = died
    ) %>% tbl_df()
  
  dat_test <- dat_test %>%
    rownames_to_column(var = "ID") %>%
    mutate(ID = as.numeric(ID),
           group = ifelse(ID < 2201, 
                          ifelse(ID > 2000, "best", D1), "optim"),
           group = factor(group))
  
  dat_test <- dat_test %>%
    gather(key, value, -ID, -group) %>%
    extract(col = key,
            into = c("var", "month"),
            regex = "(.)\\.(.)") %>%
    mutate(month = as.numeric(month)) %>%
    spread(var, value) %>%
    select(
      ID,
      month,
      group,
      dose = D,
      tumor_mass = M,
      toxicity = W,
      reward = r,
      died = d
    ) %>% group_by(ID) %>% mutate(
      tot_reward = sum(reward, na.rm = T)
    )
  dat_test
}


# best possible -----------------------------------------------------------

maxMonth <- function(dat, len = 101) {
  dat %>% ungroup() %>% 
    mutate(dose = map(1:nrow(.), ~ seq(0, 1, length.out = len))) %>%
    unnest() %>%
    mutate(
      M_next = updateM(tumor_mass, toxicity, dose),
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

best_possible <- function(dat_long) {
  dat <- dat_long %>% filter(month == 0) %>% mutate(month = as.numeric(month))
  d <- maxMonth(dat)
  out <- d
  for (i in 1:6) {
    d <- d %>% mutate(month = i) %>%
      select(ID,
             month,
             tumor_mass = M_next,
             toxicity = W_next,
             dose,
             pdied,
             reward) %>% maxMonth()
    out <- bind_rows(out, d)
  }
  out
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
      M_next = updateM(tumor_mass, toxicity, dose),
      W_next = updateW(tumor_mass, toxicity, dose),
      d_next = ifelse(group == "optim", optimD,
                      ifelse(group == "best", bestD, dose)), 
      pdeath = pDeath(M_next, W_next),
      reward = reward_est(M_next, tumor_mass, W_next, toxicity, pdeath)
    )
}

sim_test <- function(Q, npergroup = 200, ngroups = 12, Ttot = 6) {
  M0 <- runif(npergroup, min = 0, max = 2)
  W0 <- runif(npergroup, min = 0, max = 2)
  
  dat <- tibble(
    ID = 1:npergroup,
    month = rep(0, npergroup),
    tumor_mass = M0,
    toxicity = W0,
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