# source("qlearn-functs.R")

library(broom)


# Q-Learning  -------------------------------------------------------------

set.seed(20161027)
dat_long <- dat_long %>% mutate(
  # noise = runif(7000)
  Q_hat = ifelse(month == 6, NA,
                 ifelse(month == 5, reward, 0)),
  best = ifelse(month == 6, NA, 999)
)

mod_list <- list()

Q <- Qlearn(dat_long, Q_hat ~ 
              # noise + 
              tumor_mass + dose + toxicity, "dose")

dat <- Q$data %>%
  select(ID,
         month,
         dose,
         best,
         reward,
         Q_hat,
         tumor_mass,
         toxicity,
         died) %>%
  mutate(reward = lag(reward),
         Q_hat = lag(Q_hat),
         best = ifelse(died != 1 | is.na(died), best, NA))
