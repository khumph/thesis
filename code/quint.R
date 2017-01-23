library(quint)

#### applied to our data ####

source("~/Google Drive/active/setup.R")
p_load(gridExtra, rpart, quint)

source("./code/sim2.R")
set.seed(20161116)
dat_long <- sim()

# subset only to last month, remove people who have died
dat <- dat_long %>% filter(month == 1 & !is.na(tumor_mass)) %>% mutate(
  dose_cat = ifelse(dose < 0.5, "low", "high") %>% as.factor())

set.seed(20170102)

quint1 <- quint(reward ~ dose_cat | tumor_mass + toxicity + X, data = dat,
                control = quint.control(B = 5))
quint1 %>% summary()
quint1 %>% plot()

q1_p <- prune(quint1)
q1_p %>% summary() 
q1_p %>% plot()

dat2 <- dat_long %>% filter(month == 2 & !is.na(tumor_mass)) %>% mutate(
  dose_cat = ifelse(dose < 0.5, "low", "high") %>% as.factor())

set.seed(20170120)

quint2 <- quint(reward ~ dose_cat | tumor_mass + toxicity + X, data = dat2,
                control = quint.control(B = 5))
quint2 %>% summary()
quint2 %>% plot()
