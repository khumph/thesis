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
                control = quint.control(Bootstrap = F))

quint1 %>% summary()

quint1 %>% plot()
