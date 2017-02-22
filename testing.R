rm(list = ls())
int <- F
noise <- T
library(pacman)
p_load(tidyverse, rms, caret, gridExtra, stringr)
source("surv-sim-functs.R")
source("surv-sim.R")
source("q-functs.R")
source("surv-sim-test.R")
source("results-functs.R")

set.seed(20161116)
Ttot = 2
dat_long <- sim(int = int, noise = noise)

vars <- c("ID", "month", "dead", "M_next", "W_next", "d_next", "reward", "best", "surv_time")
nms <- dat_long %>% select(-one_of(vars)) %>% names()
response <- "Qhat"
treatment <- "dose"
pred_nms <- nms[!(nms %in% response)]
form_base <- paste(response,
                   "~",
                   paste0(pred_nms, collapse = " + "))
ints <- paste0(pred_nms[!(pred_nms %in% treatment)], ":", treatment)
ints <- paste(ints, collapse = " + ")
formula <- paste(c(form_base, ints), collapse = " + ")
form <- as.formula(form_base)
form_int <- as.formula(formula)

set.seed(20170207)
ex_ID <- sample(dat_long$ID, 1)
indPlot(dat_long, ex_ID)

set.seed(20170128)
Q <- Qlearn(
  data = dat_long,
  formula = form,
  treatment = "dose",
  mod_type = "rcs"
)

indPlot(optDat(Q$data, F, F), ex_ID)

meow <- sim_test(Q, F, F)
sim_test(Q, F, F) %>% plots_tab()
