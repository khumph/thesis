rm(list = ls())
ptm <- proc.time()
library(pacman)
p_load(tidyverse, rms, caret, gridExtra, stringr)
source("sim.R")
source("sim-test2.R")
source("q-functs.R")
source("results-functs.R")

samps <- 1
nstages <- 3
i <- 1

x <- list(
  orig = rep(F, 3) #,
  # noise = c(F, T, F),
  # noise_pred = c(F, T, T),
  # int = c(T, F, F),
  # int_noise = c(T, T, F),
  # int_noise_pred = rep(T, 3)
)

# boot_samp <- function(dat_long, seed) {
#   set.seed(seed)
#   ids <- sample(1:1000, 1000, replace = T)
#   map_df(0:5, ~ filter(dat_long, month == .x)[ids, ] %>%
#            mutate(ID_real = ID, ID = 1:1000))
# }

for (i in seq_along(x)) {
  int = x[[i]][1]
  noise = x[[i]][2]
  noise_pred = x[[i]][3]

# set.seed(20161116)
# dat_long <- sim(int = int, noise = noise, noise_pred = noise_pred, Ttot = 6)

dat_long <- map_df(1:samps, ~ sim(N = 100,
                                  int = int, noise = noise,
                                  noise_pred = noise_pred, Ttot = nstages,
                                  seed = 20170411 + .x) %>% mutate(samp = .x))

# bootSamps <- map_df(1:samps, ~ boot_samp(dat_long, seed = .x) %>%
#                       mutate(samp = .x))
# bootSamps <- bind_rows(dat_long %>% mutate(samp = 0), bootSamps)

if (int & !noise) {
  nms <- dat_long %>%
    select(tumor_mass, toxicity, dose, Qhat, starts_with("X")) %>% names()
} else if (noise & !int) {
  nms <- dat_long %>% select(tumor_mass, toxicity, dose, Qhat,
                             starts_with("Z"),
                             starts_with("V")) %>% names()
} else if (noise & int) {
  nms <- dat_long %>% select(tumor_mass, toxicity, dose, Qhat,
                             starts_with("X"),
                             starts_with("Z"),
                             starts_with("V")) %>% names()
} else {
  nms <- dat_long %>% select(tumor_mass, toxicity, dose, Qhat) %>% names()
}
response <- "Qhat"
treatment <- "dose"
pred_nms <- nms[!(nms %in% c(response, treatment))]
form <- paste(response,
              "~",
              paste0(c(pred_nms, treatment), collapse = " + ")) %>%
  as.formula()

ptm_rpart <- proc.time()

set.seed(20170128)
Q_rpart <- map(
  1:samps,
  ~ Qlearn(
    form = form,
    dat_long = filter(dat_long, samp == .x),
    nstages = nstages,
    # tuneGrid = expand.grid(cp = 10^(seq(-7, -1, 1))),
    # trControl = trainControl(method = "boot"),
    method = 'rpart1SE',
    control = rpart.control(cp = 1e-5,
                            maxsurrogate = 0,
                            minbucket = 5,
                            minsplit = 15)
  )
) %>% set_names(paste0("CART", 1:samps))

t_rpart <- proc.time() - ptm_rpart

# MARS --------------------------------------------------------------------

#--- predictors in PREDICTORS are allowed to interact with predictors in PARENTS
#--- but no other interactions are allowed

PREDICTORS <- pred_nms
PARENTS <-  c("dose")
allowFunct <- function(degree, pred, parents, namesx) {
  if (degree <= 1)
    return(TRUE)
  predictor <- namesx[pred]
  parents <- namesx[parents != 0]
  if ((any(predictor %in% PREDICTORS) &&
       any(parents %in% PARENTS)) ||
      (any(predictor %in% PARENTS) &&
       any(parents %in% PREDICTORS))) {
    return(TRUE)
  }
  FALSE
}

ptm_mars <- proc.time()

set.seed(20170128)
Q_mars <- map(
  1:samps, ~ Qlearn(
    form,
    filter(dat_long, samp == .x),
    nstages = nstages,
    method = 'gcvEarth',
    trControl = trainControl(method = "none"),
    nk = 200,
    allowed = allowFunct,
    tuneGrid = expand.grid(degree = 2)
  )) %>% set_names(paste0("mars", 1:samps))

(t_mars <- proc.time() - ptm_mars)


# # BMARS -------------------------------------------------------------------
# 
# ptm_bmars <- proc.time()
# 
# set.seed(20170128)
# Q_bmars <- map(
#   1:samps, ~ Qlearn(
#     form,
#     filter(dat_long, samp == .x),
#     method = 'bagEarthGCV',
#     trControl = trainControl(method = "none"),
#     nk = 200,
#     allowed = allowFunct,
#     tuneGrid = expand.grid(degree = 2)
#   )) %>% set_names(paste0("bmars", 1:samps))
# 
# t_bmars <- proc.time() - ptm_bmars

# RF ----------------------------------------------------------------------

ptm_rf <- proc.time()

# grid <- expand.grid(mtry = floor(seq(2 * length(pred_nms) / 3, length(pred_nms), length = 3)))

set.seed(20170128)
Q_rf <- map(1:samps, ~ Qlearn(
  form,
  filter(dat_long, samp == .x),
  method = 'ranger',
  nstages = nstages,
  # tuneGrid = grid,
  tuneGrid = expand.grid(mtry = length(pred_nms)),
  trControl = trainControl(method = "none"),
  min.node.size = 5,
  importance = "impurity",
  always.split.variables = c("dose"),
  num.trees = 250
)) %>% set_names(paste0("rf", 1:samps))

t_rf <- proc.time() - ptm_rf

# ptm_rf <- proc.time()
# 
# set.seed(20170128)
# Q_rf <- map(0:samps, ~ Qlearn(
#   form,
#   filter(dat_long, samp == .x),
#   method = 'ranger',
#   tuneGrid = grid,
#   trControl = trainControl(method = "oob"),
#   min.node.size = 5,
#   importance = "impurity",
#   always.split.variables = c("dose"),
#   num.trees = 500
# )) %>% set_names(paste0("rf", 0:samps))
# 
# t_rf <- proc.time() - ptm_rf

# test --------------------------------------------------------------------
ptm_test <- proc.time()

Q_list <- c(Q_rpart, Q_mars,  Q_rf)
Q_list <- Q_mars
dat_test <- sim_test(Q_list, int = int, noise = noise, noise_pred = noise_pred,
                     npergroup = 1000, Ttot = nstages, seed = 20170410)

t_test <- proc.time() - ptm_test

# save --------------------------------------------------------------------

type_nm <- paste0(
  ifelse(int, "int", ""),
  ifelse(int & noise, "_", ""),
  ifelse(noise, "noise", ""),
  ifelse(noise_pred, "_pred", ""),
  ifelse(!noise_pred & !noise & !int, "orig", "")
)
prefix <- c("dat_test_", "Q_list_")
var_nms <- map_chr(prefix, ~ paste0(.x, type_nm))
for (i in 1:length(var_nms)) {
  assign(var_nms[[i]], list(dat_test, Q_list)[[i]])
}
# save(list = var_nms, file = paste0(type_nm, ".Rdata"))
}

(t <- proc.time() - ptm)
