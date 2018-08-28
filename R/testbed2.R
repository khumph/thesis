rm(list = ls())
ptm <- proc.time()
library(pacman)
p_load(tidyverse, rms, caret, rpart, gridExtra, stringr)
source("sim.R")
source("sim-test-best-deter.R")
source("q-functs.R")
source("results-functs.R")

samps <- 100
nstages <- 3
best <- F
rpart <- F
MARS <- F
BMARS <- F
RF <- T

x <- list(
  # orig = rep(F, 3),
  # noise = c(F, T, F),
  # noise_pred = c(F, T, T),
  # int = c(T, F, F),
  int_noise = c(T, T, F)
  # int_noise_pred = rep(T, 3)
)

for (i in seq_along(x)) {

int = x[[i]][1]
noise = x[[i]][2]
noise_pred = x[[i]][3]

dat_long <- map_df(21:samps, ~ sim(N = 1000,
                                  int = int, noise = noise,
                                  noise_pred = noise_pred, Ttot = nstages,
                                  seed = 20170411 + .x) %>% mutate(samp = .x))

form <- makeForm(dat_long, int = int, noise = noise)
pred_nms <- form$pred_nms
form <- form$form

if (rpart) {
  ptm_rpart <- proc.time()
  
  set.seed(20170128)
  Q_rpart <- map(
    21:samps,
    ~ Qlearn(
      form = form, 
      dat_long = filter(dat_long, samp == .x),
      nstages = nstages,
      method = 'rpart',
      control = rpart.control(cp = 1e-5,
                              maxcompete = 0,
                              maxsurrogate = 0,
                              minbucket = 5))
  ) %>% set_names(paste0("CART", 21:samps))
  (t_rpart <- proc.time() - ptm_rpart)
}

if (MARS) {
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
    21:samps, ~ Qlearn(
      form,
      filter(dat_long, samp == .x),
      nstages = nstages,
      method = 'gcvEarth',
      trControl = trainControl(method = "none"),
      nk = 200,
      allowed = allowFunct,
      tuneGrid = expand.grid(degree = 2)
    )) %>% set_names(paste0("mars", 21:samps))
  
  (t_mars <- proc.time() - ptm_mars)
}

if (BMARS) {
  # BMARS -------------------------------------------------------------------
  
  ptm_bmars <- proc.time()
  
  set.seed(20170128)
  Q_bmars <- map(
    1:samps, ~ Qlearn(
      form,
      filter(dat_long, samp == .x),
      method = 'bagEarthGCV',
      B = 20,
      nstages = nstages,
      trControl = trainControl(method = "none"),
      nk = 200,
      allowed = allowFunct,
      tuneGrid = expand.grid(degree = 2)
    )) %>% set_names(paste0("bmars", 1:samps))
  
  (t_bmars <- proc.time() - ptm_bmars)
  
}

if (RF) {
  # RF ----------------------------------------------------------------------
  
  ptm_rf <- proc.time()
  
  set.seed(20170128)
  Q_rf <- map(21:samps, ~ Qlearn(
    form,
    filter(dat_long, samp == .x),
    method = 'ranger',
    nstages = nstages,
    tuneGrid = expand.grid(mtry = length(pred_nms)),
    trControl = trainControl(method = "none"),
    min.node.size = 5,
    importance = "impurity",
    always.split.variables = c("dose"),
    num.trees = 250
  )) %>% set_names(paste0("rf", 21:samps))
  
  (t_rf <- proc.time() - ptm_rf)
  
  # ptm_rf <- proc.time()
  # grid <- expand.grid(mtry = floor(seq(2 * length(pred_nms) / 3, length(pred_nms), length = 3)))
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
}

# test --------------------------------------------------------------------
ptm_test <- proc.time()
# Q_list <- c(Q_rpart, Q_mars, Q_bmars, Q_rf)
Q_list <- Q_rf
dat_test <- sim_test(Q_list, int = int, noise = noise, noise_pred = noise_pred,
                     npergroup = 2000, Ttot = nstages, seed = 20170410,
                     best = best)
(t_test <- proc.time() - ptm_test)

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
save(list = var_nms, file = paste0(type_nm, "_rf.Rdata"))

rm(list = c(var_nms, "dat_test", "dat_long", "Q_list",
            "Q_mars", "Q_rf", "Q_rpart", "Q_bmars"))
}

# Save all in one file ----------------------------------------------------

x <- list(
  orig = rep(F, 3),
  noise = c(F, T, F),
  noise_pred = c(F, T, T),
  int = c(T, F, F),
  int_noise = c(T, T, F),
  int_noise_pred = rep(T, 3)
)


for (i in seq_along(x)) {
  int = x[[i]][1]
  noise = x[[i]][2]
  noise_pred = x[[i]][3]
  type_nm <- paste0(
    ifelse(int, "int", ""),
    ifelse(int & noise, "_", ""),
    ifelse(noise, "noise", ""),
    ifelse(noise_pred, "_pred", ""),
    ifelse(!noise_pred & !noise & !int, "orig", "")
  )
  load(paste0(type_nm, "_rf.Rdata"))
  prefix <- c("dat_test_", "Q_list_")
  var_nms <- map(prefix, ~ paste0(.x, type_nm)) %>%
    set_names(c("dat_test", "Q_list"))
  dat <- var_nms$dat_test %>% parse(text = .) %>% eval()
  if (!int & !noise & !noise_pred) {
    dat_all <- dat %>% mutate(scenario = type_nm)
  } else {
    dat_all <- bind_rows(dat_all, dat %>% mutate(scenario = type_nm))
  }
  rm(list = var_nms %>% flatten_chr())
}
dat_all <- dat_all %>% select(-starts_with("V"), -starts_with("Z"), -preds)
dat_all <- dat_all %>% filter(scenario == "int", group == "greedy") %>%
  mutate(scenario = "int_noise") %>% bind_rows(dat_all)
dat_all <- dat_all %>% filter(scenario == "orig", group == "greedy") %>%
  mutate(scenario = "noise") %>% bind_rows(dat_all)
dat_all <- dat_all %>% group_by(group) %>% mutate(
  mod = str_extract(group,"^[a-zA-Z]+|^\\d\\.\\d|^1$|^1on1off"),
  samp = ifelse(str_detect(group, "^[a-zA-Z]+\\d+$"),
                str_extract(group, "\\d+$"), "0")) %>% ungroup()
save(dat_all, file = "all_rf.Rdata")

library(feather)
write_feather(dat_all, 'all_rf.feather')

(t <- proc.time() - ptm)
