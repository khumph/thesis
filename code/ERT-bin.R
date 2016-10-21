# Run all chunks in crt-simulation.rmd

# single step Q-learning with ERT and binary treatment -------------------------
library(caret)

# to make things simple to start, dichotomize treatment:
dat_bin <- dat %>% mutate_at(vars(num_range("D.", 0:5)),
                             funs(ifelse(. >= 0.5, 1, 0)))

# Tuning parameters:
#   mtry (# Randomly Selected Predictors) - also K ?
#   numRandomCuts (# Random Cuts) - K in CRT paper ?

grid <- expand.grid(mtry = 3,
                    numRandomCuts = 1:3)

set.seed(1)
mod_ert <- train(
  r.5 ~ M.5 + W.5 + D.5,
  data = dat_bin,
  method = 'extraTrees',
  tuneGrid = grid,
  ntree = 50, # G in CRT paper
  nodesize = 2 # nmin in CRT paper
)
mod <- mod_ert
plot(mod_ert)
mod$finalModel

all_trt <- dat_bin %>% select(M.5, W.5) %>% cbind(D.5 = rep(1, 1000))
all_cntrl <- dat_bin %>% select(M.5, W.5) %>% cbind(D.5 = rep(0, 1000))

preds_trt <- predict(mod, all_trt)
preds_cntrl <- predict(mod, all_cntrl)

best_ert <- data_frame(preds_trt, preds_cntrl) %>%
  mutate(best = as.numeric(preds_trt > preds_cntrl))
best_ert


# compare to output from DTRlearn funct ----------------------------------------

library(DTRlearn)
select <- dplyr::select
# examples from crt-simulation.rmd (need to run that first)
# to walk through function and see what it's doing
H <- cbind(M[, 6], W[, 6])
R <- r[, 6]
A <- ifelse(D[, 6] >= 0.5, 1, -1) # dichotimize treatment


set.seed(1)
dtr <- Qlearning_Single(H, A, R, pentype = "lasso", m = 4)
preds_dtr <- predict.qlearn(dtr, H)
head(preds_dtr, 10)


# compare to lasso implemented in caret ----------------------------------------

grid <- expand.grid(alpha = 1, # for lasso
                    lambda = seq(.0000001, 1, length.out = 100))

# lambda.min.ratio = 1.0e-6

tc <- trainControl(method = "cv",
             number = 4
             )

set.seed(1)
mod_lasso <- train(
  r.5 ~ M.5 + D.5 + W.5 + M.5 * D.5 + W.5 * D.5,
  data = dat_bin,
  method = 'glmnet',
  tuneGrid = grid,
  trControl = tc
)
mod <- mod_lasso
plot(mod)
summary(mod$finalModel)

all_trt <- dat_bin %>% select(M.5, W.5) %>% cbind(D.5 = rep(1, 1000))
all_cntrl <- dat_bin %>% select(M.5, W.5) %>% cbind(D.5 = rep(0, 1000))

preds_trt <- predict(mod, all_trt)
preds_cntrl <- predict(mod, all_cntrl)

best_lasso <- data_frame(preds_trt, preds_cntrl) %>%
  mutate(
    best = as.numeric(preds_trt > preds_cntrl)
    ) %>% cbind(dat_bin) %>% tbl_df()

best_lasso


