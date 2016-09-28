# Run all chunks in crt-simulation.rmd

# single step Q-learning with ERT and continuous treatment ---------------------
library(caret)

# Tuning parameters:
#   mtry (# Randomly Selected Predictors) - also K ?
#   numRandomCuts (# Random Cuts) - K in CRT paper ?

grid <- expand.grid(mtry = 3,
                    numRandomCuts = 1:3)

set.seed(1)
mod_ert <- train(
  r.5 ~ M.5 + W.5 + D.5,
  data = dat,
  method = 'extraTrees',
  tuneGrid = grid,
  ntree = 50, # G in CRT paper
  nodesize = 2 # nmin in CRT paper
)
mod <- mod_ert
plot(mod_ert)
mod$finalModel

preds <- predict(mod)  
