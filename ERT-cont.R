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

x <- seq(0, 1, by = 0.01)

get_preds <- function(x) {
  dat <- dat %>% select(M.5, W.5) %>% cbind(D.5 = rep(x, 1000))
  predict(mod, dat)
}

col_summ <- function(df, f) {
  df %>%
    keep(is_numeric) %>%
    map_dbl(f) 
}

x_char <- as.character(x)
preds <- map(x,  get_preds) %>%
  set_names(x_char) %>% tbl_df()

t5 <- dat %>% select(M.5, W.5) %>% bind_cols(preds) %>%
  gather(dose, pred, one_of(x_char)) %>% 
  mutate(dose = as.numeric(dose)) %>% 
  group_by(M.5, W.5) %>% nest() %>%
  bind_cols(
    map(t5$data, ~ col_summ(select(., pred), max)) %>%
      as_vector() %>% tbl_df() %>% set_names(nm = "max")
  ) %>%
  bind_cols(
    map(t5$data, ~ col_summ(select(., pred), ~ (which.max(.) - 1)/100)) %>%
      as_vector() %>% tbl_df() %>% set_names(nm = "optim_dose")
  )

names(t5$data) <- dat$ID

ggplot(t5$data$`1`, aes(x = dose, y = pred)) + geom_line()

# I chose the max, we could get more sophisticated in how to pick which is "best"
