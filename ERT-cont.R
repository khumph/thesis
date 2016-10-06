# Run all chunks in crt-simulation.rmd
library(caret)
library(rms)

# fit ERT model ---------------------------------------------------------------

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
# mod <- mod_ert
# plot(mod_ert)
# mod$finalModel


# fit LASSO model ---------------------------------------------------------

tc <- trainControl(method = "cv",
                   number = 4
)

grid <- expand.grid(alpha = 1, # for lasso
                    lambda = seq(.0000001, 1, length.out = 100))

set.seed(1)
mod_lasso <- train(
  r.5 ~ M.5 + D.5 + W.5 + M.5 * D.5 + W.5 * D.5,
  data = dat,
  method = 'glmnet',
  tuneGrid = grid,
  trControl = tc
)
# mod <- mod_lasso


# Single step Q-learning - OLS with RCS -----------------------------------

mod_rcs <- ols(
  r.5 ~ rcs(M.5) + rcs(D.5) + rcs(W.5) +
    rcs(M.5) %ia% rcs(D.5) + rcs(W.5) %ia% rcs(D.5),
  x = T,
  y = T,
  data = dat
)
mod <- mod_rcs

# residual plot
ggplot(dat, aes(x = predict(mod), y = residuals(mod, "dffits"))) + geom_point()

ggplot(dat, aes(x = M.5, y = r.5)) + geom_point() +
  geom_point(aes(x = M.5, y = predict(mod)), color = "blue")

ggplot(dat, aes(x = W.5, y = r.5)) + geom_point() +
  geom_point(aes(x = W.5, y = predict(mod)), color = "blue")

ggplot(dat, aes(x = predict(mod), y = r.5)) + geom_point()


## Get predictions, pick best

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

t5 <- dat %>% select(ID, M.5, W.5) %>% bind_cols(preds) %>%
  gather(dose, pred, one_of(x_char), -ID) %>% 
  mutate(dose = as.numeric(dose)) %>% 
  group_by(ID, M.5, W.5) %>% nest() %>% 
  # make a new variable from the maximium predicted rewards
  bind_cols(
    map(.$data, ~ col_summ(select(., pred), max)) %>%
      as_vector() %>% tbl_df() %>% set_names(nm = "max")
  ) %>%
  # make a new variable from the dose corresponding to the max reward
  bind_cols(
    map(.$data, ~ col_summ(select(., pred), ~ (which.max(.) - 1)/100)) %>%
      as_vector() %>% tbl_df() %>% set_names(nm = "optim_dose")
  )

names(t5$data) <- dat$ID

# plot predicted rewards across doses for a patient
ggplot(t5$data$`11`, aes(x = dose, y = pred)) + geom_point()

# plot predicted rewards across doses for random sample of patients
set.seed(1)
ggplot(data = filter(unnest(t5), ID %in% sample(1:1000, 50))) +
  geom_line(mapping = aes(x = dose, y = pred, group = ID))
