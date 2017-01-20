
dose <- numeric(200)
reward <- numeric(200)
for (i in 1:200) {
  f <- function(x) {
    predict(mod, data.frame(dat, dose = x)[i, ])
  }
  opt <- optimize(f, interval = c(0,1), tol = 0.01, maximum = T)
  dose[i] <- opt$maximum
  reward[i] <- opt$objective
}

x <- seq(0, 1, by = 0.01)
dat %>%
  mutate(dose = map(1:nrow(.), ~ x)) %>%
  unnest() %>%
  mutate(preds = predict(mod, .)) %>%
  group_by(ID) %>%
  mutate(max = max(preds),
         best = ifelse(is.na(max),
                       NA,
                       which.max(preds) - 1) / 100) %>%
  select(-dose, -preds) %>% unique()
