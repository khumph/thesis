library(rpart)
p_load(rpart.plot, caret)
dat_long <- dat_long %>% mutate(
  Q_hat = ifelse(month == 6, NA,
                 ifelse(month == 5, reward, 0)),
  Q_hat = ifelse(Q_hat == -65, -60, Q_hat) %>% factor(),
  best = ifelse(month == 6, NA, 999)
)

dat_long$Q_hat %>% plot()

parameter <-
  rpart.control(
    # minsplit = 10,
    minbucket = 5,
    xval = 20,
    cp = 1e-6
  )

month5 <- dat_long %>% filter(month == 5)

set.seed(20161105)
mod_rpart <- rpart(Q_hat ~ tumor_mass + toxicity + dose,
                control = parameter,
                data = month5)

printcp(mod_rpart)
plotcp(mod_rpart)

pruned_min_sd <- mod_rpart$cptable %>%
  as_data_frame() %>%
  mutate(
    xerror_min = min(xerror),
    xerror_min_sd = ifelse(xerror == xerror_min, xstd, 0) %>% sum(na.rm = T)
  ) %>% filter(near(xerror, xerror_min + xerror_min_sd, tol = 0.01)) %>%
  filter(nsplit == min(nsplit)) %>% 
  select(CP) %>% flatten_dbl() %>%
  prune(mod_rpart, .)

pruned_min <- mod_rpart$cptable %>%
  as_data_frame() %>%
  filter(xerror == min(xerror)) %>% 
  filter(nsplit == min(nsplit)) %>% 
  select(CP) %>% flatten_dbl() %>%
  prune(mod_rpart, .)

prp(pruned_min)

# predict(pruned_min)[1, ] %>% which.max() %>% names() 

# by_row(as_data_frame(predict(pruned_min)), ~ which.max(.) %>% names() %>% as.numeric(), .to = "preds") %>% unnest() %>% select(preds)

predict(pruned_min) %>% apply(1, function(x) factor(names(which.max(x))))

mod <- pruned_min

month5 %>%
  select(ID, tumor_mass, toxicity, -dose) %>%
  mutate(dose = map(1:nrow(.), ~ seq(0, 1, 0.01))) %>%
  unnest() %>%
  mutate(preds = predict(mod, .) %>% apply(1, function(x) as.numeric(names(which.max(x))))) %>%
  group_by(ID) %>%
  mutate(max = max(preds),
         best = ifelse(is.na(max),
                       NA,
                       which.max(preds) - 1) / 100) %>%
  select(-dose, -preds) %>% unique()



mod_rf <- train(
  Q_hat ~ tumor_mass + toxicity + dose,
  data = month5,
  method = "rf",
  na.action = na.exclude
)
