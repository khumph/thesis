library(caret)
library(randomForest)

set.seed(20161107)
mod_rf <- randomForest(Q_hat ~ tumor_mass + toxicity + dose,
             data = filter(dat_long, month == 5),
             na.action = na.omit,
             importance = T)

mod_rf$importance

# dose is least important

# trying changing reward to a factor
dat_rf <- dat_long %>% filter(month == 5) %>% 
  mutate(
    Q_hat = ifelse(Q_hat == -65, -60, Q_hat) %>% factor()
  )

dat_long %>% filter(month == 5) %>%
  select(Q_hat) %>% flatten_dbl() %>% factor() %>% summary()

mod_rf_fact <- randomForest(Q_hat ~ tumor_mass + toxicity + dose,
                       data = dat_rf,
                       na.action = na.omit,
                       importance = T)

mod_rf_fact$importance

mod_rf_fact

# predicts 0 pretty well, others not so much

mod <- mod_rf_fact

levels(predict(mod))[predict(mod)] %>% as.numeric()

dat_long %>% filter(month == 5) %>%
  select(ID, tumor_mass, toxicity, -dose) %>%
  mutate(dose = map(1:nrow(.), ~ seq(0, 1, 0.01))) %>%
  unnest() %>%
  mutate(p = predict(mod, .),
         p2 = levels(p)[p] %>% as.numeric(),
         preds = ifelse(is.na(p2), 0, p2)) %>%
  select(-p2, -p) %>% 
  group_by(ID) %>%
  mutate(max = max(preds),
         # picks the lowest dose that corresponds to the best reward
         best = (which.max(preds) - 1) / 100) %>% 
  select(-preds, -dose) %>% unique()




# Using caret -------------------------------------------------------------

dat_rf <- dat_long %>% filter(month == 5) %>% 
  mutate(
    Q_hat = ifelse(Q_hat == -65, -60, Q_hat) %>%
      factor(levels = c(-60, -5, 0, 5, 10), labels = c("n60", "n5", "p0", "p5", "p10"))
  )


tc <- trainControl(
  savePredictions = T,
  classProbs = T
)

mod_rf_c <- train(
  Q_hat ~ tumor_mass + toxicity + dose,
  method = "rf",
  data = dat_rf,
  na.action = na.omit,
  importance = T,
  trControl = tc
)

predict(mod_rf_c) %>% parse_number()

# ,
#          best = ifelse(is.na(max),
#                        NA,
#                        which.max(preds) - 1) / 100) %>% View()
#   select(-dose, -preds) %>% unique()



# ERT ---------------------------------------------------------------------

  mod_ert <- train(
    Q_hat ~ tumor_mass + toxicity + dose,
    method = "extraTrees",
    data = dat_rf,
    na.action = na.omit,
    # trControl = tc
  )

