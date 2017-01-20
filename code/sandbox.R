# redefining loss function in rpart ----------------------------------------

library(pacman)
p_load(randomForest, rattle, rpart.plot, RColorBrewer)

fancyRpartPlot(mytree)

data <- dat_long %>% filter(month == 5) %>%
  mutate(Q_hat = ifelse(Q_hat == -65 | Q_hat == -55, -60,
                        Q_hat) %>% factor())
formula <- Q_hat ~ tumor_mass + dose + toxicity
method <- "rcs"
treatment <- "dose"
form <- makeRCS(formula = formula, treatment = "dose", method = "rpart")
model <- fit_rpart(form$formula, data = data)
idvar <- "ID"
method = "rpart"
x <- seq(0,1,0.01)
max_df(data = data, model = model, form, method = "rpart")

fit <- rpart(
  formula,
  data = data,
  method = "class",
  parms = list(
    split = "information",
    loss = rbind(
      c(0, 55, 60, 65, 70, 75),
      c(0, 55, 60, 65, 70, 75),
      c(0, 55, 60, 65, 70, 75),
      c(0, 55, 60, 65, 70, 75),
      c(0, 55, 60, 65, 70, 75),
      c(0, 55, 60, 65, 70, 75)
    )
  ),
  control = rpart.control(usesurrogate = 0, maxsurrogate = 0)
)

data$Q_hat

fit_rpart <- function(formula, data, cpmethod = "min", ...) {
  mod_rpart <- rpart(formula, 
                     data = data, ...)
  if (cpmethod == "min_sd") {
    pruned <- mod_rpart$cptable %>%
      as_data_frame() %>%
      mutate(
        xerror_min = min(xerror),
        xerror_min_sd = ifelse(xerror == xerror_min, xstd, 0) %>% sum(na.rm = T)
      ) %>% filter(near(xerror, xerror_min + xerror_min_sd, tol = 0.01)) %>%
      filter(nsplit == min(nsplit)) %>% 
      select(CP) %>% flatten_dbl() %>%
      prune(mod_rpart, .)
  }
  if (cpmethod == "min") {
    pruned <- mod_rpart$cptable %>%
      as_data_frame() %>%
      filter(xerror == min(xerror)) %>% 
      filter(nsplit == min(nsplit)) %>% 
      select(CP) %>% flatten_dbl() %>%
      prune(mod_rpart, .)
  }
  pruned
}

rpart(Q$formula$formula, data = filterdat_long, )



# testing -----------------------------------------------------------------


# 
# if (is.null(idvar)) {
#   data <- data %>% mutate(ID = 1:nrow(.))
#   idvar <- "ID"
# }
# data_long <- data %>%
#   select(matches(idvar),
#          one_of(form$covariates)) %>%
#   mutate(dose = map(1:nrow(.), ~ x)) %>%
#   unnest() 
# if (method == "rcs") {
#   data_preds <- data_long %>%
#     mutate(preds = predict(model, .))
# }
# if (method == "rpart") {
#   data_preds <- data_long %>%
#     mutate(preds = predict(model, .) %>%
#              apply(1, function(x) as.numeric(names(which.max(x)))))
# }
# 
# predict(model, data_long)
# 
# 
# data_preds %>%
#   group_by_(idvar) %>%
#   mutate(max = max(preds),
#          best = ifelse(is.na(max),
#                        NA,
#                        which.max(preds) - 1) / 100) %>%
#   select(-dose, -preds) %>% unique()


# 
# # Q1 <- one_step_Q(
# #   Q_hat ~ noise + tumor_mass + dose + toxicity,
# #   data = mutate(dat_long, noise = runif(7000), Q_hat = reward) %>% filter(month == 5),
# #   treatment = "dose"
# # )
# 
#   form_char <- as.character(formula)
#   response <- form_char[2]
#   predictors <- form_char[3]
#   predictor_names <- strsplit(predictors, " \\+ ")[[1]]
#   if (method == "rcs") {
#     form_base_rcs <- paste(response,
#                            "~",
#                            paste0("rcs(", predictor_names, ")", collapse = " + "))
#     ints <-
#       paste0("rcs(", predictor_names, ")", " %ia% ", "rcs(", treatment, ")")
#     trtbytrt <-
#       paste0("rcs(", treatment, ")", " %ia% ", "rcs(", treatment, ")")
#     ints <- ints[ints != trtbytrt]
#     ints <- paste(ints, collapse = " + ")
#     formula <- paste(c(form_base_rcs, ints), collapse = " + ")
#   }
#   covariates <- predictor_names[!(predictor_names %in% treatment)]
#   list(formula = as.formula(formula), covariates = covariates, treatment = treatment)
# 
# data %>%
#   select(matches(idvar),
#          one_of(form$covariates)) %>%
#   mutate(dose = map(1:nrow(.), ~ x)) %>%
#   unnest() %>%
#   mutate(preds = predict(model, .)) %>%
#   group_by_(idvar) %>%
#   mutate(max = max(preds),
#          best = ifelse(is.na(max),
#                        NA,
#                        which.max(preds) - 1) / 100) %>%
#   select(-dose, -preds) %>% unique()



# random forest -----------------------------------------------------------




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

nested_df <- dat_long %>% filter(month == 5) %>%
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
         best1 = (which.max(preds) - 1) / 100,
         best2 = ifelse(preds == max, dose, NA) %>% median(na.rm = T))
  # select(-preds, -dose) %>% unique()


set.seed(2)
i <- sample(nested_df$ID, 1)
ggplot(filter(nested_df, ID == i)) +
  geom_point(aes(x = dose, y = preds), shape = 2) +
  geom_point(aes(x = dose, y = best1), color = "red")

set.seed(5)
ggplot(data = filter(nested_df, ID %in% sample(1:1000, 50))) +
  geom_jitter(mapping = aes(x = dose, y = preds, color = factor(ID), alpha  = 0.5), show.legend = F) +
  geom_jitter(mapping = aes(x = best1, y = preds, group = factor(ID)), color = "blue", shape = 2, show.legend = F) +
  geom_jitter(mapping = aes(x = best2, y = preds, group = factor(ID)), color = "red", shape = 2, show.legend = F)

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

