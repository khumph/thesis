source("~/Google Drive/active/setup.R")
p_load(earth)

source("./code/sim2.R")
set.seed(20161116)
dat_long <- sim()

# subset only to last month, remove people who have died
dat <- dat_long %>% filter(month == 5) %>% select(reward, tumor_mass, toxicity, X, dose) 

mars1 <- earth(reward ~ tumor_mass * dose + toxicity * dose +  X * dose, data = na.omit(dat))

summary(mars1)
plot(mars1)
predict(mars1)

x = seq(0, 1, 0.01)
idvar = "ID"
model = mars1
data <- dat

data <- data %>% mutate(ID = 1:nrow(.))
idvar <- "ID"

data_long <- data %>%
  select(matches(idvar),
         one_of(form$covariates)) %>%
  mutate(dose = map(1:nrow(.), ~ x)) %>%
  unnest()

data_preds <- data_long %>%
  mutate(preds = predict(model, .),
         preds = ifelse(is.na(preds), 0, preds)) %>%
  group_by_(idvar) %>%
  mutate(max = max(preds),
         best = (which.max(preds) - 1) / 100)

dat_long <- dat_long %>% mutate(
  Q_hat = ifelse(month == 6, NA,
                 ifelse(month == 5, reward, 0)),
  best = ifelse(month == 6, NA, 999)
)
Q <- Qlearn(dat_long, Q_hat ~ tumor_mass + toxicity + X + dose, "dose", method = "mars")

max_df(dat_long %>% filter(month == 5), model = mars1, form = form, method = "mars")


form <- makeRCS(Q_hat ~ tumor_mass + toxicity + X + dose, "dose", method = "mars")
data <- na.omit(dat_long)
model <- earth(form$formula, data)
idvar = "ID"
x = seq(0, 1, by = 0.01)
method = "mars"
nested = F
if (is.null(idvar)) {
  data <- data %>% mutate(ID = 1:nrow(.))
  idvar <- "ID"
}
data_long <- data %>%
  select(matches(idvar),
         one_of(form$covariates)) %>%
  mutate(dose = map(1:nrow(.), ~ x)) %>%
  unnest()
if (method == "rcs" | method == "mars") {
  data_preds <- data_long %>%
    mutate(preds = predict(model, .)) %>%
    group_by_(idvar) %>%
    mutate(max = max(preds),
           best = (which.max(preds) - 1) )
}
data_preds %>%
    select(-dose, -preds) %>% unique()

formula <- reward ~ 
treatment <- "dose"
method = "mars"

  form_char <- as.character(formula)
  response <- form_char[2]
  predictors <- form_char[3]
  predictor_names <- strsplit(predictors, " \\+ ")[[1]]
  if (method == "mars") {
    formula <- paste(response,
                  "~",
                  paste0(
                    predictor_names[predictor_names != treatment],
                    " * ", treatment, collapse = " + ")
                  )
  }
  covariates <- predictor_names[!(predictor_names %in% treatment)]
  list(formula = as.formula(formula), covariates = covariates, treatment = treatment)
  