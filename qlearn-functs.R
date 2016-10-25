# Q-learning functions ----------------------------------------------------

makeRCS <- function(formula, treatment) {
  form_char <- as.character(formula)
  response <- form_char[2]
  predictors <- form_char[3]
  predictor_names <- strsplit(predictors, " \\+ ")[[1]]
  form_base_rcs <- paste(response,
                         "~",
                         paste0("rcs(", predictor_names, ")", collapse = " + "))
  ints <-
    paste0("rcs(", predictor_names, ")", " %ia% ", "rcs(", treatment, ")")
  trtbytrt <-
    paste0("rcs(", treatment, ")", " %ia% ", "rcs(", treatment, ")")
  ints <- ints[ints != trtbytrt]
  ints <- paste(ints, collapse = " + ")
  formula_rcs <- paste(c(form_base_rcs, ints), collapse = " + ")
  list(formula = as.formula(formula_rcs), predictors = predictor_names, treatment = treatment)
}

subset_df <- function(occasion, data) {
  filter(data, month == occasion)
}

fit_rcs <- function(formula, data, ...) {
  ols(
    formula,
    x = T,
    y = T,
    data,
    ...
  )
}

get_preds_df <- function(x, df, mod) {
  x <- seq(0, 1, by = 0.01)
  get_preds <- function(df, mod) {
    df <- df %>%
      select(tumor_mass, toxicity) %>%
      mutate(dose = rep(x, nrow(df)))
    predict(mod, df)
  }
  map(x,  ~ get_preds(., d, mod)) %>%
    set_names(as.character(x)) %>% tbl_df()
}


nest_max_df <- function(df, preds) {
  df %>%
    select(ID, tumor_mass, toxicity) %>%
    bind_cols(preds) %>%
    gather(dose, pred, one_of(names(preds)), -ID) %>% 
    mutate(dose = as.numeric(dose)) %>% 
    group_by(ID, tumor_mass, toxicity) %>%
    nest() %>% 
    mutate(
      max = map_dbl(.$data, ~ max(.$pred)),
      dose_optim = ifelse(
        is.na(max),
        NA,
        map_dbl(.$data,
                ~ min(
                  (which.max(
                    ifelse(is.na(.$pred), 999, .$pred)
                  ) - 1) / 100)
        )
      )
    )
}

max_df <- function(occasion, mod) {
  df <- subset_df(occasion)
  # mod <- fit_rcs(df)
  preds <- get_preds_df(x, df = subset_df(occasion), mod)
  nest_max_df(df, preds)
}

maxQ <- function(occasion, mod) {
  max_df(occasion, mod)$max
}

one_step_Q <- function(occasion, mod) {
  subset_df(occasion - 1)$reward + maxQ(occasion, mod)
}


# Qlearn <- function(formula, treatment, data) {
#   
# }


# testing -----------------------------------------------------------------



d <- dat_long %>% mutate(noise = runif(7000), Q_hat = reward) %>% filter(month == 5)
form <- makeRCS(Q_hat ~ noise + tumor_mass + dose + toxicity, "dose")

mod <- fit_rcs(form$formula, data = d)

get_preds <- function(data, mod) {

}

x <- seq(0, 1, by = 0.01)
d %>%
  select(ID, one_of(form$predictors), -matches(form$treatment)) %>%
  mutate(dose = map(1:nrow(.), ~ x)) %>%
  unnest() %>%
  mutate(preds = predict(mod, .)) %>%
  group_by(ID) %>%
  mutate(max = max(preds),
         dose_optim = ifelse(is.na(max),
                             NA,
                             which.max(preds) - 1) / 100) 
