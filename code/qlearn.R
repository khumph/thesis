# source("sim.R")

# Q-learning functions ----------------------------------------------------

subset_df <- function(occasion, data = dat_long) {
  filter(data, month == occasion)
}

fit_rcs <- function(df) {
  ols(
    Q_hat ~
      rcs(tumor_mass) +
      rcs(dose) +
      rcs(toxicity) +
      rcs(tumor_mass) %ia% rcs(dose) +
      rcs(toxicity) %ia% rcs(dose),
    x = T,
    y = T,
    data = df
  )
}

get_preds_df <- function(x, df, mod) {
  get_preds <- function(x, df, mod) {
    df <- df %>%
      select(tumor_mass, toxicity) %>%
      mutate(dose = rep(x, nrow(df)))
    predict(mod, df)
  }
  map(x,  ~ get_preds(., df, mod)) %>%
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
  x <- seq(0, 1, by = 0.01)
  preds <- get_preds_df(x, df = subset_df(occasion), mod)
  nest_max_df(df, preds)
}

maxQ <- function(occasion, mod) {
  max_df(occasion, mod)$max
}

one_step_Q <- function(occasion, mod) {
    subset_df(occasion - 1)$reward + maxQ(occasion, mod)
}


# Q-Learning  -------------------------------------------------------------


dat_long <- dat_long %>% mutate(
  Q_hat = ifelse(month == 6, NA,
                 ifelse(month == 5, reward, 0)),
  dose_optim = ifelse(month == 6, NA, 999)
)

mod_list <- list()

for (i in 4:0) {
  subsetted_df <- subset_df(i + 1)
  mod <- fit_rcs(subsetted_df)
  mod_list[[i + 2]] <- mod
  df_max <- max_df(i + 1, mod = mod)
  dat_long[dat_long$month == (i), ]$Q_hat <- subset_df(i)$reward + df_max$max
  dat_long[dat_long$month == (i + 1), ]$dose_optim <- df_max$dose_optim
  if (i == 0) {
    subsetted_df <- subset_df(i)
    mod <- fit_rcs(subsetted_df)
    mod_list[[i + 1]] <- mod
    df_max <- max_df(i, mod = mod)
    dat_long[dat_long$month == i, ]$dose_optim <- df_max$dose_optim
  }
}

dat <- dat_long %>%
  select(ID,
         month,
         dose,
         dose_optim,
         reward,
         Q_hat,
         tumor_mass,
         toxicity,
         died) %>%
  mutate(reward = lag(reward),
         Q_hat = lag(Q_hat),
         dose_optim = ifelse(died != 1 | is.na(died), dose_optim, NA))

# dat %>% View()

ggplot(filter(dat, month == 4)) +
  geom_point(aes(x = tumor_mass, y = reward)) +
  geom_point(aes(x = tumor_mass, y = Q_hat), color = "blue")

