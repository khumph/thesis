set.seed(1)
noise_pred
noise
int
meow <- sim_test(Q, int = int, noise = noise, noise_pred = noise_pred)
# dat <- meow %>% filter(group == "optim")
# dat1 <- meow %>% filter(group == "1")
# dat2 <- meow %>% filter(group == "best")

d <- tibble(
  x = c(1, 2, 3, NA, NA)
)

d %>% mutate(
  y = sum(x, na.rm = T)
)

(d <- meow %>% filter(ID == 18) %>% select(ID, surv_time))

d %>% group_by(ID) %>% mutate(
  sum(surv_time[1], na.rm = T)
)

sum(d$surv_time, na.rm = T)


# 
# dat <- meow %>% filter(group == "optim" | group == "best") %>%
#   mutate(month = factor(month)) %>% 
#   group_by(group, month) %>%
#   mutate(median_dose = median(dose, na.rm = T)) %>% ungroup()
# 
# 
# meow <- meow %>% ungroup() %>% 
#   group_by(group) %>% mutate(mean_reward = mean(tot_reward))
# 
# ggplot(data = meow) +
#   # coord_trans(y = "log") +
#   # geom_boxplot(mapping = aes(
#   #   x = group,
#   #   y = tot_reward,
#   #   color = group
#   # ), notch = T) +
#   geom_point(aes(
#     x = group,
#     y = mean_reward,
#     color = group
#   ), shape = 5) +
#   labs(y = "log months of survival") 

# setdiff(dat2 %>% select(ID, dose),
#         dat1 %>% select(ID, dose))
# 
# 
# , dat1)
# 
# dat[1, ] 
# 
# maxMonth(dat, int = int, noise_pred = noise_pred, nested = T) %>% View()

# max_best(dat_long, int = int, noise_pred = noise_pred) %>% View()
# 
# best <- meow %>% filter(group == "best")
# # best1 <- meow %>% filter(group == "best", month == 0)
# # # maxMonth(best1, int, noise_pred, nested = F) %>% View()
# maxbest <- max_best(best, int, noise_pred) 
# best %>% filter(!is.na(dose)) 
# maxbest %>% filter(!is.na(dose))


# meow %>% filter(group == "best", X1 < 0.5, X2 > 0.5) %>% View()
# ids.1 <- meow %>% filter(group == "optim") %>% select(ID) %>% flatten_dbl()
# ids1 <- meow %>% filter(group == "1on1off") %>% select(ID) %>% flatten_dbl()
# idsbest <- meow %>% filter(group == "best") %>% select(ID) %>% flatten_dbl()
# indPlot(meow, sample(idsbest, 1))

# 
# dat <- filter(meow, month == 0, ID %in% c(2001:2010, 2201:2210, 2401:2410))
# 
# 
# d <- simMonthT(dat, Q, int = int, noise = noise) %>% mutate(
#   reward = ifelse(dead, log(surv_time), 0)
# )
# out <- d
# for (i in 1:(Ttot - 1)) {
#   d <- d %>% mutate(month = i,
#                     tumor_mass = M_next,
#                     toxicity = W_next) %>%
#     simMonthT(Q, int = int, noise = noise) %>%
#     mutate(reward = ifelse(!dead, 0, log(i + 1 + surv_time)))
#   out <- bind_rows(out, d)
# }
# out <- out %>% group_by(ID) %>% mutate(
#   reward = ifelse(month == Ttot - 1 & !dead, log(surv_time + Ttot - 1), reward),
#   pdeath = ifelse(is.na(pdeath), 1, pdeath),
#   tot_reward = sum(reward, na.rm = T),
#   Qhat = reward,
#   best = NA
# ) %>% arrange(ID)

# npergroup = 200
# ngroups = 13
# Ttot = 6
# set.seed(1)
# M0 <- runif(npergroup, min = 0, max = 2)
# W0 <- runif(npergroup, min = 0, max = 2)
# 
# dat <- tibble(
#   ID = 1:npergroup,
#   month = rep(0, npergroup),
#   tumor_mass = M0,
#   toxicity = W0,
#   dead = rep(F, npergroup)
# )
# 
# dat <- genIntNoise(dat, int, noise)
# 
# D1 <- rep(seq(from = 0.1, to = 1, by = 0.1), each = npergroup)
# 
# D0 <-
#   max_df(
#     data = dat,
#     model = Q$mod_list[[1]],
#     form = Q$formula,
#     mod_type = Q$mod_type,
#     pred = T
#   )$best
# 
# Dbest <- maxMonth(dat, int = int, noise = noise)$dose
# 
# D1on1off <- rep(1, npergroup)
# 
# groups <- c(
#   seq(from = 0.1, to = 1, by = 0.1) %>% as.character(),
#   "best",
#   "optim", "1on1off")
# 
# dat <- dat[rep(seq_len(nrow(dat)), ngroups), ] %>%
#   mutate(
#     ID = rep(1:(npergroup * ngroups)),
#     group = rep(groups, each = npergroup),
#     dose = c(D1, Dbest, D0, D1on1off),
#     best_dose = ifelse(group == "optim", Dbest, NA)
#   )
# 
# d <- simMonthT(dat, Q, int = int, noise = noise) %>% mutate(
#   reward = ifelse(dead, log(surv_time), 0)
# )
# out <- d
# for (i in 1:(Ttot - 1)) {
#   d <- d %>% mutate(month = i,
#                     tumor_mass = M_next,
#                     toxicity = W_next) %>%
#     simMonthT(Q, int = int, noise = noise) %>%
#     mutate(reward = ifelse(!dead, 0, log(i + 1 + surv_time)))
#   out <- bind_rows(out, d)
# }
# out <- out %>% group_by(ID) %>% mutate(
#   reward = ifelse(month == Ttot - 1 & !dead, log(surv_time + Ttot - 1), reward),
#   pdeath = ifelse(is.na(pdeath), 1, pdeath),
#   tot_reward = sum(reward, na.rm = T),
#   Qhat = reward,
#   best = NA
# ) %>% arrange(ID) %>% ungroup()
# 
# b1 <- filter(out, group == "best")
# 
# b <- max_best(filter(out, group == "best"), int = F, noise = F)


# data <- dat %>% mutate(Qhat = expected_surv)
# 
# form1 <- makeForm(form, "dose", mod_type = "caret")
# 
# d <- sample_frac(data, size = 1, replace = T)
# Q1 <- one_step_Q(form1, d,
#                  formula = form,
#                  mod_type = "caret",
#                  method = 'gcvEarth',
#                  tuneGrid = expand.grid(degree = 2))
# out <- bind_cols(Q1$data$ID, best = Q1$best)
# 
# boot <- function(data, i, mod_type, ...) {
#   
# }
# 
# if (boot) {
#   
# }