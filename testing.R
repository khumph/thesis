set.seed(1)
meow <- sim_test(Q, int = int, noise = noise) %>% filter(group == "best")
meow  %>% View()

ids.1 <- meow %>% filter(group == "optim") %>% select(ID) %>% flatten_dbl()
ids1 <- meow %>% filter(group == "1on1off") %>% select(ID) %>% flatten_dbl()
idsbest <- meow %>% filter(group == "best") %>% select(ID) %>% flatten_dbl()
indPlot(meow, sample(idsbest, 1))