library(quint)
data("bcrp")
head(bcrp)
# help(bcrp)

ex_data <- subset(bcrp, cond < 3)

formula_dep <- I(cesdt1 - cesdt3) ~ cond | cesdt1 + negsoct1 + uncomt1 + disopt1 + comorbid + age + wcht1 + nationality + marital + trext

formula_phys <- I(physt3 - physt1) ~ cond | cesdt1 + negsoct1 + uncomt1 + disopt1 + comorbid + age + wcht1 + nationality + marital + trext

# Uisng default values for tuning parameters

set.seed(47)
quint1 <- quint(formula_dep, data = ex_data)

summary(quint1)

quint1$fi # fit information
quint1$si # split information
quint1$li # leaf information

quint1pr <- prune(quint1) # does 1 se pruning

plot(quint1pr)

round(quint1pr$li, digits = 2)

set.seed(48)
quint2 <- quint(formula_phys, data = ex_data)

round(quint2$li, digits = 2)
plot(quint2)


# use difference in means instead of effect size
control3 <- quint.control(crit = "dm") 

set.seed(48)
quint3 <- quint(formula = formula_phys, data = ex_data, control = control3)
# same as before

# stop at two leaves
control4 <- quint.control(maxl = 2)

set.seed(48)
quint4 <- quint(formula_dep, data = ex_data, control = control4)
round(quint4$li, digits = 2)


# increase d_min because of small sample size
control5 <- quint.control(crit = "dm", dmin = 0.40)
set.seed(48)
quint5 <- quint(formula = formula_phys, data = ex_data, control = control5)


