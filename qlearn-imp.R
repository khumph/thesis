# Run all chunks in crt-simulation.rmd

#### DTRlearn package ####
X <- list()
for (i in 1:(ncol(M) - 1)) {
  X[[i]] <- cbind(M[, i], W[, i])
}
D_df <- data.frame(D)
r_df <- data.frame(r)

# install.packages("DTRlearn")
library(DTRlearn)
Q <- Qlearning(X = X, AA = D_df, RR = r_df, K = 6, pentype = "lasso", m = 4)
plot.qlearn(Q[[1]])
preds <- predict.qlearn(Q[[1]], X[[1]])

# The optimal treatment option is the the sign of the interaction term
# which maximize the predicted value from the regression model.

# Not sure this really helps'




#### 

# 1. Inputs: a set of training data consists of attributes x (states st, actions at) and index y (rewards rt), i.e. {(st, at, rt)i, t = 0,..., T, i = 1,..., N}.
# 2. Initialization: Let t = T + 1 and Q ̂T+1 be a function equal to zero on St × At.

# 3. Iterations: repeat computations until stopping conditions are reached (t = 0).
for (t in Ttot:1) {
  # a. t <- t − 1 .
  # b. Qt is fitted with support vector regression (SVR) or extremely randomized trees (ERT) through the following recursive equation:
  # Qt(st, at) = rt + max_{at+1}Qhat_t+1(st+1, at+1) + err 
  Qh[t, ] <- r[t] + 
    # c. Use cross-validation to choose tuning parameters C and ζ if fitting Qt via SVR with Gaussian kernel; choose plausible values of parameters K, G, nmin if fitting Qt via ERT (K = 3,G = 50, nmin = 2 in our simulation).
}
# 4. Given the sequential estimates of {Q̂0, Q̂1,..., Q̂5}, the sequential individualized optimal polices {π̂0,..., π̂5} for application to the virtual phase III trial are computed.




# Initialize Q(s, a), for all s in \mathcal{S}, a in \mathal{A}(s), arbitrarily, and Q(terminal-state, ·) = 0 
Q <- matrix(c(numeric(6000)), ncol = 6)
# Repeat (for each episode):
# Initialize S
# Repeat (for each step of episode):
for (i in 1:(Ttot)) {
# Choose A from S using policy derived from Q (e.g., epsilon-greedy) Take     action A, observe R, S'
# Q(S, A) <- Q(S, A) + alpha * (R + gamma * max_a Q(S',a) - Q(S, A))
# S <- S'
}
# until S is terminal





library(rms)
Q5 <- ols(r.5 ~ rcs(D.5) + rcs(M.5) + rcs(W.5) + rcs(D.5) %ia% rcs(M.5) + rcs(D.5) %ia% rcs(W.5), data = dat)

predict(Q5)
summary(Q5)

ggplot(dat, aes(x = D.0, y = r.0)) + geom_point() + geom_smooth()

#### other packages: ####

# MDPtoolbox
# iqlearn
# qLearn


preds <- map(x,  get_preds) %>% set_names(paste("d", x, sep = '')) %>% tbl_df()
best <- by_row(preds, max, .collate = "rows", .to = "max") %>%
  by_row(which.max, .collate = "rows", .to = "which_max") %>% 
  mutate(
    which_max = (which_max - 1)/100
  )

best %>% select(max, which_max)
