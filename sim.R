library(tidyverse)

# simulate 100 patients
N <- 1000

# over six months of treatment
Ttot <- 6

# initialize toxicities, doses
set.seed(1)
W0 <- runif(N, min = 0, max = 2)
set.seed(1)
D0 <- runif(N, min = 0.5, max = 1)

set.seed(1)
D <- replicate(5, runif(N, min = 0, max = 1))
D <- cbind(D0, D)

data_frame(
  id = 1:1000,
  toxicity = W0,
  dose = D0,
  month = 0
) %>% 
  add_column(
    D[, 1:5]
  )

# define functions for how tumor mass & toxicity change as a result of dose, each other
Wdot <- function(M, D, a1 = 0.1, b1 = 1.2, d1 = 0.5) { 
  a1 * M + b1 * (D - d1)
}

Mdot <- function(M, W, D, a2 = 0.15, b2 = 1.2, d2 = 0.5) { 
  (a2 * W - b2 * (D - d2)) * ifelse(M > 0, 1, 0)
}

W <- matrix(c(W0, numeric(6000)), ncol = 7)
M <- matrix(c(M0, numeric(6000)), ncol = 7)

for (i in 1:(Ttot)) {
  M_next <- ifelse(M[, i] == 0, 0, Wdot(M[, i], D[, i]) + W[, i])
  W[, i + 1] <- Wdot(M[, i], D[, i]) + W[, i]
  M[, i + 1] <- ifelse(M_next <= 0, 0, M_next)
}

R2 <- function(W1, W0) {
  ifelse(W1 - W0 <= -0.5, 5,
         ifelse(W1 - W0 >= 0.5, -5, 0))
}

R3 <- function(M1, M0) {
  ifelse(M1 == 0 & M0 == 0, 0,
         ifelse(M1 == 0, 15,
                ifelse(M1 - M0 <= -0.5, 5,
                       ifelse(M1 - M0 >= 0.5, -5, 0)))
  )
  
}

r2 <- matrix(c(numeric(6000)), ncol = 6)
r3 <- matrix(c(numeric(6000)), ncol = 6)
for (i in 1:(Ttot)) {
  r2[, i] <- R2(W[, i + 1], W[, i])
  r3[, i] <- R3(M[, i + 1], M[, i])
}


lambda <- function(W, M, mu0, mu1, mu2) {
  exp(mu0 + mu1 * W + mu2 * M)
}

mu0 <- -7 ### MADE UP. NO VALUE GIVEN IN TEXT
mu1 <- 1
mu2 <- 1


lams <- matrix(numeric(6000), ncol = 6)
for (i in 1:(Ttot)) {
  lams[, i] <- lambda(W[, i + 1], M[, i + 1], mu0, mu1, mu2)
}
deltaF <- exp(-lams)


p <- 1 - deltaF

set.seed(1)
died <- apply(p, c(1, 2), function(p) rbinom(1, 1, p)) 
r1 <- ifelse(died == 1, -60, 0)

r <- r1 + r2 + r3

zeroRewards <- function(x, died = died, i) { 
  # zeroR function zeros rewards after death of a patient
  # inputs:
  #   matrix of total rewards (r) 
  #   matrix of deaths (died)
  #   an index of a particular patient of interest (i)
  # returns:
  #   row vector of rewards for the patient with rewards following death = 0
  if (sum(died[i, ]) >= 1) {
    first <- which(died[i, ] == 1)[1]
    if (first < ncol(died)) {x[i, (first + 1):ncol(died)] <- 0} 
    return(x[i, ])
  } 
  return(x[i, ])
}

for (i in 1:nrow(r)) {
  r[i, ] <- zeroRewards(r, died, i)
}


colnames(D) <- c(0:5)
colnames(M) <- c(0:6)
colnames(W) <- c(0:6)
colnames(r) <- c(0:5)
dat <- data.frame(D = D, M = M, W = W, r = r) %>% tbl_df() %>%
  rownames_to_column(var = "ID") %>% mutate(ID = factor(ID))

dat_long <- dat %>%
  gather(key, value, -ID) %>%
  extract(col = key,
          into = c("var", "month"),
          regex = "(.)\\.(.)") %>%
  spread(var, value) %>%
  select(
    ID,
    month,
    dose = D,
    tumor_mass = M,
    toxicity = W,
    reward = r
  )