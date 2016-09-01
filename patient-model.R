library(pacman)
p_load(dplyr)

# p 11

# We generate a simulated clinical reinforcement trial with N = 1000 patients (replicates)
N <- 1000

# with each simulated patient experiencing 6 months (T = 6) of treatment based on this ODE model. 
t <- 6

# The initial values W0 and M0 for each patient are generated from independent uniform (0, 2) deviates.
# W0 indicates the initial value of patient's wellness
W0 <- runif(N, min = 0, max = 2)

# M0 indicates the value of tumor size when the patient is at the beginning of the study
M0 <- runif(N, min = 0, max = 2)

# The treatment set consists of doses of a chemotherapy agent with an acceptable dose range of [0, 1], where the value 1 corresponds to the maximum acceptable dose. The values chosen for chemotherapy drug level D0 are simulated from the uniform (0.5, 1) distribution,
D0 <- runif(N, min = 0.5, max = 1)

# moreover, D1,..., D5 are drawn according to a uniform distribution in the interval (0, 1). 
D <- replicate(5, runif(N, min = 0, max = 1))

D <- cbind(D0, D)

# Thus our treatment set is restricted differently at decision time t = 0 than at other decision times to reflect a requirement that patients receive at least some drug at onset of treatment. 

# Thus, the model we present must exhibit: 
# (1) tumor growth in the absence of chemotherapy; 
# (2) patients' negative wellness outcomes in response to chemotherapy; 
# (3) the drug's capability for killing tumor cells while also increasing toxicity; and 
# (4) an interaction between tumor cells and patient wellness. 
# To obtain data which satisfy these requirements, we propose using a system of ordinary difference equations (ODE) modeled as follows:

# where time (with month as unit) t = 0, 1,..., T − 1, and Ẇt and Ṁt indicate transition functions. Note that these changing rates yield a piecewise linear model over time. Without loss of trade-off between toxicity and efficacy, the piecewise linear model can be implemented very easily. For simplicity, we here consider tumor size instead of number of tumor cells. 
# Mt denotes the tumor size at time t. 
# Wt measures the negative part of wellness (toxicity). 
# Dt denotes the chemotherapy agent dose level. 


# The value of other different parameters for the model are fixed as: a1 = 0.1, a2 = 0.15, b1 = 1.2, b2 = 1.2, d1 = 0.5 and d2 = 0.5. 

a1 <- 0.1
a2 <- 0.15
b1 <- 1.2
b2 <- 1.2
d1 <- 0.5
d2 <- 0.5

Wdot <- function(M, D, a1 = 0.1, b1 = 1.2, d1 = 0.5) { 
  a1 * M + b1 * (D - d1)
}

Mdot <- function(M, W, D, a2 = 0.15, b2 = 1.2, d2 = 0.5) { 
  (a2 * W - b2 * (D - d2)) * ifelse(M > 0, 1, 0)
}

# The indicator function term 1{Mt > 0} in (3) represents the feature that when the tumor size is absorbed at 0, the patient has been cured, and there is no future recurrence of the tumor.

W <- matrix(c(W0, numeric(5000)), ncol = 6)
M <- matrix(c(M0, numeric(5000)), ncol = 6)

for (i in 1:(t - 1)) {
  W[, i + 1] <- Wdot(M[, i], D[, i]) + W[, i]
  M[, i + 1] <- Mdot(M[, i], W[, i], D[, i]) + M[, i]
}

colnames(D) <- c(0:5)
colnames(M) <- c(0:5)
colnames(W) <- c(0:5)
dat <- data.frame(D = D, M = M, W = W) %>% tbl_df()

# We decompose this reward function Rt into three parts: Rt,1(Dt, Wt+1, Mt+1) due to survival status, Rt,2(Wt, Dt, Wt+1) due to wellness effects, and Rt,3(Mt, Dt, Mt+1) due to tumor size effects. It can be described by:
#   otherwise,

R2 <- matrix(c(numeric(6000)), ncol = 6)
R3 <- matrix(c(numeric(6000)), ncol = 6)
for (i in 1:(t - 1)) {
  R2[, i] <- ifelse(W[, i + 1] - W[, i] <= -0.5, 5,
                           ifelse(W[, i + 1] - W[, i] >= 0.5, -5, 0))
  R3[, i] <- ifelse(M[, i + 1] == 0, 15,
                    ifelse(M[, i + 1] - M[, i] <= -0.5, 5,
                           ifelse(M[, i + 1] - M[, i] >= 0.5, -5, 0)))
}



# R1 <- function(D, W, M) {
#   if return(-60)