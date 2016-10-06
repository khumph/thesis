library(DTRlearn)

# examples from crt-simulation.rmd (need to run that first)
# to walk through function and see what it's doing
H <- cbind(M[, 6], W[, 6])
R <- r[, 6]
A <- ifelse(D[, 6] >= 0.5, 1, -1) # dichotimize treatment so it plays nice with the function

qlearn <- Qlearning_Single(H, A, R, pentype = "lasso", m = 4)

# Complete funciton (comments mine)
Qlearn_Single <- function(H, A, R, pentype = "lasso", m = 4) 
{
  # H is matrix of state variables
  # A is a vector of actions
  # R is a vector of rewards
  # pentype is estimation (either "lasso" or "LSE" - for ordinary least squares)
  # m is number of folds of cross validation in lasso
  
  n = length(A)
  
  # Create a matrix with the state variables, action variable
  # (in this case treatment: 1 or -1), and the state variables again
  # but negative for those who were in -1 group (control?)
  X = cbind(H, A, diag(A) %*% H) # matrix with H, A, -H for those in -1 group
  
  # estimate coefficients for the above
  if (pentype == "lasso") {
    cvfit = cv.glmnet(X, R, nfolds = m)
    co = as.matrix(predict(cvfit, s = "lambda.min", type = "coeff"))
  }
  else if (pentype == "LSE") {
    co = coef(lm(R ~ X))
  }
  else stop(gettextf("'pentype' is the penalization type for the regression step of Olearning, the default is 'lasso',\nit can also be 'LSE' without penalization"))
  
  # Create two matrices
  # one with a vector of ones (treatment) then state vars twice in a row 
  XX1 = cbind(rep(1, n), H, rep(1, n), H)
  # and another which is the same, but the second repeat is negative (not on treatment)
  XX2 = cbind(rep(1, n), H, rep(-1, n), -H)
  # estimate rewards from each matrix just created
  Q1 = XX1 %*% co
  Q2 = XX2 %*% co
  # figure out which matrix gives largest rewards
  Q = apply(cbind(Q1, Q2), 1, max)
  # return the coefficients and the vector of largest rewards
  Qsingle = list(co = co, Q = Q)
  class(Qsingle) = "qlearn"
  Qsingle
}

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
