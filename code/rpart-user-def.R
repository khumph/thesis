# initialization function

itemp <- function(y, offset, parms, wt) {
  if (is.matrix(y) && ncol(y) > 1)
    stop("Matrix response not allowed")
  if (!missing(parms) && length(parms) > 0)
    warning("parameter argument ignored")
  if (length(offset))
    y <- y - offset
  sfun <- function(yval, dev, wt, ylevel, digits) {
    paste(" mean=",
          format(signif(yval, digits)),
          ", MSE=" ,
          format(signif(dev / wt, digits)),
          sep = '')
  }
  environment(sfun) <- .GlobalEnv
  list(
    y = c(y),
    parms = NULL,
    numresp = 1,
    numy = 1,
    summary = sfun
  )
}

# evaluation function
etemp <- function(y, wt, parms) {
  wmean <- sum(y * wt) / sum(wt)
  rss <- sum(wt * (y - wmean) ^ 2)
  list(label = wmean, deviance = rss)
}


# splitting function 
stemp <- function(y, wt, x, parms, continuous)
{
  # Center y
  n <- length(y)
  y <- y - sum(y * wt) / sum(wt)
  if (continuous) {
    # continuous x variable
    temp <- cumsum(y * wt)[-n]
    left.wt  <- cumsum(wt)[-n]
    right.wt <- sum(wt) - left.wt
    lmean <- temp / left.wt
    rmean <- -temp / right.wt
    goodness <- (left.wt * lmean ^ 2 + right.wt * rmean ^ 2) / sum(wt * y ^ 2)
    list(goodness = goodness, direction = sign(lmean))
  } else {
    # Categorical X variable
    ux <- sort(unique(x))
    wtsum <- tapply(wt, x, sum)
    ysum  <- tapply(y * wt, x, sum)
    means <- ysum / wtsum
    
    # For anova splits, we can order the categories by their means
    #  then use the same code as for a non-categorical
    ord <- order(means)
    n <- length(ord)
    temp <- cumsum(ysum[ord])[-n]
    left.wt  <- cumsum(wtsum[ord])[-n]
    right.wt <- sum(wt) - left.wt
    lmean <- temp / left.wt
    rmean <- -temp / right.wt
    list(
      goodness = (left.wt * lmean ^ 2 + right.wt * rmean ^ 2) / sum(wt * y ^ 2),
      direction = ux[ord]
    )
  }
}

library(rpart)
mystate <- data.frame(state.x77, region = state.region)
names(mystate) <- casefold(names(mystate)) # remove mixed case
ulist <- list(eval = etemp, split = stemp, init = itemp)
fit1 <- rpart(
  murder ~ population + illiteracy + income + life.exp +
    hs.grad + frost + region,
  data = mystate,
  method = ulist,
  minsplit = 10
)
fit2 <- rpart(
  murder ~ population + illiteracy + income + life.exp +
    hs.grad + frost + region,
  data = mystate,
  method = 'anova',
  minsplit = 10,
  xval = 0
) 

all.equal(fit1$frame, fit2$frame)

all.equal(fit1$splits, fit2$splits)

all.equal(fit1$csplit, fit2$csplit)

all.equal(fit1$where, fit2$where)

all.equal(fit1$cptable, fit2$cptable)

