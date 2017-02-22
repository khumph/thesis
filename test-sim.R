library(testthat)
library(tidyverse)
library(magrittr)

updateW <- function(M, W, D, a1 = 0.1, b1 = 1.2, d1 = 0.5) {
  W_next <- a1 * M + b1 * (D - d1) # + W
  ifelse(W_next > 0, W_next, 0)
}


tdat <- expand.grid(M = seq(0, 6, 0.25),
                    W = seq(0, 6, 0.25),
                    D = seq(0, 1, 0.05)) %>% mutate(
                      W_next = updateW(M, W, D)
                    )

ggplot(tdat) + geom_raster(aes(x = M, y = D, fill = W_next))


test_that("updateW ", {
  expect_equal(updateW(0, 0, 0), -0.6)
  expect_equal(updateW(1, 1, 1), 0.1 * 1 + 1.2 * (1 - 0.5) + 1)
  expect_equal(str_length("abc"), 3)
})

test_that("str_length of factor is length of level", {
  expect_equal(str_length(factor("a")), 1)
  expect_equal(str_length(factor("ab")), 2)
  expect_equal(str_length(factor("abc")), 3)
})

test_that("str_length of missing is missing", {
  expect_equal(str_length(NA), NA_integer_)
  expect_equal(str_length(c(NA, 1)), c(NA, 1))
  expect_equal(str_length("NA"), 2)
})

test_that("", {
  expect_equal(str_length(NA), NA_integer_)
  expect_equal(str_length(c(NA, 1)), c(NA, 1))
  expect_equal(str_length("NA"), 2)
})