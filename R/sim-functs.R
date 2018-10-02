#' Update tumor mass, toxicity, and expected survival for the next stage
#'
#' Creates vectors of toxicities and tumor masses for the next month of
#' treatment as well as the rate constant for an exponential survival curve,
#' based off of the current tumor mass, toxiciity, dose, and patient
#' characteristics.
#'
#' See writeup for a full description
updateW <- function(M, W, D, a = 0.1, b = 1.2, c = 1, d = 0.5) {
  W_next <- a * M + b * (c * D - d) + W
  return(replace(W_next, W_next < 0, 0))
}

updateM <- function(M, W, D, a = 0.15, b = 1.2, c = 1, d = 0.5) {
  M_next <- a * W - b * (c * D - d) + M
  return(replace(M_next, M <= 0 | M_next < 0, 0))
}

lambda <- function(M, W, Z = 0, mu0 = -5.5, mu1 = 1, mu2 = 1.2, mu3 = 0.75) {
  return(exp(mu0 + mu1 * W + mu2 * M + mu3 * W * M + Z))
}
