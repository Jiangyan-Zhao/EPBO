#' @title AHP: Analytic Hierarchy Process
#' @description AHP is a multi-criteria decision analysis method developed
#' by Saaty, which can also be used to
#' determine indicator weights.
#' @param A a numeric matrix, i.e. pairwise comparison matrix
#' @return a list object that contains: W (Weight vector), CR (Consistency ratio),
#' Lmax (Maximum eigenvalue), CI (Consistency index)
#' @export
#' @examples
#' A = matrix(c(1, 1/2, 4, 3, 3,
#' 2, 1, 7, 5, 5,
#' 1/4, 1/7, 1, 1/2, 1/3,
#' 1/3, 1/5, 2, 1, 1,
#' 1/3, 1/5, 3, 1, 1), byrow = TRUE, nrow = 5)
#' AHP(A)
AHP <- function(A) {
  rlt <- eigen(A)
  Lmax <- Re(rlt$values[1]) # Maximum eigenvalue
  # Weight vector
  W <- Re(rlt$vectors[,1]) / sum(Re(rlt$vectors[,1]))
  # Consistency index
  n <- nrow(A)
  CI <- (Lmax - n) / (n - 1)
  # Consistency ratio
  # Saaty's random Consistency indexes
  RI <- c(0,0,0.58,0.90,1.12,1.24,1.32,1.41,1.45,1.49,1.51)
  CR <- CI / RI[n]
  list(W = W, CR = CR, Lmax = Lmax, CI = CI)
}