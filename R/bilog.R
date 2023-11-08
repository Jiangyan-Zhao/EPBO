#' @title Bilog-transform outcomes
#' 
#' @description
#' The Bilog transform is useful for modeling outcome constraints 
#' as it magnifies values near zero and flattens extreme values.
#' 
#' @param Y A training targets
#' 
#' @returns The transformed outcome observations.
#' 
#' @seealso \code{\link[EPBO]{bilog}}
#'
#' @author Jiangyan Zhao \email{zhaojy2017@126.com}
#' 
#' @references Eriksson, D., & Poloczek, M. (2021, March). Scalable constrained Bayesian optimization. 
#' In \emph{International Conference on Artificial Intelligence and Statistics} (pp. 730-738). PMLR.
#' 
#' @export
#' 

bilog = function(Y)
{
  Y_tf = sign(Y) * log(abs(Y) + 1)
  return(Y_tf)
}