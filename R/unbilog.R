#' @title Un-transform bilog-transformed outcomes
#' 
#' @description
#' The Bilog transform is useful for modeling outcome constraints 
#' as it magnifies values near zero and flattens extreme values.
#' 
#' @param Y A bilog-transfomred targets
#' 
#' @returns The un-transformed outcome observations.
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

unbilog = function(Y)
{
  Y_utf = sign(Y) * (exp(abs(Y)) - 1)
  return(Y_utf)
}









