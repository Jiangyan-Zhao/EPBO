#' @title Min-max normalize and un-normalize function
#' 
#' @description Min-max normalize X w.r.t. the provided bounds.
#' 
#' @param X \code{matrix} the candidate points
#' 
#' @param bounds 2-column \code{matrix} describing the bounding box. 
#' the first column gives lower bounds and the second gives upper bounds
#' 
#' @returns Min-max normalize or un-normalized X w.r.t. the provided bounds.
#' The normalized values will be contained within \code{[0, 1]^d}, 
#' if all elements of \code{X} are contained within \code{bounds}.
#' the un-normalized values will be contained within \code{bounds},
#' if all elements of \code{X} are contained within \code{[0, 1]^d}, 
#' 
#' @author Jiangyan Zhao \email{zhaojy2017@126.com}
#' 
#' @export


normalize = function(X, bounds)
{
  if(is.null(nrow(X))) X = matrix(X, nrow=1)
  Xnorm = (t(X) - bounds[,1]) / (bounds[,2] - bounds[,1])
  Xnorm = t(Xnorm)
  return(Xnorm)
}

unnormalize = function(X, bounds)
{
  if(is.null(nrow(X))) X = matrix(X, nrow=1)
  Xnorm = t(X) * (bounds[,2] - bounds[,1]) + bounds[,1]
  Xnorm = t(Xnorm)
  return(Xnorm)
}
