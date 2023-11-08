#' @title Un-normalize min-max function.
#' 
#' @description Un-normalizes X w.r.t. the provided bounds.
#' 
#' @param X \code{matrix} 
#' 
#' @param bounds 2-column \code{matrix} describing the bounding box. 
#' the first column gives lower bounds and the second gives upper bounds
#' 
#' @returns Un-normalized X w.r.t. the provided bounds.
#' If all elements of \code{X} are contained within \code{[0, 1]^d}, 
#' the un-normalized values will be contained within \code{bounds}.
#' 
#' @seealso \code{\link[EPBO]{normalize}}
#' 
#' @author Jiangyan Zhao \email{zhaojy2017@126.com}
#' 
#' 
#' @export
#' 
#' 
#' @examples 
#' X_normalized = matrix(runif(6), ncol=2)
#' bounds = cbind(c(-1,-1), c(1,1))
#' X = unnormalize(X_normalized, bounds)

unnormalize = function(X, bounds)
{
  if(is.null(nrow(X))) X = matrix(X, nrow=1)
  Xnorm = t(X) * (bounds[,2] - bounds[,1]) + bounds[,1]
  Xnorm = t(Xnorm)
  return(Xnorm)
}
