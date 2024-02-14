#' @title PoF
#' 
#' @description PoF
#' 
#' @param x a vector containing a single candidate point
#' @param Cgpi the GP surrogate models of the constraints
#' 
#' @returns The PoF value(s) at \code{x}. 
#' 
#' @seealso \code{\link[EPBO]{AF_ScaledEI}}, \code{\link[EPBO]{AF_EY}}, \code{\link[EPBO]{AF_OOSS}}
#' 
#' @author Jiangyan Zhao \email{zhaojy2017@126.com}
#' 
#' @import laGP
#' @importFrom stats pnorm 
#' 
#' @export


AF_PoF = function(x, Cgpi)
{
  if(is.null(nrow(x))) x = matrix(x, nrow=1)
  ncand = nrow(x)
  nc = length(Cgpi) # number of the constraint
  
  ## PoF 
  PoF = rep(1, ncand)
  for (j in 1:nc) {
    pred_C = predGPsep(Cgpi[j], x, lite=TRUE)
    mu_C = pred_C$mean
    sigma_C = sqrt(pred_C$s2)
    PoF = PoF * pnorm(0, mu_C, sigma_C)
  }
  
  return(PoF)
}