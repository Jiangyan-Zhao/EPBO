#' @title Predictive standard deviation acquisition function
#' 
#' @description The predictive standard deviation acquisition function of the PEPBO method
#' 
#' @param x a vector containing a single candidate point; or a \code{matrix} with 
#' multiple candidate points
#' @param fgpi the GP surrogate model of the objective function
#' @param fmean the mean of the objective value
#' @param fsd the standard deviation of the objective value
#' @param Cgpi the GP surrogate models of the constraints
#' @param rho the penalty parameters
#' @param equal an optional vector containing zeros and ones, whose length equals the number of
#' constraints, specifying which should be treated as equality constraints (\code{1}) and 
#' which as inequality (\code{0}) 
#' 
#' @returns The Predictive Standard Deviation at \code{x}. 
#' 
#' @seealso \code{\link[EPBO]{AF_ScaledEI}}, \code{\link[EPBO]{AF_OOSS}}, \code{\link[EPBO]{AF_AE}}
#' 
#' @author Jiangyan Zhao \email{zhaojy2017@126.com}
#' 
#' @import laGP
#' @importFrom stats dnorm 
#' @importFrom stats pnorm 


AF_EYvsSDY = function(x, fgpi, fmean, fsd, Cgpi, rho, equal)
{
  if(!is.null(nrow(x))) stop("x must be a vector, that is, one test point at a time")
  x = matrix(x, nrow=1)
  
  ## the predictive mean of the exact penalty surrogate
  EY = AF_EY(x, fgpi, fmean, fsd, Cgpi, rho, equal)
  
  ## the predictive standard deviation of the exact penalty surrogate
  logSDY = AF_SDY(x, fgpi, fmean, fsd, Cgpi, rho, equal) # log scale
  
  ## both criterion are minimized
  AF = c(-logSDY, EY)
  
  return(AF) 
}
