#' @title EI versus PoF
#' 
#' @description EI versus PoF
#' 
#' @param x a vector containing a single candidate point
#' @param fgpi the GP surrogate model of the objective function
#' @param Cgpi the GP surrogate models of the constraints
#' @param fmin the best feasible objective value obtained so far
#' 
#' @returns The AE value(s) at \code{x}. 
#' 
#' @seealso \code{\link[EPBO]{AF_ScaledEI}}, \code{\link[EPBO]{AF_EY}}, \code{\link[EPBO]{AF_OOSS}}
#' 
#' @author Jiangyan Zhao \email{zhaojy2017@126.com}
#' 
#' @references Parr, J. M., A. J. Keane, A. I. Forrester, and C. M. Holden (2012). Infill sampling criteria
#' for surrogate-based optimization with constraint handling. \emph{Engineering Optimization} 44(10), 1147–1166.
#' 
#' @import laGP
#' @importFrom stats dnorm 
#' @importFrom stats pnorm 
#' @importFrom stats quantile
#' 
#' @export


AF_EIvsPoF = function(x, fgpi, Cgpi, fmin)
{
  if(!is.null(nrow(x))) stop("x must be a vector, that is, one test point at a time")
  x = matrix(x, nrow=1)
  
  ## the EI part
  pred_f = predGPsep(fgpi, x, lite=TRUE)
  mu_f = pred_f$mean
  sigma_f = sqrt(pred_f$s2)
  if(!is.finite(fmin)) fmin = 1e6
  d = (fmin - mu_f)/sigma_f
  EI = sigma_f*(d*pnorm(d) + dnorm(d))
  EI = max(EI, .Machine$double.eps)
  logEI = log(EI)
  
  ## the PoF part
  nc = length(Cgpi) # number of the constraint
  # logPoF = 0
  PoF = 1
  for (j in 1:nc) {
    pred_C = predGPsep(Cgpi[j], x, lite=TRUE)
    mu_C = pred_C$mean
    sigma_C = sqrt(pred_C$s2)
    # logPoF = logPoF + log(pnorm(0, mu_C, sigma_C)) # (log scale)
    PoF = PoF * pnorm(0, mu_C, sigma_C)
  }
  # logPoF = max(logPoF, log(.Machine$double.eps))
  

  ## EI versus PoF 
  # AF = c(logEI, logPoF)
  AF = c(logEI, PoF)
  
  return(-AF) # minimize
}