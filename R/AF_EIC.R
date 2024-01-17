#' @title EI versus PoF
#' 
#' @description EI versus PoF
#' 
#' @param x a vector containing a single candidate point
#' @param fgpi the GP surrogate model of the objective function
#' @param fmean the mean of the objective value
#' @param fsd the standard deviation of the objective value
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


AF_EIC = function(x, fgpi, fmean, fsd, Cgpi, fmin)
{
  if(!is.null(nrow(x))) stop("x must be a vector, that is, one test point at a time")
  x = matrix(x, nrow=1)
  
  ## the EI part
  pred_f = predGPsep(fgpi, x, lite=TRUE)
  mu_f = pred_f$mean * fsd + fmean
  sigma_f = sqrt(pred_f$s2) * fsd
  if(!is.finite(fmin)) fmin = quantile(mu_f, p=0.9) # adopt the recommendation of laGP package
  d = (fmin - mu_f)/sigma_f
  EI = sigma_f*(d*pnorm(d) + dnorm(d))
  
  ## the PoF part
  PoF = 1
  for (j in 1:length(Cgpi)) {
    pred_C = predGPsep(Cgpi[j], x, lite=TRUE)
    mu_C = pred_C$mean
    sigma_C = sqrt(pred_C$s2)
    PoF = PoF*pnorm(0, mu_C, sigma_C)
  }
  
  ## EI versus PoF 
  AF = -c(log(EI), PoF) # minimize
  return(AF)
}