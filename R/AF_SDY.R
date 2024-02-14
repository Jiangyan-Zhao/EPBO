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
#' @keywords acquistionFunction
#' 
#' @import laGP
#' @importFrom stats dnorm 
#' @importFrom stats pnorm 
#' 
#' @export


AF_SDY = function(x, fgpi, fmean, fsd, Cgpi, rho, equal)
{
  if(is.null(nrow(x))) x = matrix(x, nrow=1)
  ncand = nrow(x) # number of the candidate points
  
  ## objective
  pred_f = predGPsep(fgpi, x, lite=TRUE)
  # mu_f = pred_f$mean * fsd + fmean
  # sigma_f = sqrt(pred_f$s2) * fsd
  
  ## constraints
  nc = length(Cgpi) # number of the constraint
  EV = VV = matrix(NA, nc, ncand)
  for (j in 1:nc) {
    pred_C = predGPsep(Cgpi[j], x, lite=TRUE)
    mu_C = pred_C$mean
    sigma_C = sqrt(pred_C$s2)
    dC = mu_C/sigma_C
    dC_cdf = pnorm(dC)
    dC_pdf = dnorm(dC)
    if(equal[j]){ 
      EV[j,] = mu_C * (2 * dC_cdf - 1) + 2 * sigma_C * dC_pdf
      VV[j,] = (mu_C^2 + sigma_C^2) - EV[j,]^2
    }else{
      EV[j,] = mu_C * dC_cdf + sigma_C * dC_pdf
      VV[j,] = (mu_C^2 + sigma_C^2) * dC_cdf + mu_C * sigma_C * dC_pdf- EV[j,]^2
    }
  }
  VV = pmax(0, VV)
  
  ## the predictive standard deviation of the exact penalty surrogate
  SDY = sqrt(pred_f$s2 + rho^2 %*% VV)
  SDY = pmax(.Machine$double.xmin, SDY)
  
  return(log(SDY)) # log scale
}
