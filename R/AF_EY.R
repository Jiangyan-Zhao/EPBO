#' @title Predictive mean acquisition function
#' 
#' @description The predictive mean acquisition function of the (P)EPBO method
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
#' @returns The Predictive Mean at \code{x}. 
#' 
#' @seealso \code{\link[EPBO]{AF_ScaledEI}}, \code{\link[EPBO]{AF_OOSS}}, \code{\link[EPBO]{AF_AE}}
#' 
#' @author Jiangyan Zhao \email{zhaojy2017@126.com}
#' 
#' @references Jiangyan Zhao and Jin Xu (07 Feb 2024): Bayesian Optimization via Exact Penalty, 
#' \emph{Technometrics}, DOI: 10.1080/00401706.2024.2315937
#' 
#' @keywords acquistionFunction
#' 
#' @import laGP
#' @importFrom stats dnorm 
#' @importFrom stats pnorm 
#' 
#' @export


AF_EY = function(x, fgpi, fmean, fsd, Cgpi, rho, equal)
{
  if(is.null(nrow(x))) x = matrix(x, nrow=1)
  ncand = nrow(x) # number of the candidate points

  ## objective
  pred_f = predGPsep(fgpi, x, lite=TRUE)
  mu_f = pred_f$mean * fsd + fmean
  
  ## constraints
  nc = length(Cgpi) # number of the constraint
  EV = matrix(NA, nc, ncand)
  for (j in 1:nc) {
    pred_C = predGPsep(Cgpi[j], x, lite=TRUE)
    mu_C = pred_C$mean
    sigma_C = sqrt(pred_C$s2)
    dC = mu_C/sigma_C
    EV[j,] = mu_C*((equal[j]+1)*pnorm(dC) - equal[j]) + (equal[j]+1)*sigma_C*dnorm(dC)
  }
  
  ## the predictive mean of the exact penalty surrogate
  EY = mu_f + rho%*%EV
  
  return(EY)
}
