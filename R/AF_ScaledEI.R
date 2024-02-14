#' @title ScaledEI acquisition function
#' 
#' @description The Scaled expected improvement acquisition function of the EPBO method
#' 
#' @param x a vector containing a single candidate point; or a \code{matrix} with 
#' multiple candidate points
#' @param fgpi the GP surrogate model of the objective function
#' @param fmean the mean of the objective value
#' @param fsd the standard deviation of the objective value
#' @param Cgpi the GP surrogate models of the constraints
#' @param epbest the best exact penalty value obtained so far
#' @param rho the penalty parameters
#' @param equal an optional vector containing zeros and ones, whose length equals the number of
#' constraints, specifying which should be treated as equality constraints (\code{1}) and 
#' which as inequality (\code{0}) 
#' 
#' @returns The ScaledEI value(s) at \code{x}.
#' 
#' @seealso \code{\link[EPBO]{AF_OOSS}}, \code{\link[EPBO]{AF_EY}}, \code{\link[EPBO]{AF_AE}}
#' 
#' @author Jiangyan Zhao \email{zhaojy2017@126.com}
#' 
#' @references Noe, U. and D. Husmeier (2018). On a new improvement-based acquisition function 
#' for Bayesian optimization. \emph{arXiv:1808.06918}.
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

AF_ScaledEI = function(x, fgpi, fmean, fsd, Cgpi, epbest, rho, equal)
{
  if(is.null(nrow(x))) x = matrix(x, nrow=1)
  ncand = nrow(x) # number of the candidate points
  
  ## objective
  pred_f = predGPsep(fgpi, x, lite=TRUE)
  mu_f = pred_f$mean * fsd + fmean
  sigma_f = sqrt(pred_f$s2) * fsd
  
  ## constraints
  nc = length(Cgpi) # number of the constraint
  mu_C = sigma_C = omega = matrix(NA, nc, ncand)
  for (j in 1:nc) {
    pred_C = predGPsep(Cgpi[j], x, lite=TRUE)
    mu_C[j,] = pred_C$mean
    sigma_C[j,] = sqrt(pred_C$s2)
    omega[j,] = (equal[j]+1)*pnorm(mu_C[j,]/sigma_C[j,]) - equal[j]
  }
  
  ## Acquaisition function
  mu_ep = mu_f + rho%*%(omega*mu_C)
  sigma_ep = sqrt(sigma_f^2 + (rho^2)%*%((omega*sigma_C)^2))
  d = (epbest - mu_ep)/sigma_ep
  EI = d*pnorm(d) + dnorm(d) # expected improvement (remove sigma_ep)
  VI = pnorm(d) + d*EI - EI^2 # variance of the improvement (remove sigma_ep)
  VI = pmax(.Machine$double.xmin, VI)
  ScaledEI = EI/sqrt(VI) # Scaled expected improvement
  # ScaledEI[d < -6] = 0
  ScaledEI[is.nan(ScaledEI)] = 0
  return(ScaledEI)
}
