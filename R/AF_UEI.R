#' @title UEI acquisition function
#' 
#' @description Scaled expected improvement
#' 
#' @param x description
#' 
#' @param fgpi description
#' 
#' @param Cgpi description
#' 
#' @param epbest description
#' 
#' @param rho description
#' 
#' @param equal an optional vector containing zeros and ones, whose length equals the number of
#' constraints, specifying which should be treated as equality constraints (\code{1}) and 
#' which as inequality (\code{0}) 
#' 
#' 
#' @returns AF 
#' 
#' @author Jiangyan Zhao \email{zhaojy2017@126.com}
#' 
#' @references Noe, U. and D. Husmeier (2018). On a new improvement-based acquisition function 
#' for Bayesian optimization. \emph{arXiv:1808.06918}.
#' 
#' @import laGP
#' @importFrom stats dnorm 
#' @importFrom stats pnorm 
#' 
#' 
#' @export
#' 
#' @examples 
#' B = rbind(c(0, 1), c(0, 1)) 
#' 
#' 

AF_UEI = function(x, fgpi, Cgpi, epbest, rho, equal)
{
  if(is.null(nrow(x))) x = matrix(x, nrow=1)
  ncand = nrow(x) # number of the candidate points
  
  ## objective
  pred_f = predGPsep(fgpi, x, lite=TRUE)
  mu_f = pred_f$mean
  sigma_f = sqrt(pred_f$s2)
  
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
  EI = sigma_ep * (d*pnorm(d) + dnorm(d)) # expected improvement
  VI = sigma_ep^2 * ((d^2+1)*pnorm(d) + d*dnorm(d)) - EI^2 # variance of the improvement (remove sigma_ep)
  VI = pmax(0, VI)
  UEI = EI + 2*sqrt(VI) # Scaled expected improvement
  return(UEI)
}
