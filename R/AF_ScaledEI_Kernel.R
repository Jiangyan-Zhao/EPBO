#' @title ScaledEI acquisition function
#' 
#' @description Scaled expected improvement
#' 
#' @param x description
#' 
#' @param fgpi description
#' 
#' @param fmean description
#' 
#' @param fsd description
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

AF_ScaledEI_Kernel = function(x, fgpi, fmean, fsd, Cgpi, epbest, rho, equal, type="UK")
{
  if(is.null(nrow(x))) x = matrix(x, nrow=1)
  ncand = nrow(x) # number of the candidate points
  
  ## objective
  pred_f = predict(object=fgpi, newdata=x, type=type, checkNames = FALSE, light.return = TRUE)
  mu_f = pred_f$mean * fsd + fmean
  sigma_f = pred_f$sd * fsd
  
  ## constraints
  nc = length(Cgpi) # number of the constraint
  mu_C = sigma_C = omega = matrix(NA, nc, ncand)
  for (j in 1:nc) {
    pred_C = predict(object=Cgpi[[j]], newdata=x, type=type, checkNames = FALSE, light.return = TRUE)
    mu_C[j,] = pred_C$mean
    sigma_C[j,] = pred_C$sd
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
  # ScaledEI[d <= -6] = 0
  return(ScaledEI)
}
