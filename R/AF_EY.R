#' @title Predictive mean acquisition function
#' 
#' @description Predictive mean
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

AF_EY = function(x, fgpi, fmean, fsd, Cgpi, rho, equal)
{
  if(is.null(nrow(x))) x = matrix(x, nrow=1)
  ncand = nrow(x) # number of the candidate points

  ## objective
  pred_f = predGPsep(fgpi, x, lite=TRUE)
  mu_f = pred_f$mean * fsd + fmean
  
  ## constraints
  nc = length(Cgpi) # number of the constraint
  EV = matrix(NA, nc, nrow(x))
  for (j in 1:nc) {
    pred_C = predGPsep(Cgpi[j], x, lite=TRUE)
    mu_C = pred_C$mean
    sigma_C = sqrt(pred_C$s2)
    dC = mu_C/sigma_C
    EV[j,] = mu_C*((equal[j]+1)*pnorm(dC) - equal[j]) + (equal[j]+1)*sigma_C*dnorm(dC)
  }
  
  ## Acquaisition function
  AF = mu_f + rho%*%EV
  
  return(AF)
}