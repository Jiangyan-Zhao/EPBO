#' @title UEI acquisition function
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
#' @param N description
#' 
#' @param beta description
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

AF_UEI_MC = function(x, fgpi, fmean, fsd, Cgpi, epbest, rho, equal, N=1000, beta=3)
{
  if(is.null(nrow(x))) x = matrix(x, nrow=1)
  ncand = nrow(x) # number of the candidate points
  
  ## objective
  pred_f = predGPsep(fgpi, x, lite=TRUE)
  mu_f = pred_f$mean * fsd + fmean
  sigma_f = sqrt(pred_f$s2) * fsd
  
  ## constraints
  nc = length(Cgpi) # number of the constraint
  mu_C = sigma_C = matrix(NA, nc, ncand)
  for (j in 1:nc) {
    pred_C = predGPsep(Cgpi[j], x, lite=TRUE)
    mu_C[j,] = pred_C$mean
    sigma_C[j,] = sqrt(pred_C$s2)
  }
  
  ## expected improvement
  EP = matrix(NA, nrow = N, ncol = ncand)
  for (n in 1:N) {
    EP[n,] = rnorm(ncand, mu_f, sigma_f)
    for (j in 1:nc) {
      if(equal[j]){
        EP[n,] = EP[n,] + rho[j]*abs(rnorm(ncand, mu_C[j,], sigma_C[j,]))
      }else{
        EP[n,] = EP[n,] + rho[j]*pmax(0, rnorm(ncand, mu_C[j,], sigma_C[j,]))
      }
    }
  }
  improvement = matrix(pmax(0, epbest-EP), nrow = N)
  EI = colMeans(improvement) # expected improvement
  SI = colSds(improvement)   # Standard deviation of the improvement
  UEI = EI + beta * SI
  EI[is.nan(UEI)] = 0
  return(UEI)
}
