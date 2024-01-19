#' @title Asymmetric entropy acquisition function
#' 
#' @description The asymmetric entropy (AE) acquisition function of the AE method
#' 
#' @param x a vector containing a single candidate point; or a \code{matrix} with 
#' multiple candidate points
#' @param fgpi the GP surrogate model of the objective function
#' @param fnorm the maximum of the objective 
#' @param Cgpi the GP surrogate models of the constraints
#' @param Cnorm the maxima of the constraints
#' @param fmin the best objective value obtained so far
#' @param alpha1 a specified weight for the EI part
#' @param alpha2 a specified weight for the AE part
#' @param omega a mode location parameter
#' 
#' @returns The AE value(s) at \code{x}. 
#' 
#' @seealso \code{\link[EPBO]{AF_ScaledEI}}, \code{\link[EPBO]{AF_EY}}, \code{\link[EPBO]{AF_OOSS}}
#' 
#' @author Jiangyan Zhao \email{zhaojy2017@126.com}
#' 
#' @references Lindberg, D. V. and H. K. Lee (2015). Optimization under constraints by applying an 
#' asymmetric entropy measure. \emph{Journal of Computational and Graphical Statistics} 24(2), 379-393.
#' 
#' @import laGP
#' @importFrom stats dnorm 
#' @importFrom stats pnorm 
#' @importFrom stats quantile


AF_AE = function(x, fgpi, fnorm, Cgpi, Cnorm, fmin,
                 alpha1=1, alpha2=5, omega=2/3)
{
  if(is.null(nrow(x))) x = matrix(x, nrow=1)
  ncand = nrow(x)
  nc = length(Cgpi) # number of the constraint
  
  ## calculate EI part
  pred_f = predGPsep(fgpi, x, lite=TRUE)
  mu_f = pred_f$mean * fnorm
  sigma_f = sqrt(pred_f$s2) * fnorm
  if(!is.finite(fmin)) fmin = quantile(mu_f, p=0.9) # adopt the recommendation of laGP package
  d = (fmin - mu_f)/sigma_f
  ei = sigma_f*(d*pnorm(d) + dnorm(d))
  
  ## asymmetric entropy
  PoF = rep(1, ncand)
  for (j in 1:nc) {
    pred_C = predGPsep(Cgpi[j], x, lite=TRUE)
    mu_C = pred_C$mean * Cnorm[j]
    sigma_C = sqrt(pred_C$s2) * Cnorm[j]
    PoF = PoF * pnorm(0, mu_C, sigma_C)
  }
  Sa = 2*PoF*(1-PoF)/(PoF-2*omega*PoF+omega^2)
  
  ## acquisition function
  AF = ei^alpha1 * Sa^alpha2
  return(AF)
}