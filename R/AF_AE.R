#' @title Asymmetric Entropy
#' 
#' @description Asymmetric Entropy
#' 
#' @param x description
#' 
#' @param fgpi description
#' 
#' @param fnorm description
#' 
#' @param Cgpi description
#' 
#' @param Cnorm description
#' 
#' @param fmin description
#' 
#' @param alpha1 description
#' 
#' @param alpha2 description
#' 
#' @param omega description
#' 
#' @returns AF 
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
#' 
#' 
#' 
#' @export
#' 
#' @examples 
#' B = rbind(c(0, 1), c(0, 1)) 
#' 
#' 

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
  PoF = matrix(NA, nrow = nrow(x), ncol = nc)
  for (j in 1:nc) {
    pred_C = predGPsep(Cgpi[j], x, lite=TRUE)
    mu_C = pred_C$mean * Cnorm[j]
    sigma_C = sqrt(pred_C$s2) * Cnorm[j]
    PoF[,j] = pnorm(0, mu_C, sigma_C)
  }
  p = apply(PoF, 1, prod)
  Sa = 2*p*(1-p)/(p-2*omega*p+omega^2)
  
  ## acquisition function
  AF = ei^alpha1 * Sa^alpha2
  return(AF)
}