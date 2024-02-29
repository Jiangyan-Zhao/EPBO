#' @title PEIC
#' 
#' @description PEIC
#' 
#' @param x a vector containing a single candidate point
#' @param fgpi the GP surrogate model of the objective function
#' @param Cgpi the GP surrogate models of the constraints
#' @param fmin the best feasible objective value obtained so far
#' @param df the lengthscale parameters of the \code{fgpi} surrogate model
#' @param point_update the updating points in a single iteration in parallel
#' 
#' @returns The AE value(s) at \code{x}. 
#' 
#' @seealso \code{\link[EPBO]{AF_ScaledEI}}, \code{\link[EPBO]{AF_EY}}, \code{\link[EPBO]{AF_OOSS}}
#' 
#' @author Jiangyan Zhao \email{zhaojy2017@126.com}
#' 
#' @references Qian, J., Y. Cheng, J. Zhang, J. Liu, and D. Zhan (2021). A parallel constrained efficient global 
#' optimization algorithm for expensive constrained optimization problems. \emph{Engineering Optimization} 53(2), 300–320.
#' 
#' @import laGP
#' @importFrom plgp covar.sep
#' @importFrom stats dnorm 
#' @importFrom stats pnorm 
#' @importFrom stats quantile
#' 
#' @export


AF_PEIC = function(x, fgpi, Cgpi, fmin, df, point_update)
{
  if(is.null(nrow(x))) x = matrix(x, nrow=1)
  ncand = nrow(x)
  nc = length(Cgpi) # number of the constraint
  
  ## PoF 
  PoF = rep(1, ncand)
  for (j in 1:nc) {
    pred_C = predGPsep(Cgpi[j], x, lite=TRUE)
    mu_C = pred_C$mean
    sigma_C = sqrt(pred_C$s2)
    PoF = PoF * pnorm(0, mu_C, sigma_C)
  }
  
  ## EI 
  EI = rep(1, ncand)
  if(is.finite(fmin)){
    pred_f = predGPsep(fgpi, x, lite=TRUE)
    mu_f = pred_f$mean
    sigma_f = sqrt(pred_f$s2)
    d = (fmin - mu_f)/sigma_f
    EI = sigma_f*(d*pnorm(d) + dnorm(d))
  }
  
  ## EIC (=PoF if there is no feasible point in the design)
  EIC = EI * PoF
  
  ## Influence function
  IF = rep(1, ncand) # the first updating point is selected by the standard EIC
  if(nrow(point_update) > 0){
    corr = covar.sep(x, point_update, d=df, g=1e-6) # correlation function
    IF = apply(1 - corr, 1, prod)
  }
  
  ## PEIC (=PPoF if there is no feasible point in the design)
  PEIC = EIC * IF
  
  return(PEIC)
}