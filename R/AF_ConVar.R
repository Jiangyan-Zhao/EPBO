#' @title Constained variance acquisition function
#' 
#' @description Predictive mean
#' 
#' @param x description
#' 
#' @param fgpi description
#' 
#' @param Cgpi description
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

AF_ConVar = function(x, fgpi, Cgpi, equal)
{
  if(is.null(nrow(x))) x = matrix(x, nrow=1)
  ncand = nrow(x) # number of the candidate points
  
  ## objective
  pred_f = predGPsep(fgpi, x, lite=TRUE)
  sigma_f = sqrt(pred_f$s2)
  
  ## constraints
  nc = length(Cgpi) # number of the constraint
  LCB = matrix(NA, ncand, nc)
  for (j in 1:nc) {
    pred_C = predGPsep(Cgpi[j], x, lite=TRUE)
    mu_C = pred_C$mean
    sigma_C = sqrt(pred_C$s2)
    if(equal[j]){
      LCB[,j] = abs(mu_C) - 6*sigma_C
    }else{
      LCB[,j] = mu_C - 6*sigma_C
    }
  }
  
  ## Acquaisition function
  infeasible = apply(LCB > 0, 1, any)
  if(all(infeasible)){
    AF = sigma_f
  }else{
    sigma_f[infeasible] = 0
    AF = sigma_f
  }
  
  return(AF)
}
