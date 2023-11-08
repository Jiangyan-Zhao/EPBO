#' @title ScaledEI acquisition function
#' 
#' @description Scaled expected improvement
#' 
#' @param x description
#' 
#' @param fgpi description
#' 
#' @param Cgpi description
#' 
#' @param rho description
#' 
#' @param equal an optional vector containing zeros and ones, whose length equals the number of
#' constraints, specifying which should be treated as equality constraints (\code{1}) and 
#' which as inequality (\code{0}) 
#' 
#' @param batch_size aaa
#' 
#' @returns ScaledEI 
#' 
#' @author Jiangyan Zhao \email{zhaojy2017@126.com}
#' 
#' @references Noe, U. and D. Husmeier (2018). On a new improvement-based acquisition function 
#' for Bayesian optimization. \emph{arXiv:1808.06918}.
#' 
#' @import laGP
#' @importFrom stats rnorm
#' 
#' 
#' @export
#' 
#' @examples 
#' B = rbind(c(0, 1), c(0, 1)) 
#' 
#' 

AF_TS = function(x, fgpi, Cgpi, rho, equal, batch_size=1)
{
  if(is.null(nrow(x))) x = matrix(x, nrow=1)
  ncand = nrow(x) # number of the candidate points
  
  ## objective
  pred_f = predGPsep(fgpi, x)
  mu_f = pred_f$mean
  sigma_f = pred_f$Sigma
  df_f = pred_f$df
  
  ## constraints
  nc = length(Cgpi) # number of the constraint
  mu_C = matrix(NA, nc, ncand)
  sigma_C = vector("list", length = nc)
  df_C = rep(NA, nc)
  for (j in 1:nc) {
    pred_C = predGPsep(Cgpi[j], x)
    mu_C[j,] = pred_C$mean
    sigma_C[[j]] = pred_C$Sigma
    df_C[j] = pred_C$df
  }
  
  ## expected improvement
  AF = matrix(NA, batch_size, ncand)
  for (b in 1:batch_size) {
    AF[b,] = rmvt(1, sigma_f, df_f) + mu_f
    for (j in 1:nc) {
      if(equal[j]){
        AF[b,] = AF[b,] + rho[j]*abs(rmvt(1, sigma_C[[j]], df_C[j])+mu_C[j,])
      }else{
        AF[b,] = AF[b,] + rho[j]*pmax(0, rmvt(1, sigma_C[[j]], df_C[j])+mu_C[j,])
      }
    }
  }
  return(AF)
}
