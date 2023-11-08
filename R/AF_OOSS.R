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
#' 
#' @returns AF 
#' 
#' @author Jiangyan Zhao \email{zhaojy2017@126.com}
#' 
#' @references Pourmohamad, T. and H. K. H. Lee (2022). Bayesian optimization via barrier functions. 
#' \emph{Journal of Computational and Graphical Statistics} 31(1), 74-83.
#' 
#' @import laGP
#' 
#' 
#' 
#' @export
#' 
#' @examples 
#' B = rbind(c(0, 1), c(0, 1)) 
#' 
#' 

AF_OOSS = function(x, fgpi, fnorm, Cgpi, Cnorm)
{
  if(is.null(nrow(x))) x = matrix(x, nrow=1)
  ncand = nrow(x)
  nc = length(Cgpi) # number of the constraint
  
  ## objective part
  pred_f = predGPsep(fgpi, x, lite=TRUE)
  mu_f = pred_f$mean * fnorm
  sigma_f2 = pred_f$s2 * fnorm^2
  AF = - mu_f
  
  ## constraints part
  noValid = FALSE
  for (j in 1:nc) {
    pred_C = predGPsep(Cgpi[j], x, lite=TRUE)
    mu_C = pred_C$mean * Cnorm[j]
    if(all(mu_C >= 0)){
      noValid = TRUE
      break
    }
    sigma_C2 = pred_C$s2 * Cnorm[j]^2
    AF = AF + sigma_f2*(log(-mu_C) + sigma_C2/(2*mu_C^2))
  }
  
  ## 
  if(noValid || all(is.na(AF))){
    AF = 1
    for (j in 1:nc) {
      pred_C = predGPsep(Cgpi[j], x, lite=TRUE)
      mu_C = pred_C$mean * Cnorm[j]
      sigma_C = sqrt(pred_C$s2) * Cnorm[j]
      AF = AF * pnorm(0, mu_C, sigma_C)
    }
  }
  return(AF)
}