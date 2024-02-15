#' @title min-max normalization
#' 
#' @description min-max normalization
#' 
#' @param x a \code{matrix} or \code{vector}
#' 
#' @returns The scaled values.
#'
#' @author Jiangyan Zhao \email{zhaojy2017@126.com}
#' 
#' @export

min_max_scale = function(x, byrow = FALSE){
  if(byrow){
    if(is.null(nrow(x))) x = matrix(x, nrow=1)
    x_scaled = apply(x, 1, function(x){(x-min(x))/(max(x)-min(x))})
  }else{
    if(is.null(ncol(x))) x = matrix(x, ncol=1)
    x_scaled = apply(x, 2, function(x){(x-min(x))/(max(x)-min(x))})
  }
  return(x_scaled)
}