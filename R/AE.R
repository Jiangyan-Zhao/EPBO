#' @title Asymmetric Entropy
#' 
#' @description Black-box optimization under inequality constraints via an asymmetric entropy.
#' 
#' @param blackbox blackbox of an input (\code{x}), facilitating vectorization on a 
#' \code{matrix} \code{X} thereof,  returning a \code{list} 
#' with elements \code{$obj} containing the (scalar) objective value and \code{$c} 
#' containing a vector of evaluations of the (multiple) constraint function at \code{x}.
#' @param B 2-column \code{matrix} describing the bounding box.  The number of rows
#' of the \code{matrix} determines the input dimension (\code{length(x)} in \code{blackbox(x)}); 
#' the first column gives lower bounds and the second gives upper bounds
#' @param Xstart optional matrix of starting design locations in lieu of, or in addition to,
#' \code{start} random ones;  we recommend \code{nrow(Xstart) + start >= 6}; also must
#' have \code{ncol(Xstart) = nrow(B)}
#' @param start  positive integer giving the number of random starting locations before 
#' sequential design (for optimization) is performed; \code{start >= 6} is
#' recommended unless \code{Xstart} is non-\code{NULL}; in the current version
#' the starting locations come from a space-filling design via \code{\link[tgp]{dopt.gp}}
#' @param end positive integer giving the total number of evaluations/trials in the 
#' optimization; must have \code{end > start}
#' @param alpha1 weight for expected improvement (default: 1) 
#' @param alpha2 weight for asymmetric entropy (default: 5)  
#' @param omega mode location parameter (default: 2/3) 
#' @param urate positive integer indicating  how many optimization trials should pass before
#' each MLE/MAP update is performed for GP correlation lengthscale 
#' parameter(s) 
#' @param ncandf function taking a single integer indicating the optimization trial number \code{t}, where
#' \code{start < t <= end}, and returning the number of search candidates (e.g., for
#' expected improvement calculations) at round \code{t}; the default setting
#' allows the number of candidates to grow linearly with \code{t}
#' @param dg.start 2-vector giving starting values for the lengthscale and nugget parameters
#' of the GP surrogate model(s) for constraints
#' @param ab prior parameters; see \code{\link{darg}} describing the prior used on the
#' lengthscale parameter during emulation(s) for the constraints
#' @param dlim 2-vector giving bounds for the lengthscale parameter(s) under MLE/MAP inference,
#' thereby augmenting the prior specification in \code{ab}
#' @param verb a non-negative integer indicating the verbosity level; the larger the value the
#' more that is printed to the screen
#' @param ... additional arguments passed to \code{blackbox}
#'
#' @returns The output is a \code{list} summarizing the progress of the evaluations of the 
#' blackbox under optimization:
#' \item{prog }{ vector giving the best feasible (\code{g(x) <= 0 && |h(x)| <= ethresh}) value of the objective over the trials }
#' \item{xbest }{ vector giving the recommended solution}
#' \item{obj }{ vector giving the value of the objective for the input under consideration at each trial }
#' \item{C }{ \code{matrix} giving the value of the constraint function for the input under consideration at each trial}
#' \item{X }{ \code{matrix} giving the input values at which the blackbox function was evaluated }
#' 
#' @seealso \code{\link[laGP]{optim.auglag}}, \code{\link[laGP]{optim.efi}}, \code{\link[EPBO]{optim.BM}}, \code{\link[EPBO]{optim.EP}}
#'  
#' @author Jiangyan Zhao \email{zhaojy2017@126.com}
#' 
#' @references Lindberg, D. V. and H. K. Lee (2015). Optimization under constraints by applying an 
#' asymmetric entropy measure. \emph{Journal of Computational and Graphical Statistics} 24(2), 379-393.
#' 
#' @keywords optimize
#' @keywords design
#' 
#' @import laGP
#' @import tgp
#' @importFrom stats dnorm 
#' @importFrom stats optim 
#' @importFrom stats pnorm 
#' @importFrom stats sd
#' @importFrom utils tail
#' 
#' @export
#'
#' @examples 
#' # search space
#' B = rbind(c(-2.25, 2.5), c(-2.5, 1.75))
#' # Objective and constraint function for use 
#' MTP = function(x){
#'   f = -(cos((x[1]-0.1)*x[2]))^2 - x[1]*sin(3*x[1]+x[2])
#'   t = atan2(x[1], x[2])
#'   c = x[1]^2 + x[2]^2 -((2*cos(t)-1/2*cos(2*t)-1/4*cos(3*t)-1/8*cos(4*t))^2) - ((2*sin(t))^2)
#'   return(list(obj=f, c=c))
#' }
#' AE = optim.AE(MTP, B, ncandf = function(k){ 1e3 }, start = 20, end = 120, verb = 0)
#' # progress, best feasible value of the objective over the trials
#' AE$prog
#' # the recommended solution
#' AE$xbest


optim.AE = function(
    blackbox, B, start=10, end=100, 
    Xstart=NULL, urate=10, ncandf=function(k) { k }, 
    alpha1=1, alpha2=5, omega=2/3,
    dg.start=c(0.1,1e-6), ab=c(3/2,8), dlim=sqrt(ncol(B))*c(1/100,10), 
    verb=2, ...)
{
  ## check start
  if(start >= end) stop("must have start < end")
  ## get initial design
  X = dopt.gp(start, Xcand=lhs(10*start, B))$XX
  X = rbind(Xstart, X)
  start = nrow(X)
  
  ## first run to determine dimensionality of the constraint
  out = blackbox(X[1,], ...)
  nc = length(out$c)
  
  ## allocate progress objects, and initialize
  prog = obj = rep(NA, start)
  C = matrix(NA, nrow=start, ncol=nc)
  obj[1] = out$obj; C[1,] = out$c
  prog[1] = ifelse(all(C[1,] <= 0), out$obj, Inf)
  
  ## remainder of starting run
  for(k in 2:start) {
    out = blackbox(X[k,], ...); obj[k] = out$obj; C[k,] = out$c
    ## update best so far
    prog[k] = ifelse(obj[k] < prog[k-1] && all(C[k,] <= 0), obj[k], prog[k-1]) 
  }
  
  ## best valid so far
  m2 = prog[start]
  if(is.finite(m2)){ # best solution so far
    xbest = X[which.min(prog),]
  }else{xbest = B[,2]} 
  
  ## initialize objective surrogate
  fnorm = max(abs(obj))
  fgpi = newGPsep(X, obj/fnorm, d=dg.start[1], g=dg.start[2], dK=TRUE)
  df = mleGPsep(fgpi, param = "d", tmin = dlim[1], tmax = dlim[2], ab = ab, 
                verb=verb-1)$d
  dfs = matrix(df, nrow = 1)
  
  ## initializing constraint surrogates
  Cgpi = rep(NA, nc)
  dc = matrix(NA, nrow=nc, ncol=nrow(B))
  Cnorm = apply(abs(C), 2, max)
  for (j in 1:nc) {
    Cgpi[j] = newGPsep(X, C[,j]/Cnorm[j], d=dg.start[1], g=dg.start[2], dK=TRUE)
    dc[j,] = mleGPsep(Cgpi[j], param = "d", tmin = dlim[1], tmax = dlim[2], ab = ab, 
                      verb=verb-1)$d
  }
  dcs <- matrix(rowMeans(dc, na.rm=TRUE), nrow=1)
  
  ## keeping track
  maes = c()
  
  ## iterating over the black box evaluations
  since = 0
  for(k in (start+1):end) {
    ## rebuild surrogates periodically under new normalized responses
    if(k > (start+1) && (k %% urate == 0)) {
      ##objective surrogate
      deleteGPsep(fgpi)
      df[df<dlim[1]] = 10*dlim[1]
      df[df>dlim[2]] = dlim[2]/10
      fnorm = max(abs(obj))
      fgpi = newGPsep(X, obj/fnorm, d=df, g=dg.start[2], dK=TRUE)
      df = mleGPsep(fgpi, param = "d", tmin = dlim[1], tmax = dlim[2], ab = ab,
                    verb=verb-1)$d
      dfs = rbind(dfs, df)
      
      ## constraint surrogates 
      Cnorm = apply(abs(C), 2, max)
      for(j in 1:nc) {
        deleteGPsep(Cgpi[j])
        dc[j, dc[j,]<dlim[1]] = 10*dlim[1]
        dc[j, dc[j,]>dlim[2]] = dlim[2]/10
        Cgpi[j] = newGPsep(X, C[,j]/Cnorm[j], d=dc[j,], g=dg.start[2], dK=TRUE)
        dc[j,] = mleGPsep(Cgpi[j], param = "d",  tmin = dlim[1], tmax = dlim[2], 
                          ab = ab, verb=verb-1)$d
      }
      dcs = rbind(dcs, rowMeans(dc, na.rm=TRUE))
    }
    
    ## random candidate grid
    ncand = ncandf(k)
    cands = lhs(ncand, B)
    
    ## evaluate acquisition function
    AF = AF_AE(cands, fgpi, fnorm, Cgpi, Cnorm, m2, alpha1, alpha2, omega) 
    m = which.max(AF)
    maes = c(maes, AF[m])
    
    ## calculate next point
    xnext = matrix(cands[m,], nrow = 1)
    X = rbind(X, xnext)
    
    if(verb > 0) {
      cat("k=", k, " ", sep="")
      cat("AsyEn=", AF[m], sep="")
      cat("; xnext=[", paste(signif(xnext,3), collapse=" "), sep="")
      cat("]; xbest=[", paste(signif(xbest,3), collapse=" "), sep="")
      cat("]; ybest=", m2, ", since=", since, "\n", sep="")
    }
    
    ## new run
    out = blackbox(xnext, ...)
    fnext = out$obj; obj = c(obj, fnext); C = rbind(C, out$c)
    
    ## check if best valid has changed
    since = since + 1
    valid <- apply(C, 1, function(x) { all(x <= 0) })
    if(all(out$c <= 0) && fnext < tail(prog,1)) {
      m2 = fnext; since = 0 
    } ## otherwise m2 unchanged; should be the same as prog[length(prog)]
    prog = c(prog, m2)
    xbest = X[which.min(prog),] # best solution so far
    
    ## update GP fits
    updateGPsep(fgpi, xnext, fnext/fnorm)
    for(j in 1:nc) updateGPsep(Cgpi[j], xnext, out$c[j]/Cnorm[j])
  }
  
  ## delete GP surrogates
  deleteGPsep(fgpi)
  for(j in 1:nc) deleteGPsep(Cgpi[j])
  
  return(list(prog = prog, xbest = xbest, obj = obj, C=C, X = X))
}
