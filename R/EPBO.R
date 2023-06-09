#' @title Exact Penalty Bayesian Optimization
#' 
#' @description Black-box optimization under mixed equality and inequality constraints via an exact penalty.
#' 
#' @param blackbox blackbox of an input (\code{x}), facilitating vectorization on a 
#' \code{matrix} \code{X} thereof,  returning a \code{list} 
#' with elements \code{$obj} containing the (scalar) objective value and \code{$c} 
#' containing a vector of evaluations of the (multiple) constraint function at \code{x}.
#' @param B 2-column \code{matrix} describing the bounding box.  The number of rows
#' of the \code{matrix} determines the input dimension (\code{length(x)} in \code{blackbox(x)}); 
#' the first column gives lower bounds and the second gives upper bounds
#' @param equal an optional vector containing zeros and ones, whose length equals the number of
#' constraints, specifying which should be treated as equality constraints (\code{1}) and 
#' which as inequality (\code{0}) 
#' @param ethresh a threshold used for equality constraints to determine validity for 
#' progress measures; ignored if there are no equality constraints
#' @param Xstart optional matrix of starting design locations in lieu of, or in addition to,
#' \code{start} random ones;  we recommend \code{nrow(Xstart) + start >= 6}; also must
#' have \code{ncol(Xstart) = nrow(B)}
#' @param start  positive integer giving the number of random starting locations before 
#' sequential design (for optimization) is performed; \code{start >= 6} is
#' recommended unless \code{Xstart} is non-\code{NULL}; in the current version
#' the starting locations come from a space-filling design via \code{\link[tgp]{dopt.gp}}
#' @param end positive integer giving the total number of evaluations/trials in the 
#' optimization; must have \code{end > start}
#' @param urate positive integer indicating  how many optimization trials should pass before
#' each MLE/MAP update is performed for GP correlation lengthscale 
#' parameter(s) 
#' @param rho positive vector initial exact penalty parameter in the exact penalty function; 
#' the default setting of \code{rho = NULL} causes an automatic starting value to be chosen
#' @param ncandf function taking a single integer indicating the optimization trial number \code{t}, where
#' \code{start < t <= end}, and returning the number of search candidates (e.g., for
#' expected improvement calculations) at round \code{t}; the default setting
#' allows the number of candidates to grow linearly with \code{t}
#' @param ey.tol  a scalar proportion indicating how many of the EIs
#' at \code{ncandf(t)} candidate locations must be non-zero to \dQuote{trust}
#' that metric to guide search, reducing to an EY-based search instead 
#' (choosing that proportion to be one forces EY-based search)
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
#'  \item{prog }{ vector giving the best feasible (\code{g(x) <= 0 && |h(x)| <= ethresh}) value of the objective over the trials }
#'  \item{xbest }{ vector giving the recommended solution}
#'  \item{obj }{ vector giving the value of the objective for the input under consideration at each trial }
#'  \item{C }{ \code{matrix} giving the value of the constraint function for the input under consideration at each trial}
#'  \item{X }{ \code{matrix} giving the input values at which the blackbox function was evaluated }
#'  
#'  @seealso \code{\link[laGP]{optim.auglag}}, \code{\link[laGP]{optim.efi}}, \code{\link[EPBO]{optim.AE}}, \code{\link[EPBO]{optim.BM}}
#'  
#' @author Jiangyan Zhao \email{zhaojy2017@126.com}
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
#' ## Inequality constrained problem
#' # search space
#' B = rbind(c(0, 1), c(0, 1)) 
#' # Objective and constraint function for use 
#' HSQ = function(x, B=rbind(c(0, 1), c(0, 1))){
#'   x = pmin(pmax(B[,1], x), B[,2])
#'   t = function(z)
#'     return(exp(-(z-1)^2) + exp(-0.8*(z+1)^2) - 0.05*sin(8*(z+0.1)))
#'   herbtooth = -t(4*x[1]-2)*t(4*x[2]-2)
#'   c1 = 1.5 - x[1] - 2*x[2] - 0.5*sin(2*pi*(x[1]^2 - 2*x[2]))
#'   c2 = x[1]^2 + x[2]^2 - 1.5
#'   return(list(obj=herbtooth, c=cbind(c1,c2)))
#' }
#' EPBO = optim.EP(HSQ, B, ncandf = function(k){ 1e3 }, start = 20, end = 120, verb = 0)
#' # progress, best feasible value of the objective over the trials
#' EPBO$prog
#' # the recommended solution
#' EPBO$xbest
#' 
#' ## Mixed constrained problem
#' #' # search space
#' B = rbind(c(0, 1), c(0, 1))
#' # Objective and constraint function for use 
#' goldstein.price = function(X) 
#' {
#'   if(is.null(nrow(X))) X <- matrix(X, nrow=1)
#'   m <- 8.6928
#'   s <- 2.4269
#'   x1 <- 4 * X[,1] - 2
#'   x2 <- 4 * X[,2] - 2
#'   a <- 1 + (x1 + x2 + 1)^2 * (19 - 14 * x1 + 3 * x1^2 - 14 * x2 + 6 * x1 * x2 + 3 * x2^2)
#'   b <- 30 + (2 * x1 - 3 * x2)^2 * (18 - 32 * x1 + 12 * x1^2 + 48 * x2 - 36 * x1 * x2 + 27 * x2^2)
#'   f <- log(a * b)
#'   f <- (f - m)/s
#'   return(f)
#' }
#' 
#' toy.c1 <- function(X) {
#'   if (is.null(dim(X))) X <- matrix(X, nrow=1)
#'   c1 = 3/2 - X[,1] - 2*X[,2] - .5*sin(2*pi*(X[,1]^2 - 2*X[,2]))
#' }
#' 
#' branin.c = function(X){
#'   if (is.null(dim(X))) X <- matrix(X, nrow=1)
#'   x1 = X[,1] * 15 - 5
#'   x2 = X[,2] * 15
#'   f = (x2 - 5/(4 * pi^2) * (x1^2) + 5/pi * x1 - 6)^2 + 10 * (1 - 1/(8 * pi)) * cos(x1) + 10
#'   return(-f+25)
#' }
#' 
#' parr <- function(X){
#'   if (is.null(dim(X))) X <- matrix(X, nrow=1)
#'   x1 <- (2 * X[,1] - 1)
#'   x2 <- (2 * X[,2] - 1)
#'   g <- (4-2.1*x1^2+1/3*x1^4)*x1^2 + x1*x2 + (-4+4*x2^2)*x2^2 + 
#'     3*sin(6*(1-x1)) + 3*sin(6*(1-x2))
#'   return(-g+4)
#' }
#' 
#' gsbp.constraints <- function(x){
#'   return(cbind(toy.c1(x), branin.c(x), parr(x)))
#' }
#' 
#' # problem definition
#' gsbpprob <- function(X, known.only=FALSE)
#' { 
#'   if(is.null(nrow(X))) X <- matrix(X, nrow=1)
#'   if(known.only) stop("known.only not supported for this example")
#'   f <- goldstein.price(X)
#'   C <- gsbp.constraints(X)
#'   return(list(obj=f, c=cbind(C[,1], C[,2]/100, C[,3]/10)))
#' }
#' EPBO = optim.EP(gsbpprob, B, ncandf = function(k){ 1e3 }, start = 20, end = 120, 
#'                 ethresh = 1e-2, equal = c(0,1,1), verb = 0)
#'                 
#' # progress, best feasible value of the objective over the trials
#' EPBO$prog
#' # the recommended solution
#' EPBO$xbest

optim.EP = function(blackbox, B, 
                    equal=FALSE, ethresh=1e-2, 
                    Xstart=NULL, start=10, end=100, 
                    urate=10, rho=NULL, 
                    ncandf=function(k) { k }, 
                    ey.tol=0.01,
                    dg.start=c(0.1, 1e-6), ab=c(3/2, 8), 
                    dlim=sqrt(ncol(B))*c(1/100, 10), verb=2, ...)
{
  ## check start
  if(start >= end) stop("must have start < end")
  ## get initial design
  X = dopt.gp(start, Xcand=lhs(10*start, B))$XX
  X = rbind(Xstart, X)
  start = nrow(X)
  
  ## first run to determine the number of constraints
  out = blackbox(X[1,], ...)
  nc = length(out$c)
  
  ## check equal argument
  if(length(equal) == 1) {
    if(!(equal %in% c(0,1))) 
      stop("equal should be a logical scalar or vector of lenghth(blackbox(x)$c)")
    if(equal) equal = rep(1, nc)
    else equal = rep(0, nc)
  } else {
    if(length(equal) != nc)
      stop("equal should be a logical scalar or vector of lenghth(blackbox(x)$c)")
    if(! any(equal != 0 | equal != 1))
      stop("equal should be a logical scalar or vector of lenghth(blackbox(x)$c)")
  }
  equal = as.logical(equal)
  
  ## allocate progress objects, and initialize
  prog = obj = rep(NA, start)
  C = CV = matrix(NA, nrow=start, ncol=nc) # constraints and constraint violations
  obj[1] = out$obj; C[1, ] = out$c
  CV[1, !equal] = pmax(0, C[1, !equal]) # constraints violation for inequality
  CV[1, equal] = abs(C[1, equal]) # constraints violation for equality
  
  prog[1] = ifelse(all(C[1,!equal] <= 0) && all(abs(C[1,equal]) <= ethresh), obj[1], Inf)
  
  ## remainder of starting run
  for(k in 2:start) { ## now that problem is vectorized we can probably remove for
    out = blackbox(X[k,], ...)
    obj[k] = out$obj; C[k,] = out$c
    CV[k, !equal] = pmax(0, C[k, !equal]) # constraints violation for inequality
    CV[k, equal] = abs(C[k, equal]) # constraints violation for equality
    ## update best so far
    prog[k] = ifelse(obj[k] < prog[k-1] && all(C[k,!equal] <= 0) && all(abs(C[k,equal]) <= ethresh), obj[k], prog[k-1]) 
  }
  
  ## handle initial rho value
  idx_v = c() # the index of validity points
  for (k in 1:start) {
    if(all(C[k,!equal] <= 0) && all(abs(C[k,equal]) <= ethresh)) idx_v = c(idx_v, k)
  }
  
  ## handle initial rho values
  if(is.null(rho)){
    ECV = colMeans(CV) # the constraint violation averaged over the current design
    if(all(ECV == 0)){ # 
      rho = rep(0, nc)
    }else {
      rho = mean(abs(obj)) * ECV/sum(ECV^2)
    }
    if(any(equal)) rho[equal] = pmax(1/ethresh/sum(equal), rho[equal])
  }else{
    if(length(rho) != nc || any(rho < 0)) stop("rho should be a non-negative vector whose length is equal to the number of constraints.")
  }
  ## printing initial
  if(verb > 0){ cat("initial rho=(", paste(signif(rho,3), collapse=", "), ")\n", sep="") }
  
  ## calculate EP for data seen so far
  ep = obj + CV%*%rho
  ybest = min(ep) # best EP seen so far
  since = 0
  m2 = prog[start] # best valid solution so far
  if(is.finite(m2)){ xbest = X[which.min(prog),] }else{ xbest = B[,2] } # best solution so far
  
  ## initialize objective surrogate
  fmean = mean(obj); fsd = sd(obj)
  fgpi = newGPsep(X, (obj-fmean)/fsd, d = dg.start[1], g = dg.start[2], dK = TRUE)
  df = mleGPsep(fgpi, param = "d", tmin = dlim[1], tmax = dlim[2], ab = ab, verb = verb-1)$d
  dfs = matrix(df, nrow = 1)
  
  ## initializing constraint surrogates
  Cgpi = rep(NA, nc)
  dc = matrix(NA, nrow=nc, ncol=nrow(B))
  Cmean = apply(abs(C), 2, mean); Csd = apply(abs(C), 2, sd)
  for (j in 1:nc) {
    Cgpi[j] = newGPsep(X, (C[,j]-Cmean[j])/Csd[j], d=dg.start[1], g=dg.start[2], dK=TRUE)
    dc[j,] = mleGPsep(Cgpi[j], param = "d", tmin = dlim[1], tmax = dlim[2], ab = ab, verb=verb-1)$d
  }
  dcs = matrix(rowMeans(dc, na.rm=TRUE), nrow=1)
  
  ## keeping track
  mseis = c()
  rhos = rho
  
  ## iterating over the black box evaluations
  for (k in (start+1):end) {
    ## rebuild surrogates periodically under new normalized responses
    if(k > (start+1) && (k %% urate == 0)) {
      ## objective surrogate
      deleteGPsep(fgpi)
      df[df<dlim[1]] = 10*dlim[1]
      df[df>dlim[2]] = dlim[2]/10
      fmean = mean(obj); fsd = sd(obj)
      fgpi = newGPsep(X, (obj-fmean)/fsd, d=df, g=dg.start[2], dK=TRUE)
      df = mleGPsep(fgpi, param = "d", tmin = dlim[1], tmax = dlim[2], ab = ab, verb=verb-1)$d
      dfs = rbind(dfs, df)
      
      ## constraint surrogates 
      Cmean = apply(abs(C), 2, mean); Csd = apply(abs(C), 2, sd)
      for(j in 1:nc) {
        deleteGPsep(Cgpi[j])
        dc[j, dc[j,]<dlim[1]] = 10*dlim[1]
        dc[j, dc[j,]>dlim[2]] = dlim[2]/10
        Cgpi[j] = newGPsep(X, (C[,j]-Cmean[j])/Csd[j], d=dc[j,], g=dg.start[2], dK=TRUE)
        dc[j,] = mleGPsep(Cgpi[j], param = "d",  tmin = dlim[1], tmax = dlim[2], ab = ab, verb=verb-1)$d
      }
      dcs = rbind(dcs, rowMeans(dc, na.rm=TRUE))
    }
    
    
    ## Update the rho when the smallest EP is not feasible
    if(is.finite(m2)){
      ck = C[which.min(ep),]
      invalid = rep(NA, nc)
      invalid[!equal] = (ck[!equal] > 0); invalid[equal] = (abs(ck[equal]) > ethresh)
      nupdate = 0
      while(any(invalid) && nupdate<=5){
        rho[invalid] = rho[invalid]*2
        ep = obj + CV%*%rho
        ck = C[which.min(ep),]
        invalid[!equal] = (ck[!equal] > 0); invalid[equal] = (abs(ck[equal]) > ethresh)
        nupdate = nupdate + 1
      }
      ybest = min(ep)
    }
    
    ## random candidate grid
    ncand = ncandf(k)
    cands = lhs(ncand, B)
    
    ## calculate composite surrogate, and evaluate EI and/or EY
    AF = EP_AcqFunc(cands, fgpi, fmean, fsd, Cgpi, Cmean, Csd, ybest, rho, 
                    equal, blackbox, eiey="ei")
    nzsei = sum(AF > sqrt(.Machine$double.eps))
    if(nzsei <= ey.tol*ncand){ # minimization predictive mean approach
      AF = EP_AcqFunc(cands, fgpi, fmean, fsd, Cgpi, Cmean, Csd, ybest, rho, 
                      equal, blackbox, eiey="ey")
      m = which.min(AF)
      out = optim(par=cands[m, ], fn=EP_AcqFunc, method="L-BFGS-B",
                  lower=B[,1], upper=B[,2],
                  fgpi=fgpi, fmean=fmean, fsd=fsd,
                  Cgpi=Cgpi, Cmean=Cmean, Csd=Csd,
                  eiey="ey",ybest=ybest, rho=rho, equal=equal, 
                  blackbox=blackbox)
      mseis = c(mseis, Inf)
      by = "ey"
    }else if(nzsei <= 0.1*ncand){ # maximization expected improvement approach
      cands = rbind(cands, lhs(10*ncand, B))
      AF = c(AF, EP_AcqFunc(cands[-(1:ncand),], fgpi, fmean, fsd, Cgpi, Cmean, Csd, ybest, rho, 
                            equal, blackbox, eiey="ei"))
      nzsei = sum(AF > 0)
      if(nzsei <= ey.tol*ncand) stop("not sure")
      ncand = ncand + 10*ncand
      
      m = which.max(AF)
      out = optim(par=cands[m, ], fn=EP_AcqFunc, method="L-BFGS-B",
                  lower=B[,1], upper=B[,2], 
                  control = list(fnscale = -1), # maximization problem
                  fgpi=fgpi, fmean=fmean, fsd=fsd,
                  Cgpi=Cgpi, Cmean=Cmean, Csd=Csd,
                  eiey="ei", ybest=ybest, rho=rho, equal=equal, 
                  blackbox=blackbox)
      mseis = c(mseis, out$value)
      by = "ei"
    }else{ # maximization expected improvement approach
      m = which.max(AF)
      out = optim(par=cands[m, ], fn=EP_AcqFunc, method="L-BFGS-B",
                  lower=B[,1], upper=B[,2], 
                  control = list(fnscale = -1), # maximization problem
                  fgpi=fgpi, fmean=fmean, fsd=fsd,
                  Cgpi=Cgpi, Cmean=Cmean, Csd=Csd,
                  eiey="ei", ybest=ybest, rho=rho, equal=equal, 
                  blackbox=blackbox)
      mseis = c(mseis, out$value)
      by = "ei"
    }
    
    
    ## calculate next point
    xnext = matrix(out$par, nrow = 1)
    X = rbind(X, xnext)
    
    ## new run
    out = blackbox(xnext, ...)
    fnext = out$obj; obj = c(obj, fnext); C = rbind(C, out$c)
    CV = rbind(CV, rep(NA, nc))
    CV[k, !equal] = pmax(0, C[k, !equal]) # constraints violation for inequality
    CV[k, equal] = abs(C[k, equal]) # constraints violation for equality
    
    ## the index of validity points
    if(all(out$c[!equal] <= 0) && all(abs(out$c[equal]) <= ethresh)) idx_v = c(idx_v, k)
    
    ## check if best valid has changed
    since = since + 1
    if(all(out$c[!equal] <= 0) && all(abs(out$c[equal]) <= ethresh) && fnext < tail(prog,1)) {
      m2 = fnext; since = 0
    } # otherwise m2 unchanged; should be the same as tail(prog,1)
    prog = c(prog, m2)
    if(is.finite(m2)) xbest = X[which.min(prog),] # best solution so far
    
    ## rho update
    ECV = colMeans(CV) # the constraint violation averaged over the current population
    if(all(ECV == 0)){ # 
      rho_new = rep(0, nc)
    }else {
      rho_new = mean(abs(obj)) * ECV/sum(ECV^2)
    }
    rho = pmax(rho_new, rho)
    rhos = rbind(rhos, rho) # keep track of rho
    
    ## calculate EP for data seen so far
    ep = obj + CV%*%rho
    ybest = min(ep)
    
    ## printing progress
    if(verb > 0 && any(rho_new > rho)) cat("updating rho=(", paste(signif(rho,3), collapse=", "), ")\n", sep="")
    
    ## progress meter
    if(verb > 0) {
      cat("k=", k, " ", sep="")
      cat(by, "=", AF[m], " (", nzsei,  "/", ncand, ")", sep="")
      cat("; xnext=[", paste(signif(xnext,3), collapse=" "), sep="")
      cat("]; xbest=[", paste(signif(xbest,3), collapse=" "), sep="")
      cat("]; ybest (prog=", m2, ", ep=", ybest, ", since=", since, ")\n", sep="")
    }
    
    
    ## update GP fits
    updateGPsep(fgpi, xnext, (fnext-fmean)/fsd, verb = verb-2)
    df = mleGPsep(fgpi, param = "d", tmin = dlim[1], tmax = dlim[2], ab = ab, verb=verb-1)$d
    dfs = rbind(dfs, df)
    
    for(j in 1:nc) {
      updateGPsep(Cgpi[j], xnext, (out$c[j]-Cmean[j])/Csd[j], verb = verb-2)
      dc[j,] = mleGPsep(Cgpi[j], param = "d",  tmin = dlim[1], tmax = dlim[2], ab = ab, verb=verb-1)$d
    }
    dcs = rbind(dcs, rowMeans(dc, na.rm=TRUE))
  }
  
  ## delete GP surrogates
  deleteGPsep(fgpi)
  for(j in 1:nc) deleteGPsep(Cgpi[j])
  
  return(list(prog = prog, xbest = xbest, obj = obj, C=C, X = X))
}


EP_AcqFunc = function(x, fgpi, fmean, fsd, Cgpi, Cmean, Csd, ybest, rho, equal, blackbox=NULL, eiey="ei")
{
  if(is.null(nrow(x))) x = matrix(x, nrow=1)
  ncand = nrow(x)
  nc = length(Cgpi) # number of the constraint
  
  ## Acquaisition function
  if(eiey == "ei"){
    ## objective
    pred_f = predGPsep(fgpi, x, lite=TRUE)
    mu_f = pred_f$mean * fsd + fmean
    sigma_f =  sqrt(pred_f$s2) * fsd
    
    ## constraints
    mu_C = sigma_C = omega = matrix(NA, nc, ncand)
    
    for (j in 1:nc) {
      pred_C = predGPsep(Cgpi[j], x, lite=TRUE)
      mu_C[j,] = pred_C$mean * Csd[j] + Cmean[j]
      sigma_C[j,] = sqrt(pred_C$s2) * Csd[j]
      omega[j,] = (equal[j]+1)*pnorm(mu_C[j,]/sigma_C[j,]) - equal[j]
    }
    
    ## Acquaisition function
    mu_ep = mu_f + rho%*%(omega*mu_C)
    sigma_ep = sqrt(sigma_f^2 + (rho^2)%*%((omega*sigma_C)^2))
    d = (ybest - mu_ep)/sigma_ep
    # # expected improvement
    # AF = sigma_ep*(d*pnorm(d) + dnorm(d))
    # AF[is.nan(AF)] = 0 # AF=NaN if sigma_ep = 0
    # Scaled expected improvement
    numerator = d*pnorm(d) + dnorm(d) # numerator
    denominator = pnorm(d) + d*numerator - numerator^2 #denominator
    denominator = pmax(.Machine$double.xmin, denominator)
    AF = numerator/sqrt(denominator)
    AF[is.nan(AF)] = 0 # AF=NaN if sigma_ep = 0
  }else if(eiey == "ey"){
    ## objective
    pred_f = predGPsep(fgpi, x, lite=TRUE)
    mu_f = pred_f$mean * fsd + fmean
    
    ## constraints
    mu_C = sigma_C = EV = matrix(NA, nc, nrow(x))
    
    for (j in 1:nc) {
      pred_C = predGPsep(Cgpi[j], x, lite=TRUE)
      mu_C[j,] = pred_C$mean * Csd[j] + Cmean[j]
      sigma_C[j,] = sqrt(pred_C$s2) * Csd[j]
      dC = mu_C[j,]/sigma_C[j,]
      EV[j,] = mu_C[j,]*((equal[j]+1)*pnorm(dC) - equal[j]) + 
        (equal[j]+1)*sigma_C[j,]*dnorm(dC)
    }
    
    ## Acquaisition function
    AF = mu_f + rho%*%EV
  }else{
    stop("eiey should be 'ei' or 'ey'.")
  }
  return(AF)
}
