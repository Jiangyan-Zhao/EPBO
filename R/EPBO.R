#' @title EPBO: Bayesian Optimization via Exact Penalty
#' 
#' @description Black-box optimization under mixed equality and inequality constraints via an exact penalty.
#' 
#' @param blackbox blackbox of an input (\code{x}), facilitating vectorization on a 
#' \code{matrix} \code{X} thereof, returning a \code{list} 
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
#' @param rho positive vector initial exact penalty parameters in the exact penalty function; 
#' the default setting of \code{rho = NULL} causes the automatic starting values to be chosen
#' @param ncandf function taking a single integer indicating the optimization trial number \code{k}, where
#' \code{start < k <= end}, and returning the number of search candidates (e.g., for
#' expected improvement calculations) at round \code{k}; the default setting
#' allows the number of candidates to grow linearly with \code{k}
#' @param dg_start 2-vector giving starting values for the lengthscale and nugget parameters
#' of the GP surrogate model(s) for the objective and constraints
#' @param dlim 2-vector giving bounds for the lengthscale parameter(s) under MLE inference
#' @param plotprog \code{logical} indicating if progress plots should be made after each inner iteration;
#' the plots show four panels tracking the best feasible objective values, the ScaledEI or EY surface, 
#' the predictive mean and standard deviation of the objective function, over the first two input variables
#' @param ey_tol a scalar proportion indicating how many of the ScaledEIs 
#' at \code{ncandf(k)} candidate locations must be non-zero to \dQuote{trust}
#' that metric to guide search, reducing to an EY-based search instead 
#' (choosing that proportion to be one forces EY-based search)
#' @param verb a non-negative integer indicating the verbosity level; the larger the value the
#' more that is printed to the screen
#' @param ... additional arguments passed to \code{blackbox}
#' 
#' @details
#' To be added.
#' 
#'
#' @returns The output is a \code{list} summarizing the progress of the evaluations of the 
#' blackbox under optimization:
#'  \item{prog }{ \code{vector} giving the best feasible (\code{g(x) <= 0 && |h(x)| <= ethresh}) value of the objective over the trials }
#'  \item{xbest }{ \code{vector} giving the recommended solution}
#'  \item{obj }{ \code{vector} giving the value of the objective for the input under consideration at each trial }
#'  \item{C }{ \code{matrix} giving the value of the constraint function for the input under consideration at each trial}
#'  \item{X }{ \code{matrix} giving the input values at which the blackbox function was evaluated }
#'  \item{feasibility}{vector giving the feasibility for the input under consideration at each trial }
#'  \item{rho}{vector giving the penalty parameters }
#'  
#' @seealso \code{\link[laGP]{optim.auglag}}, \code{\link[laGP]{optim.efi}}, \code{\link[EPBO]{optim.AE}}, \code{\link[EPBO]{optim.BM}}
#'  
#' @author Jiangyan Zhao \email{zhaojy2017@126.com}
#' 
#' @references Jiangyan Zhao and Jin Xu (07 Feb 2024): Bayesian Optimization via Exact Penalty, 
#' \emph{Technometrics}, DOI: 10.1080/00401706.2024.2315937
#' 
#' @keywords optimize
#' @keywords design
#' 
#' @import laGP
#' @import tgp
#' @importFrom graphics image
#' @importFrom graphics par
#' @importFrom graphics points
#' @importFrom stats dnorm
#' @importFrom stats median
#' @importFrom stats optim 
#' @importFrom stats pnorm 
#' @importFrom stats sd
#' @importFrom stats runif
#' @importFrom utils tail
#' 
#' @export
#'
#' @examples 
#' ## Inequality constrained problem
#' # search space
#' B = rbind(c(0, 1), c(0, 1)) 
#' # Objective and constraint function for use 
#' HSQ = function(x){
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
#' # search space
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
#'   f = goldstein.price(X)
#'   C = gsbp.constraints(X)
#'   return(list(obj=f, c=cbind(C[,1], C[,2]/100, C[,3]/10)))
#' }
#' EPBO = optim.EP(gsbpprob, B, ncandf = function(k){ 1e3 }, start = 20, end = 120, 
#'                 ethresh = 1e-2, equal = c(0,1,1), verb = 0)
#'                 
#' # progress, best feasible value of the objective over the trials
#' EPBO$prog
#' # the recommended solution
#' EPBO$xbest

optim.EP = function(
    blackbox, B, equal=FALSE, ethresh=1e-2, 
    Xstart=NULL, start=10, end=100, urate=10, rho=NULL, 
    ncandf=function(k) { k }, 
    dg_start=c(1e-2*sqrt(nrow(B)), 1e-6), 
    dlim=c(1e-4, 1)*sqrt(nrow(B)), 
    plotprog=FALSE, ey_tol=1e-2,
    verb=2, ...)
{
  ## check start
  if(start >= end) stop("must have start < end")
  if(start == 0 & is.null(Xstart)) stop("must have start>0 or given Xstart")
  
  dim = nrow(B) # dimension
  # Hypercube [0, 1]^d
  Hypercube = matrix(c(rep(0,dim), rep(1,dim)), ncol=2) 
  
  ## get initial design
  X = dopt.gp(start, Xcand=lhs(10*start, B))$XX
  ## X = lhs(start, B)
  X = rbind(Xstart, X)
  X_unit = normalize(X, B)
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
  prog = rep(NA, start)                       # progress
  obj = rep(NA, start)                        # objective values
  feasibility = rep(NA, start)                # Whether the point is feasible?
  C = matrix(NA, nrow=start, ncol=nc)         # constraint values
  C_bilog = matrix(NA, nrow=start, ncol=nc)   # bi-log transformation of constraint values
  CV = matrix(NA, nrow=start, ncol=nc)        # constraint violations in terms of bi-log transformation
  
  obj[1] = out$obj; C[1, ] = out$c; C_bilog[1, ] = bilog(C[1, ])
  CV[1, !equal] = pmax(0, C_bilog[1, !equal]) # CV for inequality
  CV[1, equal] = abs(C_bilog[1, equal])       # CV for equality
  feasibility[1] = all(C[1,!equal] <= 0) && all(abs(C[1,equal]) <= ethresh)
  prog[1] = ifelse(feasibility[1], obj[1], Inf)
  
  ## remainder of starting run
  for(k in 2:start) { ## now that problem is vectorized we can probably remove for
    out = blackbox(X[k,], ...)
    obj[k] = out$obj; C[k,] = out$c; C_bilog[k, ] = bilog(C[k, ])
    CV[k, !equal] = pmax(0, C_bilog[k, !equal])
    CV[k, equal] = abs(C_bilog[k, equal])
    feasibility[k] = all(C[k,!equal] <= 0) && all(abs(C[k,equal]) <= ethresh)
    prog[k] = ifelse(feasibility[k] && obj[k] < prog[k-1], obj[k], prog[k-1]) 
  }
  
  ## handle initial rho values
  if(is.null(rho)){
    if(all(feasibility)){            # 
      rho = rep(0, nc)
    }else {
      ECV = colMeans(CV) # averaged CV
      rho = mean(abs(obj)) * ECV/sum(ECV^2)
    }
    if(any(equal)) rho[equal] = pmax(1/ethresh/sum(equal), rho[equal])
  }else{
    if(length(rho) != nc || any(rho < 0)) stop("rho should be a non-negative vector whose length is equal to the number of constraints.")
  }
  
  ## calculate EP for data seen so far
  scv = CV%*%rho        # weighted sum of the constraint violations
  ep = obj + scv        # the EP values
  epbest = min(ep)      # best EP seen so far
  m2 = prog[start]      # BOFV
  since = 0
  ## best solution so far
  if(is.finite(m2)){            # if at least one feasible solution was found
    xbest = X[which.min(prog),] 
  }else{                        # if none of the solutions are feasible
    xbest = X[which.min(scv),] 
  }
  xbest_unit = as.vector(normalize(xbest, B))
  
  ## initialize objective surrogate
  ab = darg(NULL, X_unit)$ab
  fmean = mean(obj); fsd = sd(obj) # for standard normalization on objective values
  fgpi = newGPsep(X_unit, (obj-fmean)/fsd, d = dg_start[1], g = dg_start[2], dK = TRUE)
  df = mleGPsep(fgpi, param = "d", tmin = dlim[1], tmax = dlim[2], ab = ab, verb = verb-1)$d
  deleteGPsep(fgpi)
  df[df<dlim[1]] = 10*dlim[1]
  df[df>dlim[2]] = dlim[2]/10
  fgpi = newGPsep(X_unit, (obj-fmean)/fsd, d=df, g=dg_start[2], dK=TRUE)
  df = mleGPsep(fgpi, param = "d", tmin = dlim[1], tmax = dlim[2], ab = ab, verb=verb-1)$d
  
  ## initializing constraint surrogates
  Cgpi = rep(NA, nc)
  dc = matrix(NA, nrow=nc, ncol=dim)
  for (j in 1:nc) {
    Cgpi[j] = newGPsep(X_unit, C_bilog[,j], d=dg_start[1], g=dg_start[2], dK=TRUE)
    dc[j,] = mleGPsep(Cgpi[j], param = "d", tmin = dlim[1], tmax = dlim[2], ab = ab, verb=verb-1)$d
    deleteGPsep(Cgpi[j])
    dc[j, dc[j,]<dlim[1]] = 10*dlim[1]
    dc[j, dc[j,]>dlim[2]] = dlim[2]/10
    Cgpi[j] = newGPsep(X_unit, C_bilog[,j], d=dc[j,], g=dg_start[2], dK=TRUE)
    dc[j,] = mleGPsep(Cgpi[j], param = "d",  tmin = dlim[1], tmax = dlim[2], ab = ab, verb=verb-1)$d
  }
  
  ## printing initial design
  if(verb > 0) {
    cat("The initial design: ")
    cat("ab=[", paste(signif(ab,3), collapse=", "), sep="")
    cat("]; rho=[", paste(signif(rho,3), collapse=", "), sep="")
    cat("]; xbest=[", paste(signif(xbest,3), collapse=" "), sep="")
    cat("]; ybest (prog=", m2, ", ep=", epbest, ", since=", since, ")\n", sep="")
  }
  
  # AF_time = 0 # AF running time
  
  ## iterating over the black box evaluations
  for (k in (start+1):end) {
    ## rebuild surrogates periodically under new normalized responses
    if(k > (start+1) && (k %% urate == 0)) {
      ## objective surrogate
      deleteGPsep(fgpi)
      df[df<dlim[1]] = 10*dlim[1]
      df[df>dlim[2]] = dlim[2]/10
      fmean = mean(obj); fsd = sd(obj)
      fgpi = newGPsep(X_unit, (obj-fmean)/fsd, d=df, g=dg_start[2], dK=TRUE)
      df = mleGPsep(fgpi, param = "d", tmin = dlim[1], tmax = dlim[2], ab = ab, verb=verb-1)$d
      
      ## constraint surrogates 
      for(j in 1:nc) {
        deleteGPsep(Cgpi[j])
        dc[j, dc[j,]<dlim[1]] = 10*dlim[1]
        dc[j, dc[j,]>dlim[2]] = dlim[2]/10
        Cgpi[j] = newGPsep(X_unit, C_bilog[,j], d=dc[j,], g=dg_start[2], dK=TRUE)
        dc[j,] = mleGPsep(Cgpi[j], param = "d",  tmin = dlim[1], tmax = dlim[2], ab = ab, verb=verb-1)$d
      }
    }
    
    
    ## Update the rho when the smallest EP is not feasible
    ck = C[which.min(ep),]
    invalid = rep(NA, nc)
    invalid[!equal] = (ck[!equal] > 0); invalid[equal] = (abs(ck[equal]) > ethresh)
    if(is.finite(m2) && any(invalid) && since > 2){
      rho[invalid] = rho[invalid] * 2
      scv = CV%*%rho; ep = obj + scv; epbest = min(ep)
      if(verb > 0) cat(" the smallest EP is not feasible, updating rho=(", 
                       paste(signif(rho,3), collapse=", "), ")\n", sep="")
    }
    
    ## calculate composite surrogate, and evaluate SEI and/or EY
    ncand = ncandf(k)
    cands = lhs(ncand, Hypercube) # random candidate grid
    # tic = proc.time()[3] # Start time
    AF = AF_ScaledEI(cands, fgpi, fmean, fsd, Cgpi, epbest, rho, equal)
    nzsei = sum(AF > sqrt(.Machine$double.eps))
    # Augment the candidate points
    if((ey_tol*ncand < nzsei && nzsei <= 0.1*ncand) || (since>1 && since %% 5 == 0)){
      cands = rbind(cands, lhs(10*ncand, Hypercube))
      AF = c(AF, AF_ScaledEI(cands[-(1:ncand),], fgpi, fmean, fsd, Cgpi, epbest, rho, equal))
      nzsei = sum(AF > sqrt(.Machine$double.eps))
      ncand = 11*ncand
    }
    if(nzsei <= ey_tol*ncand){ # minimize predictive mean approach
      by = "EY"
      AF = AF_EY(cands, fgpi, fmean, fsd, Cgpi, rho, equal)
      m = which.min(AF)
      out_AF = optim(par=cands[m, ], fn=AF_EY, method="L-BFGS-B",
                     lower=0, upper=1,
                     fgpi=fgpi, Cgpi=Cgpi, fmean=fmean, fsd=fsd, 
                     rho=rho, equal=equal)
    }else{# maximize scaled expected improvement approach
      by = "ScaledEI"
      m = which.max(AF)
      out_AF = optim(par=cands[m, ], fn=AF_ScaledEI, method="L-BFGS-B",
                     lower=0, upper=1,
                     fgpi=fgpi, Cgpi=Cgpi, fmean=fmean, fsd=fsd, 
                     control = list(fnscale = -1), # maximization problem
                     epbest=epbest, rho=rho, equal=equal)
    }
    # toc = proc.time()[3] # End time
    # AF_time = AF_time + toc - tic # AF running time
    
    ## calculate next point
    xnext_unit = matrix(out_AF$par, nrow = 1)
    X_unit = rbind(X_unit, xnext_unit)
    xnext = unnormalize(xnext_unit, B)
    X = rbind(X, xnext)
    
    ## new run
    out = blackbox(xnext, ...)
    fnext = out$obj; obj = c(obj, fnext); C = rbind(C, out$c)
    C_bilog = rbind(C_bilog, bilog(out$c))
    CV = rbind(CV, rep(NA, nc))
    CV[k, !equal] = pmax(0, C_bilog[k, !equal]) # constraint violations for inequality
    CV[k, equal] = abs(C_bilog[k, equal]) # constraint violations for equality
    
    ## check if best valid has changed
    feasibility = c(feasibility, all(out$c[!equal] <= 0) && all(abs(out$c[equal]) <= ethresh))
    since = since + 1
    if(feasibility[k] && fnext < prog[k-1]) {
      m2 = fnext; since = 0
    } # otherwise m2 unchanged; should be the same as prog[k-1]
    prog = c(prog, m2)
    
    ## rho update
    if(all(feasibility)){ # 
      rho_new = rep(0, nc)
    }else {
      ECV = colMeans(CV)
      rho_new = mean(abs(obj)) * ECV/sum(ECV^2)
    }
    if(verb > 0 && any(rho_new > rho)){ # printing progress
      cat("  updating rho=[", paste(signif(pmax(rho_new, rho),3), collapse=", "), "]\n", sep="")
    }
    rho = pmax(rho_new, rho)
    
    ## calculate EP for data seen so far
    scv  = CV%*%rho; ep = obj + scv; epbest = min(ep)
    if(is.finite(m2)){ # best solution so far
      xbest = X[which.min(prog),] 
    }else{ 
      xbest = X[which.min(scv),] 
    }
    xbest_unit = as.vector(normalize(xbest, B))
    
    ## progress meter
    if(verb > 0) {
      cat("k=", k, " ", sep="")
      cat(by, "=", paste(signif(out_AF$value,3), collapse=" "), ", ncand=", ncand, sep="")
      cat("; xnext ([", paste(signif(xnext,3), collapse=" "), 
          "], feasibility=", feasibility[k], ")\n", sep="")
      cat(" xbest=[", paste(signif(xbest,3), collapse=" "), sep="")
      cat("]; ybest (prog=", paste(signif(m2,3), collapse=" ")) 
      cat(", ep=", paste(signif(epbest,3), collapse=" "))
      cat(", since=", since, ")\n", sep="")
    }
    
    ## update GP fits
    updateGPsep(fgpi, xnext_unit, (obj[k]-fmean)/fsd, verb = 0)
    df = mleGPsep(fgpi, param = "d", tmin = dlim[1], tmax = dlim[2], ab = ab, verb = 0)$d
    sigma_f = sqrt(predGPsep(fgpi, X_unit, lite=TRUE)$s2) #* fsd
    if(median(sigma_f) > 0.1){# rebuild GP as sigma is large at observed points
      if(verb > 0) cat(" rebuild objective surrogate as sigma is large at observed points \n")
      deleteGPsep(fgpi)
      fgpi = newGPsep(X_unit, (obj-fmean)/fsd, d = dlim[1], g = dg_start[2], dK=TRUE)
      df = mleGPsep(fgpi, param = "d", tmin = dlim[1], tmax = dlim[2], ab = ab, verb = 0)$d
    }
    
    for(j in 1:nc) {
      updateGPsep(Cgpi[j], xnext_unit, C_bilog[k,j], verb = 0)
      dc[j,] = mleGPsep(Cgpi[j], param = "d",  tmin = dlim[1], tmax = dlim[2], ab = ab, verb = 0)$d
      sigma_c = sqrt(predGPsep(Cgpi[j], X_unit, lite=TRUE)$s2)
      if(median(sigma_c) > 0.1){ # rebuild GP as sigma is large at observed points
        if(verb > 0) cat(" rebuild constrained", j, "surrogate as sigma is large at observed points \n")
        deleteGPsep(Cgpi[j])
        Cgpi[j] = newGPsep(X_unit, C_bilog[,j], d = dlim[1], g = dg_start[2], dK=TRUE)
        dc[j,] = mleGPsep(Cgpi[j], param = "d",  tmin = dlim[1], tmax = dlim[2], ab = ab, verb = 0)$d
      }
    }
    
    ## plot progress
    if(plotprog) {
      par(ps=16, mfrow=c(2,2))
      ## progress
      if(is.finite(m2)){
        plot(prog, type="l", lwd=1.6, 
             xlab="n", ylab="BFOV", main="progress")
      }else{
        plot(prog, type="l", ylim=range(obj), lwd=1.6, 
             xlab="n", ylab="BFOV", main="progress")
      }
      ## acquisition function
      cands_unnormalize = unnormalize(cands, B)
      span = ifelse(length(AF) < 30, 0.5, 0.1)
      graphic = interp.loess(cands_unnormalize[,1], cands_unnormalize[,2], 
                             as.vector(AF), span=span)
      image(graphic, xlim=range(X[,1]), ylim=range(X[,2]), main=by)
      points(X[1:start,1:2], col=feasibility[1:start]+3)
      points(X[-(1:start),1:2, drop=FALSE], col=feasibility[-(1:start)]+3, pch=19)
      points(X[k,,drop=FALSE], col="black", pch=18, cex=1.6)
      ## mean of objective function
      pred_f = predGPsep(fgpi, cands, lite=TRUE)
      mu_f = pred_f$mean * fsd + fmean
      sigma_f = sqrt(pred_f$s2) * fsd
      graphic = interp.loess(cands_unnormalize[,1], cands_unnormalize[,2],
                             mu_f, span=span)
      image(graphic, xlim=range(X[,1]), ylim=range(X[,2]), main="mu_f")
      points(X[1:start,1:2], col=feasibility[1:start]+3)
      points(X[-(1:start),1:2, drop=FALSE], col=feasibility[-(1:start)]+3, pch=19)
      points(X[k,,drop=FALSE], col="black", pch=18, cex=1.6)
      ## standard deviation of objective function
      graphic = interp.loess(cands_unnormalize[,1], cands_unnormalize[,2],
                             sigma_f, span=span)
      image(graphic, xlim=range(X[,1]), ylim=range(X[,2]), main="sd_f")
      points(X[1:start,1:2], col=feasibility[1:start]+3)
      points(X[-(1:start),1:2, drop=FALSE], col=feasibility[-(1:start)]+3, pch=19)
      points(X[k,,drop=FALSE], col="black", pch=18, cex=1.6)
    }
  }
  
  ## delete GP surrogates
  deleteGPsep(fgpi)
  for(j in 1:nc) deleteGPsep(Cgpi[j])
  
  return(list(prog = prog, xbest = xbest, 
              obj = obj, C=C, X = X, 
              feasibility=feasibility, rho=rho)) #AF_time = AF_time
}
