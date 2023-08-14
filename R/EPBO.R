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
#' @param dg_start 2-vector giving starting values for the lengthscale and nugget parameters
#' of the GP surrogate model(s) for constraints
#' @param dlim 2-vector giving bounds for the lengthscale parameter(s) under MLE inference
#' @param trcontrol an optional list of control parameters of the trust region. See "Details".
#' @param verb a non-negative integer indicating the verbosity level; the larger the value the
#' more that is printed to the screen
#' @param ... additional arguments passed to \code{blackbox}
#' 
#' @details
#' Additional details...
#' 
#'
#' @returns The output is a \code{list} summarizing the progress of the evaluations of the 
#' blackbox under optimization:
#'  \item{prog }{ vector giving the best feasible (\code{g(x) <= 0 && |h(x)| <= ethresh}) value of the objective over the trials }
#'  \item{xbest }{ vector giving the recommended solution}
#'  \item{obj }{ vector giving the value of the objective for the input under consideration at each trial }
#'  \item{C }{ \code{matrix} giving the value of the constraint function for the input under consideration at each trial}
#'  \item{X }{ \code{matrix} giving the input values at which the blackbox function was evaluated }
#'  \item{idx_v}{index of the valid solutions}
#'  
#' @seealso \code{\link[laGP]{optim.auglag}}, \code{\link[laGP]{optim.efi}}, \code{\link[EPBO]{optim.AE}}, \code{\link[EPBO]{optim.BM}}
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
#'   B = rbind(c(0, 1), c(0, 1))
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

optim.EP = function(blackbox, B,
                    equal=FALSE, ethresh=1e-2, 
                    Xstart=NULL, start=10, end=100, 
                    urate=10, rho=NULL, 
                    ncandf=function(k) { k }, 
                    dg_start=c(0.1*sqrt(nrow(B)), 1e-6), 
                    dlim=sqrt(nrow(B))*c(1e-3, 10), 
                    trcontrol=list(GLratio = c(1,1), sigma = 0.5^(1/nrow(B))/2, 
                                   maxsigma = 0.8^(1/nrow(B)), minsigma = 1/2^8,
                                   beta=0.5, alpha=2, kappa = 1e-4),
                    verb=2, ...)
{
  ## check start
  if(start >= end) stop("must have start < end")
  if(start == 0 & is.null(Xstart)) stop("must have start>0 or given Xstart")
  
  dim = nrow(B) # dimension
   
  # Check trcontrol
  if (is.null(trcontrol$sigma)) trcontrol$sigma = 0.5^(1/dim)/2
  trcontrol$sigma0 = trcontrol$sigma
  if (is.null(trcontrol$maxsigma)) trcontrol$maxsigma = 0.8^(1/dim)
  if (is.null(trcontrol$minsigma)) trcontrol$minsigma = 1/2^8
  if (is.null(trcontrol$beta)) trcontrol$beta = 0.5
  if (is.null(trcontrol$alpha)) trcontrol$alpha =  1/trcontrol$beta
  if (is.null(trcontrol$kappa)) trcontrol$kappa = 1e-4
  if (is.null(trcontrol$GLratio)) trcontrol$GLratio = c(1,1)
  
  max_global = trcontrol$GLratio[1]
  max_local = trcontrol$GLratio[2]
  if(max_global+max_local == 0) stop("must have max_global>0 or max_local>0")
  if (max_global == 0) {
    global_step <- FALSE
  } else if (max_local == 0){
    global_step <- TRUE
  } else {
    global_step <- TRUE
  }
  n_global = n_local = 0
  
  
  ## get initial design
  Hypercube = matrix(c(rep(0,dim), rep(1,dim)), ncol=2)                     # Hypercube [0, 1]^d
  if(!is.null(Xstart)){
    Xstart_unit = t(apply(Xstart, 1, function(x){(x-B[,1])/(B[,2]-B[,1])})) # Scale inputs to [0, 1]^d
    # Xstart_unit = (Xstart - matrix(rep(B[,1], nrow(Xstart)), byrow=TRUE, ncol=dim)) / matrix(rep((B[,2]-B[,1]), nrow(Xstart)), byrow=TRUE, ncol=dim)
  }else{
    Xstart_unit = Xstart = NULL
  }
  if(start>0){
    X_unit = dopt.gp(start, Xcand=lhs(10*start, Hypercube))$XX
    X = t(apply(X_unit, 1, function(x){x * (B[,2]-B[,1]) + B[,1]}))        # unscale inputs to B
    # X = matrix(rep(B[,1], start), byrow=TRUE, ncol=dim) + X_unit * matrix(rep((B[,2]-B[,1]), start), byrow=TRUE, ncol=dim)
  }else{
    X_unit = X = NULL
  }
  X_unit = rbind(Xstart_unit, X_unit)
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
  prog = rep(NA, start)                       # progress
  obj = rep(NA, start)                        # objective values
  feasibility = rep(NA, start)                # Whether the point is feasible?
  C = matrix(NA, nrow=start, ncol=nc)         # constraint values
  C_bilog = matrix(NA, nrow=start, ncol=nc)   # bi-log transformation of constraint values
  CV = matrix(NA, nrow=start, ncol=nc)        # constraint violations in terms of bi-log transformation
  
  obj[1] = out$obj; C[1, ] = out$c; C_bilog[1, ] = sign(C[1, ])*log(1+abs(C[1, ]))
  CV[1, !equal] = pmax(0, C_bilog[1, !equal]) # CV for inequality
  CV[1, equal] = abs(C_bilog[1, equal])       # CV for equality
  feasibility[1] = all(C[1,!equal] <= 0) && all(abs(C[1,equal]) <= ethresh)
  prog[1] = ifelse(feasibility[1], obj[1], Inf)
  
  ## remainder of starting run
  for(k in 2:start) { ## now that problem is vectorized we can probably remove for
    out = blackbox(X[k,], ...)
    obj[k] = out$obj; C[k,] = out$c; C_bilog[k, ] = sign(C[k, ])*log(1+abs(C[k, ]))
    CV[k, !equal] = pmax(0, C_bilog[k, !equal])
    CV[k, equal] = abs(C_bilog[k, equal])
    feasibility[k] = all(C[k,!equal] <= 0) && all(abs(C[k,equal]) <= ethresh)
    prog[k] = ifelse(feasibility[k] && obj[k] < prog[k-1], obj[k], prog[k-1]) 
  }
  
  ## handle initial rho value
  fmean = mean(obj); fsd = sd(obj)
  obj_norm = (obj-fmean)/fsd       # standard normalization on objective values
  # obj_norm = obj/max(abs(obj))

  ## handle initial rho values
  if(is.null(rho)){
    if(all(feasibility)){            # 
      rho = rep(1e-6, nc)
    }else {
      ECV = colMeans(CV[!feasibility,,drop=FALSE]) # averaged CV over infeasibile design points
      rho = mean(abs(obj_norm)) * ECV/sum(ECV^2)
    }
    if(any(equal)) rho[equal] = pmax(1/ethresh/sum(equal), rho[equal])
  }else{
    if(length(rho) != nc || any(rho < 0)) stop("rho should be a non-negative vector whose length is equal to the number of constraints.")
  }
  ## printing initial
  if(verb > 0){ cat("initial rho=(", paste(signif(rho,3), collapse=", "), ")\n", sep="") }
  
  ## calculate EP for data seen so far
  scv = CV%*%rho        # weighted sum of the constraint violations
  ep = obj_norm + scv
  epbest = min(ep)      # best EP seen so far
  m2 = prog[start]      # best feasible solution so far
  ## best solution so far
  if(is.finite(m2)){            # if at least one feasible solution was found
    xbest = X[which.min(prog),] 
  }else{                        # if none of the solutions are feasible
    xbest = X[which.min(scv),] 
  }
  
  ## initialize objective surrogate
  ab = darg(NULL, X_unit)$ab
  fgpi = newGPsep(X_unit, obj_norm, d = dg_start[1], g = dg_start[2], dK = TRUE)
  df = mleGPsep(fgpi, param = "d", tmin = dlim[1], tmax = dlim[2], ab = ab, verb = verb-1)$d
  
  ## initializing constraint surrogates
  Cgpi = rep(NA, nc)
  dc = matrix(NA, nrow=nc, ncol=dim)
  for (j in 1:nc) {
    Cgpi[j] = newGPsep(X_unit, C_bilog[,j], d=dg_start[1], g=dg_start[2], dK=TRUE)
    dc[j,] = mleGPsep(Cgpi[j], param = "d", tmin = dlim[1], tmax = dlim[2], ab = ab, verb=verb-1)$d
  }
  
  ## iterating over the black box evaluations
  for (k in (start+1):end) {
    ## rebuild surrogates periodically under new normalized responses
    if(k > (start+1) && (k %% urate == 0)) {
      ## objective surrogate
      deleteGPsep(fgpi)
      df[df<dlim[1]] = 10*dlim[1]
      df[df>dlim[2]] = dlim[2]/10
      fmean = mean(obj); fsd = sd(obj)
      obj_norm = (obj-fmean)/fsd
      # obj_norm = obj/max(abs(obj))
      fgpi = newGPsep(X_unit, obj_norm, d=df, g=dg_start[2], dK=TRUE)
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
    if(is.finite(m2)){
      ck = C[which.min(ep),]
      invalid = rep(NA, nc)
      invalid[!equal] = (ck[!equal] > 0); invalid[equal] = (abs(ck[equal]) > ethresh)
      nupdate = 0
      while(any(invalid) && nupdate<=5){
        rho[invalid] = rho[invalid]*2
        scv = CV%*%rho
        ep = obj_norm + scv
        ck = C[which.min(ep),]
        invalid[!equal] = (ck[!equal] > 0); invalid[equal] = (abs(ck[equal]) > ethresh)
        nupdate = nupdate + 1
      }
      epbest = min(ep)
    }
    
    
    ## random candidate grid
    ncand = ncandf(k)
    if(global_step){
      if(verb > 0){ cat("k=", k, " global step ", "\n", sep="") }
      cands = lhs(ncand, Hypercube)
      TRlower = Hypercube[, 1]
      TRupper = Hypercube[, 2]
    }else{
      if(verb > 0) { cat("k=", k, " local step ", "\n", sep="") }
      ## Trust region space
      # weights = df / (prod(df))^(1/length(df))
      # trcontrol$sigma = trcontrol$sigma * weights 
      xbest_unit = (xbest-B[,1])/(B[,2]-B[,1])
      TRlower = pmax(xbest_unit - trcontrol$sigma, 0)
      TRupper = pmin(xbest_unit + trcontrol$sigma, 1)
      TRspace = cbind(TRlower, TRupper)
      
      ## random candidate grid
      if(dim <= 20){
        cands = lhs(ncand, TRspace)
      }else{
        pert = lhs(ncand, TRspace)
        # create a perturbation mask
        prob_perturb = min(20/dim, 1)
        mask = matrix(runif(ncand*dim), ncol = dim) <= prob_perturb
        ind_mask = which(rowSums(mask) == 0)
        mask[ind_mask, sample.int(dim, size = length(ind_mask), replace = TRUE)] = 1
        cands = matrix(rep(xbest, ncand), nrow = ncand, byrow = TRUE)
        cands[which(mask==1)] = pert[which(mask==1)]
      }
    }
    
    if(verb > 0 && !global_step) { 
      cat("trcontrol: sigma=", trcontrol$sigma, " ", sep="") 
      cat("TRlower=[", paste(signif(TRlower*(B[,2]-B[,1])+B[,1],3), collapse=" "), sep="")
      cat("]; TRupper=[", paste(signif(TRupper*(B[,2]-B[,1])+B[,1],3), collapse=" "), "]\n", sep="")
    }
    
    ## calculate composite surrogate, and evaluate SEI and/or EY
    by = "sei"
    AF = EP_AcqFunc(cands, fgpi, Cgpi, epbest, rho, equal, eiey=by)
    m = which.max(AF)
    if(AF[m] > 1e-6){ # maximization expected improvement approach
      out = optim(par=cands[m, ], fn=EP_AcqFunc, method="L-BFGS-B",
                  lower=TRlower, upper=TRupper,
                  fgpi=fgpi, Cgpi=Cgpi, 
                  control = list(fnscale = -1), # maximization problem
                  eiey=by, epbest=epbest, rho=rho, equal=equal)
    }else{# Restart optimization with minimization predictive mean approach
      by = "ey"
      AF = EP_AcqFunc(cands, fgpi, Cgpi, epbest, rho, equal, eiey=by)
      m = which.min(AF)
      out = optim(par=cands[m, ], fn=EP_AcqFunc, method="L-BFGS-B",
                  lower=TRlower, upper=TRupper,
                  fgpi=fgpi, Cgpi=Cgpi, 
                  eiey=by, epbest=epbest, rho=rho, equal=equal)
    }

    ## calculate next point
    xnext_unit = matrix(out$par, nrow = 1)
    X_unit = rbind(X_unit, xnext_unit)
    xnext = xnext_unit * (B[,2]-B[,1]) + B[,1]
    X = rbind(X, xnext)
    
    ## new run
    out = blackbox(xnext, ...)
    fnext = out$obj; obj = c(obj, fnext); C = rbind(C, out$c)
    obj_norm = (obj-fmean)/fsd
    # obj_norm = obj/max(abs(obj))
    C_bilog = rbind(C_bilog, sign(out$c)*log(1+abs(out$c)))
    CV = rbind(CV, rep(NA, nc))
    CV[k, !equal] = pmax(0, C_bilog[k, !equal]) # constraint violations for inequality
    CV[k, equal] = abs(C_bilog[k, equal]) # constraint violations for equality
    
    ## check if best valid has changed
    feasibility = c(feasibility, all(out$c[!equal] <= 0) && all(abs(out$c[equal]) <= ethresh))
    if(feasibility[k] && fnext < prog[k-1]) {
      m2 = fnext
    } # otherwise m2 unchanged; should be the same as prog[k-1]
    prog = c(prog, m2)
    
    ## rho update
    if(all(feasibility)){ # 
      rho_new = rep(1e-6, nc)
    }else {
      ECV = colMeans(CV[!feasibility,,drop=FALSE])
      rho_new = mean(abs(obj_norm)) * ECV/sum(ECV^2)
    }
    if(verb > 0 && any(rho_new > rho)){ # printing progress
      cat("updating rho=(", paste(signif(pmax(rho_new, rho),3), collapse=", "), ")\n", sep="")
    }
    rho = pmax(rho_new, rho)
    
    ## calculate EP for data seen so far
    scv  = CV%*%rho
    ep = obj_norm + scv
    epbest = min(ep)
    if(is.finite(m2)){ # best solution so far
      xbest = X[which.min(prog),] 
    }else{ 
      xbest = X[which.min(scv),] 
    }
    
    ## progress meter
    if(verb > 0) {
      cat(by, "=", AF[m], sep="")
      cat("; xnext=[", paste(signif(xnext,3), collapse=" "), sep="")
      cat("]; xbest=[", paste(signif(xbest,3), collapse=" "), sep="")
      cat("]; ybest (prog=", m2, ", ep=", epbest, ")\n", sep="")
    }
    
    ## sufficient decrease condition
    if(is.finite(m2)){ # if at least one valid solution was found
      if(feasibility[k] && fnext < (prog[k-1]-trcontrol$kappa*trcontrol$sigma^2)){success = TRUE}else{success = FALSE}
    }else{             # if none of the solutions are valid
      if(which.min(scv) == k){success = TRUE}else{success = FALSE}
    }
    
    
    # Update number of consecutive trial / local steps
    if(global_step){n_global = n_global + 1}else{n_local = n_local + 1}

    
    ## Update Trust-region scheme
    if (success) {
      ## Successfull step: update TR and go back to global step if needed
      if (global_step) {
        # If successful global step: keep global and reinitialize counters
        global_step = TRUE
        n_global = n_local = 0
      } else {
        # If successful local step: update TR and switch to global only if local steps exhausted
        trcontrol$sigma = min(trcontrol$alpha*trcontrol$sigma, trcontrol$maxsigma)
        if (n_local >= max_local) {
          global_step = TRUE
          n_global = n_local = 0
        }
      }
    } else {
      ## Unsuccessfull step
      if (global_step) {
        # Unsuccessfull global step: go to local step if max_global attained
        if (n_global >= max_global) {
          global_step = FALSE
          n_global = n_local = 0
        }
      } else {
        # Unsuccessfull local step: update TR and go back to global step if max_local attained
        trcontrol$sigma = trcontrol$beta*trcontrol$sigma
        if (n_local >= max_local) {
          global_step = TRUE
          n_global = n_local = 0
        }
      }
    }
    if(trcontrol$sigma < trcontrol$minsigma){
      # restart when trust region becomes too small
      trcontrol$sigma = trcontrol$sigma0
      # n_global = n_local = 0
      if(verb > 0) cat("restart when trust region becomes too small \n")
    }
    
    if(max_global == 0){global_step = FALSE}else if(max_local == 0){global_step = TRUE}
    
    ## update GP fits
    updateGPsep(fgpi, xnext_unit, obj_norm[k], verb = verb-2)
    df = mleGPsep(fgpi, param = "d", tmin = dlim[1], tmax = dlim[2], ab = ab, verb=verb-1)$d
    
    for(j in 1:nc) {
      updateGPsep(Cgpi[j], xnext_unit, C_bilog[k,j], verb = verb-2)
      dc[j,] = mleGPsep(Cgpi[j], param = "d",  tmin = dlim[1], tmax = dlim[2], ab = ab, verb=verb-1)$d
    }
  }
  
  ## delete GP surrogates
  deleteGPsep(fgpi)
  for(j in 1:nc) deleteGPsep(Cgpi[j])
  
  return(list(prog = prog, xbest = xbest, obj = obj, C=C, X = X, feasibility=feasibility, rho=rho))
}


EP_AcqFunc = function(x, fgpi, Cgpi, epbest, rho, equal, eiey="sei")
{
  if(is.null(nrow(x))) x = matrix(x, nrow=1)
  ncand = nrow(x)
  nc = length(Cgpi) # number of the constraint
  
  ## Acquaisition function
  if(eiey == "sei"){
    ## objective
    # if(is.finite(m2)){           # if at least one valid solution was found
      pred_f = predGPsep(fgpi, x, lite=TRUE)
      mu_f = pred_f$mean
      sigma_f = sqrt(pred_f$s2)
    # }else{                      # if none of the solutions are valid  
    #   mu_f = sigma_f = 0
    #   epbest = scvbest
    # }
    
    ## constraints
    mu_C = sigma_C = omega = matrix(NA, nc, ncand)
    for (j in 1:nc) {
      pred_C = predGPsep(Cgpi[j], x, lite=TRUE)
      mu_C[j,] = pred_C$mean
      sigma_C[j,] = sqrt(pred_C$s2)
      omega[j,] = (equal[j]+1)*pnorm(mu_C[j,]/sigma_C[j,]) - equal[j]
    }
    
    ## Acquaisition function
    mu_ep = mu_f + rho%*%(omega*mu_C)
    sigma_ep = sqrt(sigma_f^2 + (rho^2)%*%((omega*sigma_C)^2))
    d = (epbest - mu_ep)/sigma_ep
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
    # if(is.finite(m2)){           # if at least one valid solution was found
      pred_f = predGPsep(fgpi, x, lite=TRUE)
      mu_f = pred_f$mean
    # }else{                      # if none of the solutions are valid  
    #   mu_f = 0
    # }
    
    ## constraints
    mu_C = sigma_C = EV = matrix(NA, nc, nrow(x))
    
    for (j in 1:nc) {
      pred_C = predGPsep(Cgpi[j], x, lite=TRUE)
      mu_C[j,] = pred_C$mean
      sigma_C[j,] = sqrt(pred_C$s2)
      dC = mu_C[j,]/sigma_C[j,]
      EV[j,] = mu_C[j,]*((equal[j]+1)*pnorm(dC) - equal[j]) + (equal[j]+1)*sigma_C[j,]*dnorm(dC)
    }
    
    ## Acquaisition function
    AF = mu_f + rho%*%EV
  }else{
    stop("eiey should be 'sei' or 'ey'.")
  }
  
  return(AF)
}
