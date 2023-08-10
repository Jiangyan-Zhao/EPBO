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
#' @param ab prior parameters; see \code{\link{darg}} describing the prior used on the
#' lengthscale parameter during emulation(s) for the constraints
#' @param dlim 2-vector giving bounds for the lengthscale parameter(s) under MLE/MAP inference,
#' thereby augmenting the prior specification in \code{ab}
#' @param TR_control an optional list of control parameters of the trust region. See "Details".
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
                    ey_tol=0.01,
                    dg_start=c(0.1, 1e-6), ab=c(3/2, 8), 
                    dlim=sqrt(nrow(B))*c(1e-3, 10), 
                    # TR_control=list(length=0.8, length0=0.8,
                    #                 length_min=0.5^7, length_max=1.6, 
                    #                 failure_counter=0, failure_tolerance=nrow(B), 
                    #                 success_counter=0, success_tolerance=nrow(B)),
                    verb=2, ...)
{
  ## check start
  if(start >= end) stop("must have start < end")
  
  dim = nrow(B) # dimension
   
  # if (is.null(TR_control$length)) TR_control$length = 0.8
  # if (is.null(TR_control$length0)) TR_control$length0 = 0.8
  # if (is.null(TR_control$length_min)) TR_control$length_min = 0.5^7
  # if (is.null(TR_control$length_max)) TR_control$length_max = 1.6
  # if (is.null(TR_control$failure_counter)) TR_control$failure_counter = 0
  # if (is.null(TR_control$failure_tolerance)) TR_control$failure_tolerance = dim
  # if (is.null(TR_control$success_counter)) TR_control$success_counter = 0
  # if (is.null(TR_control$success_tolerance)) TR_control$success_tolerance = dim
  
  
  ## get initial design
  Hypercube = matrix(c(rep(0,dim), rep(1,dim)), ncol=2)             # Hypercube
  if(start>0){
    X_unit = dopt.gp(start, Xcand=lhs(10*start, Hypercube))$XX
    X = t(apply(X_unit, 1, function(x){x * (B[,2]-B[,1]) + B[,1]}))
  }else{
    X_unit = X = NULL
  }
  if(!is.null(Xstart)){
    X_unit = rbind(t(apply(Xstart, 1, function(x){(x-B[,1])/(B[,2]-B[,1])})), X_unit)
    X = rbind(Xstart, X)
  }
  start = nrow(X)
  if(start == 0) stop("must have nrow(Xstart)+start>0")
  
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
    ECV = colMeans(CV)            # averaged CV
    if(all(ECV == 0)){            # 
      rho = rep(0, nc)
    }else {
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
  ybest = min(ep)       # best EP seen so far
  scv_best = min(scv) # best scv seen so far
  m2 = prog[start]      # best feasible solution so far
  ## best solution so far
  if(is.finite(m2)){            # if at least one feasible solution was found
    xbest = X[which.min(prog),] 
  }else{                        # if none of the solutions are feasible
    xbest = X[which.min(scv),] 
  }
  
  ## initialize objective surrogate
  fgpi = newGPsep(X_unit, obj_norm, d = dg_start[1], g = dg_start[2], dK = TRUE)
  df = mleGPsep(fgpi, param = "d", tmin = dlim[1], tmax = dlim[2], ab = ab, verb = verb-1)$d
  
  ## initializing constraint surrogates
  Cgpi = rep(NA, nc)
  dc = matrix(NA, nrow=nc, ncol=nrow(B))
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
      ybest = min(ep)
    }
    
    # ## Trust region
    # weights = df / (prod(df))^(1/length(df))
    # TR_length = TR_control$length * (B[,2]-B[,1]) * weights 
    # TR_lb = pmin(pmax(B[,1], xbest - TR_length/2), B[,2])
    # TR_ub = pmin(pmax(B[,1], xbest + TR_length/2), B[,2])
    # TR_B = cbind(TR_lb, TR_ub)
    # if(verb > 0) {
    #   cat("k=", k, " ", sep="")
    #   cat("TR_lb=[", paste(signif(TR_lb,3), collapse=" "), sep="")
    #   cat("]; TR_ub=[", paste(signif(TR_ub,3), collapse=" "), "]\n", sep="")
    # }
    
    ## random candidate grid
    ncand = ncandf(k)
    cands = lhs(ncand, Hypercube)
    # pert = lhs(ncand, TR_B)
    # # create a perturbation mask
    # prob_perturb = min(20/nrow(B), 1)
    # mask = matrix(runif(ncand*nrow(B)), ncol = nrow(B)) <= prob_perturb
    # ind_mask = which(rowSums(mask) == 0)
    # mask[ind_mask, sample.int(nrow(B), size = length(ind_mask), replace = TRUE)] = 1
    # cands = matrix(rep(xbest, ncand), nrow = ncand, byrow = TRUE)
    # cands[which(mask==1)] = pert[which(mask==1)]
    
    ## calculate composite surrogate, and evaluate SEI and/or EY
    AF = EP_AcqFunc(cands, fgpi, Cgpi, ybest, rho, equal, eiey="sei")
    nzsei = sum(AF > sqrt(.Machine$double.eps))
    if(nzsei > 0.01*ncand){ # maximization expected improvement approach in TR_B
      m = which.max(AF)
      out = optim(par=cands[m, ], fn=EP_AcqFunc, method="L-BFGS-B",
                  lower=rep(0,dim), upper=rep(1,dim),
                  fgpi=fgpi, Cgpi=Cgpi, 
                  control = list(fnscale = -1), # maximization problem
                  eiey="sei", ybest=ybest, rho=rho, equal=equal)
      by = "sei"
    }else{
      # if(verb > 0) cat("TR_B to B \n")
      # cands_big = lhs(ncand, B)
      # AF = EP_AcqFunc(cands_big, fgpi, fmean, fsd, Cgpi, Cmean, Csd, ybest, rho, 
      #                 equal, blackbox, eiey="ei")
      # nzsei = sum(AF > sqrt(.Machine$double.eps))
      # if(nzsei > 0.01*ncand){ # maximization expected improvement approach in B
      #   m = which.max(AF)
      #   out = optim(par=cands_big[m, ], fn=EP_AcqFunc, method="L-BFGS-B",
      #               lower=B[,1], upper=B[,2], 
      #               control = list(fnscale = -1), # maximization problem
      #               fgpi=fgpi, fmean=fmean, fsd=fsd,
      #               Cgpi=Cgpi, Cmean=Cmean, Csd=Csd,
      #               eiey="ei", ybest=ybest, rho=rho, equal=equal, 
      #               blackbox=blackbox)
      #   by = "ei"
      # }else{# minimization predictive mean approach in TR_B
        AF = EP_AcqFunc(cands, fgpi, Cgpi, ybest, rho, equal, eiey="ey")
        m = which.min(AF)
        out = optim(par=cands[m, ], fn=EP_AcqFunc, method="L-BFGS-B",
                    lower=rep(0,dim), upper=rep(1,dim),
                    fgpi=fgpi, Cgpi=Cgpi, 
                    eiey="ey",ybest=ybest, rho=rho, equal=equal)
        by = "ey"
      # }
    }
    
    
    ## calculate next point
    xnext_unit = matrix(out$par, nrow = 1)
    X_unit = rbind(X_unit, xnext_unit)
    xnext = xnext_unit * (B[,2]-B[,1]) + B[,1]
    X = rbind(X, xnext)
    
    ## new run
    out = blackbox(xnext, ...)
    fnext = out$obj; obj = c(obj, fnext); 
    obj_norm = (obj-fmean)/fsd
    # obj_norm = obj/max(abs(obj))
    C = rbind(C, out$c); C_bilog = rbind(C_bilog, sign(out$c)*log(1+abs(out$c)))
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
    ECV = colMeans(CV) # the constraint violation averaged over the design
    if(all(ECV == 0)){ # 
      rho_new = rep(0, nc)
    }else {
      rho_new = mean(abs(obj_norm)) * ECV/sum(ECV^2)
    }
    if(verb > 0 && any(rho_new > rho)){ # printing progress
      cat("updating rho=(", paste(signif(pmax(rho_new, rho),3), collapse=", "), ")\n", sep="")
    }
    rho = pmax(rho_new, rho)
    
    ## calculate EP for data seen so far
    scv  = CV%*%rho
    ep = obj_norm + scv
    ybest = min(ep)
    scv_best = min(scv)
    if(is.finite(m2)){ # best solution so far
      xbest = X[which.min(prog),] 
    }else{ 
      xbest = X[which.min(scv),] 
    }
    
    ## progress meter
    if(verb > 0) {
      cat(by, "=", AF[m], " (", nzsei,  "/", ncand, ")", sep="")
      cat("; xnext=[", paste(signif(xnext,3), collapse=" "), sep="")
      cat("]; xbest=[", paste(signif(xbest,3), collapse=" "), sep="")
      cat("]; ybest (prog=", m2, ", ep=", ybest, ")\n", sep="")
    }
    
    # ## update TR_control
    # if(is.finite(m2)){ # if at least one valid solution was found
    #   if(which.min(prog) == k){
    #     TR_control$success_counter = TR_control$success_counter + 1
    #     TR_control$failure_counter = 0
    #   }else{
    #     TR_control$failure_counter = TR_control$failure_counter + 1
    #     TR_control$success_counter = 0
    #   }
    # }else{             # if none of the solutions are valid
    #   if(which.min(scv) == k){
    #     TR_control$success_counter = TR_control$success_counter + 1
    #     TR_control$failure_counter = 0
    #   }else{
    #     TR_control$failure_counter = TR_control$failure_counter + 1
    #     TR_control$success_counter = 0
    #   }
    # }
    # ## update the trust region
    # if(TR_control$success_counter == TR_control$success_tolerance){
    #   # expand trust region
    #   TR_control$length = min(TR_control$length * 2, TR_control$length_max)
    #   TR_control$success_counter = 0
    # }
    # if(TR_control$failure_counter == TR_control$failure_tolerance){
    #   # shrink trust region
    #   TR_control$length = TR_control$length / 2
    #   TR_control$failure_counter = 0
    # }
    # if(TR_control$length < TR_control$length_min){
    #   # restart when trust region becomes too small
    #   TR_control$length = TR_control$length0
    #   TR_control$success_counter = 0
    #   TR_control$failure_counter = 0
    #   if(verb > 0) cat("restart when trust region becomes too small \n")
    # }
    # if(verb > 0) {
    #   cat("TR_control: length=", TR_control$length,
    #       "; success_counter=", TR_control$success_counter,
    #       "; failure_counter=", TR_control$failure_counter, "\n", sep="")
    # }
    
    
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


EP_AcqFunc = function(x, fgpi, Cgpi, ybest, rho, equal, eiey="ei")
{
  if(is.null(nrow(x))) x = matrix(x, nrow=1)
  ncand = nrow(x)
  nc = length(Cgpi) # number of the constraint
  
  ## Acquaisition function
  if(eiey == "sei"){
    ## objective
    pred_f = predGPsep(fgpi, x, lite=TRUE)
    mu_f = pred_f$mean
    sigma_f =  sqrt(pred_f$s2)
    
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
    mu_f = pred_f$mean
    
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
