#' @title PEPBO: Parallel Bayesian Optimization via Exact Penalty
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
#' @param nprl positive integer giving the number of the parallel points per iteration
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
#' @param dg_start 2-vector giving starting values for the lengthscale and nugget parameters
#' of the GP surrogate model(s) for the objective and constraints
#' @param dlim 2-vector giving bounds for the lengthscale parameter(s) under MLE inference
#' @param plotPareto \code{logical} indicating if progress plots should be made after each inner iteration;
#' the plots show four panels tracking the best feasible objective values, the ScaledEI or EY surface, 
#' the predictive mean and standard deviation of the objective function, over the first two input variables
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
#' @keywords optimize
#' @keywords design
#' 
#' @import laGP
#' @import tgp
#' @import mco
#' @importFrom pso psoptim
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

optim.PEP = function(
    blackbox, B, nprl=3, 
    equal=FALSE, ethresh=1e-2, 
    Xstart=NULL, start=10, end=100, urate=ceiling(10/nprl), rho=NULL, 
    dg_start=c(1e-2*sqrt(nrow(B)), 1e-6), 
    dlim=c(1e-4, 1)*sqrt(nrow(B)), 
    plotPareto=FALSE,verb=2, ...)
{
  ## check start
  if(nprl < 3) stop("must have nprl > 3")
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
  for(k in seq(from = start, to = end-1, by = nprl)) {
    ## rebuild surrogates periodically under new normalized responses
    if(k > start && (ceiling((k-start)/nprl) %% urate == 0)) {
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
    if(is.finite(m2) && any(invalid)){
      rho[invalid] = rho[invalid] * 2
      scv = CV%*%rho; ep = obj + scv; epbest = min(ep)
      if(verb > 0) cat(" the smallest EP is not feasible, updating rho=(", 
                       paste(signif(rho,3), collapse=", "), ")\n", sep="")
    }
    
    ## Approximate the Pareto front and Pareto set via the NSGA-II algorithm.
    AF_Pareto = nsga2(
      fn=AF_EYvsSDY, idim=dim, odim=2,
      fgpi=fgpi, fmean=fmean, fsd=fsd, Cgpi=Cgpi,
      rho=rho, equal=equal,
      generations=100, popsize=100*nprl,
      cprob=0.9, cdist=20, mprob=0.1, mdist=20,
      lower.bounds=rep(0, dim), upper.bounds=rep(1, dim))
    AF_PF = AF_Pareto$value # Pareto front (-logSDY, EY)
    AF_PS = AF_Pareto$par   # Pareto set
    PF_selected = matrix(NA, nrow = nprl, ncol = 2)
    
    ## Remove Pareto sets with probability of feasibility less than 0.5
    PS_PoF_idx = which(AF_PoF(AF_PS, Cgpi) > 0.5)
    if(length(PS_PoF_idx) > 0.2 * nrow(AF_PS)){
      AF_PF = AF_PF[PS_PoF_idx,]
      AF_PS = AF_PS[PS_PoF_idx,]
    }
    
    ## pure exploitation (minimize the predictive mean)
    min_EY_idx = which.min(AF_PF[,2])
    PF_selected[1,] = AF_PF[min_EY_idx, ]
    out_EY = optim( # enhance the EY approach
      par=AF_PS[min_EY_idx, ], fn=AF_EY, method="L-BFGS-B",
      lower=0, upper=1,
      fgpi=fgpi, Cgpi=Cgpi, fmean=fmean, fsd=fsd, 
      rho=rho, equal=equal)
    
    # calculate next point
    # xnext_unit = AF_PS[min_EY_idx, ]
    xnext_unit = matrix(out_EY$par, nrow = 1)
    X_unit = rbind(X_unit, xnext_unit)
    xnext = unnormalize(xnext_unit, B)
    X = rbind(X, xnext)
    
    ## Perform k-means clustering on the AF's Pareto front
    AF_PF_Scaled = min_max_scale(AF_PF)
    AF_cl = kmeans(x = AF_PF_Scaled, centers = nprl-1, 
                   iter.max = 100, nstart = 5*nprl)

    for (cl in 1:(nprl-1)) {
      ## a single member with highest EIC is selected from each cluster
      cl_idx = which(AF_cl$cluster == cl)
      cl_PS = AF_PS[cl_idx, , drop=FALSE]
      cl_PF = AF_PF[cl_idx, , drop=FALSE]
      
      
      # ## center of the Pareto front
      # cl_center = AF_cl$centers[cl,]
      # dist = apply(AF_PF_Scaled[cl_idx, , drop=FALSE], 1, 
      #              function(x){norm(x- cl_center, type="2")}) # Euclidean distance
      # min_dist_idx = which.min(dist)  # Find index of minimum distance
      # PF_selected[cl+1,] = cl_PF[min_dist_idx,]
      
      ## max std criterion
      max_STD_idx = which.max(-cl_PF[,1])
      PF_selected[cl+1,] = cl_PF[max_STD_idx,]
      
      
      # calculate next point
      # xnext_unit = cl_PS[min_dist_idx, ]
      xnext_unit = cl_PS[max_STD_idx, ]
      X_unit = rbind(X_unit, xnext_unit)
      xnext = unnormalize(xnext_unit, B)
      X = rbind(X, xnext)
    }
    
    ## new runs
    for(cl in 1:nprl) { ## now that problem is vectorized we can probably remove for
      out = blackbox(X[k+cl,], ...)
      fnext = out$obj; obj = c(obj, fnext); C = rbind(C, out$c)
      C_bilog = rbind(C_bilog, bilog(out$c))
      CV = rbind(CV, rep(NA, nc))
      CV[k+cl, !equal] = pmax(0, C_bilog[k+cl, !equal])
      CV[k+cl, equal] = abs(C_bilog[k+cl, equal])
      
      ## check if best valid has changed
      feasibility = c(feasibility, all(out$c[!equal] <= 0) && all(abs(out$c[equal]) <= ethresh))
      since = since + 1
      if(feasibility[k+cl] && fnext < prog[k+cl-1]) {
        m2 = fnext; since = 0
      } # otherwise m2 unchanged; should be the same as prog[k-1]
      prog = c(prog, m2)
      
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
        cat("k=", k+cl, " ", sep="")
        cat("; xnext ([", paste(signif(X[k+cl,],3), collapse=" "), 
            "], feasibility=", feasibility[k+cl], ")\n", sep="")
        cat(" xbest=[", paste(signif(xbest,3), collapse=" "), sep="")
        cat("]; ybest (prog=", paste(signif(m2,3), collapse=" ")) 
        cat(", ep=", paste(signif(epbest,3), collapse=" "))
        cat(", since=", since, ")\n", sep="")
      }
    }

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
    
    ## plot progress, Pareto Front, and Pareto set
    if(plotPareto) {
      par(ps=16, mfrow=c(2,2))
      ## progress
      if(is.finite(m2)){
        plot(prog, type="l", lwd=1.6, 
             xlab="n", ylab="BFOV", main="progress")
      }else{
        plot(prog, type="l", ylim=range(obj), lwd=1.6, 
             xlab="n", ylab="BFOV", main="progress")
      }
      
      ## legend
      plot.new()
      legend("center", bty="n",
             # title="legend in parameter space",
             legend=c(
               "Pareto set",
               "previous points", 
               paste0("selected updated (q=", nprl, ")")),
             pch=c(4,16,18), 
             pt.cex = c(1,1.2,1.5), pt.lwd = c(2,1,1))
      
      # Pareto Front
      plot(AF_PF, 
           pch = 4, col = AF_cl$cluster+1,
           lwd = 1.5, cex = 0.6,
           xlab="-logSDY", ylab="EY", main="Objective space")
      points(PF_selected, 
             col = 1:nprl, pch = 18, cex = 2)
      
      # Pareto Set
      plot(unnormalize(AF_PS[,1:2], B), 
           xlim = B[1,], ylim = B[2,], 
           pch = 4, col = AF_cl$cluster+1,
           lwd = 1.5, cex = 0.6,
           xlab = "x1", ylab = "x2", main="Parameter space")
      points(tail(X[,1:2], nprl), 
             col = 1:nprl, pch = 18, cex = 2)
      points(X[1:k,1:2], pch = 16, cex = 0.6)
      par(mfrow=c(1,1))
    }
    
    ## update GP fits
    updateGPsep(fgpi, tail(X_unit, nprl), (tail(obj,nprl)-fmean)/fsd, verb = 0)
    df = mleGPsep(fgpi, param = "d", tmin = dlim[1], tmax = dlim[2], ab = ab, verb = 0)$d
    sigma_f = sqrt(predGPsep(fgpi, X_unit, lite=TRUE)$s2) #* fsd
    if(median(sigma_f) > 0.1){# rebuild GP as sigma is large at observed points
      if(verb > 0) cat(" rebuild objective surrogate as sigma is large at observed points \n")
      deleteGPsep(fgpi)
      fgpi = newGPsep(X_unit, (obj-fmean)/fsd, d = dlim[1], g = dg_start[2], dK=TRUE)
      df = mleGPsep(fgpi, param = "d", tmin = dlim[1], tmax = dlim[2], ab = ab, verb = 0)$d
    }
    
    for(j in 1:nc) {
      updateGPsep(Cgpi[j], tail(X_unit, nprl), tail(C_bilog[,j],nprl), verb = 0)
      dc[j,] = mleGPsep(Cgpi[j], param = "d",  tmin = dlim[1], tmax = dlim[2], ab = ab, verb = 0)$d
      sigma_c = sqrt(predGPsep(Cgpi[j], X_unit, lite=TRUE)$s2)
      if(median(sigma_c) > 0.1){ # rebuild GP as sigma is large at observed points
        if(verb > 0) cat(" rebuild constrained", j, "surrogate as sigma is large at observed points \n")
        deleteGPsep(Cgpi[j])
        Cgpi[j] = newGPsep(X_unit, C_bilog[,j], d = dlim[1], g = dg_start[2], dK=TRUE)
        dc[j,] = mleGPsep(Cgpi[j], param = "d",  tmin = dlim[1], tmax = dlim[2], ab = ab, verb = 0)$d
      }
    }
  }
  
  ## delete GP surrogates
  deleteGPsep(fgpi)
  for(j in 1:nc) deleteGPsep(Cgpi[j])
  
  return(list(prog = prog, xbest = xbest, 
              obj = obj, C=C, X = X, 
              feasibility=feasibility, rho=rho)) #AF_time = AF_time
}
