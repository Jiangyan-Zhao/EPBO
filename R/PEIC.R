#' @title PEIC: Pseudo expected improvement with constraints
#' 
#' @description Black-box optimization under inequality constraints via the PEIC.
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
#' @param nprl positive integer giving the number of the parallel points per iteration
#' @param urate positive integer indicating  how many optimization trials should pass before
#' each MLE/MAP update is performed for GP correlation lengthscale 
#' parameter(s) 
#' @param dg.start 2-vector giving starting values for the lengthscale and nugget parameters
#' of the GP surrogate model(s) for constraints
#' @param dlim 2-vector giving bounds for the lengthscale parameter(s) under MLE/MAP inference
#' @param plotprog \code{logical} indicating if the Pareto plots should be made after each inner iteration;
#' the plots show two panels tracking the best feasible objective values, 
#' and the updating points in a single iteration in parallel
#' @param verb a non-negative integer indicating the verbosity level; the larger the value the
#' more that is printed to the screen
#' @param ... additional arguments passed to \code{blackbox}
#'
#' @returns The output is a \code{list} summarizing the progress of the evaluations of the 
#' blackbox under optimization:
#' \item{prog }{ \code{vector} giving the best feasible (\code{g(x) <= 0 && |h(x)| <= ethresh}) value of the objective over the trials }
#' \item{xbest }{ \code{vector} giving the recommended solution}
#' \item{obj }{ \code{vector} giving the value of the objective for the input under consideration at each trial }
#' \item{C }{ \code{matrix} giving the value of the constraint function for the input under consideration at each trial}
#' \item{X }{ \code{matrix} giving the input values at which the blackbox function was evaluated }
#' 
#' @seealso \code{\link[laGP]{optim.auglag}}, \code{\link[laGP]{optim.efi}}, \code{\link[EPBO]{optim.BM}}, \code{\link[EPBO]{optim.EP}}
#'  
#' @author Jiangyan Zhao \email{zhaojy2017@126.com}
#' 
#' @references Qian, J., Y. Cheng, J. Zhang, J. Liu, and D. Zhan (2021). A parallel constrained efficient global 
#' optimization algorithm for expensive constrained optimization problems. \emph{Engineering Optimization} 53(2), 300–320.
#' 
#' @keywords optimize
#' @keywords design
#' 
#' @import laGP
#' @import tgp
#' @importFrom pso psoptim
#' @importFrom stats dnorm
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
#' PEIC = optim.PEIC(MTP, B, verb = 0)
#' # progress, best feasible value of the objective over the trials
#' PEIC$prog
#' # the recommended solution
#' PEIC$xbest


optim.PEIC = function(
    blackbox, B, nprl=3, start=10, end=100, 
    Xstart=NULL, urate=ceiling(10/nprl),
    dg_start=c(0.1*sqrt(nrow(B)), 1e-6), 
    dlim=c(1e-3, 10)*sqrt(nrow(B)), 
    plotprog=FALSE, verb=2, ...)
{
  ## check start
  if(start >= end) stop("must have start < end")
  if(start == 0 & is.null(Xstart)) stop("must have start>0 or given Xstart")
  
  dim = nrow(B) # dimension
  
  ## get initial design
  X = dopt.gp(start, Xcand=lhs(10*start, B))$XX
  ## X = lhs(start, B)
  X = rbind(Xstart, X)
  X_unit = normalize(X, B)
  start = nrow(X)
  
  ## first run to determine the number of constraints
  out = blackbox(X[1,], ...)
  nc = length(out$c)
  
  ## allocate progress objects, and initialize
  prog = rep(NA, start)                       # progress
  obj = rep(NA, start)                        # objective values
  feasibility = rep(NA, start)                # Whether the point is feasible?
  C = matrix(NA, nrow=start, ncol=nc)         # constraint values in terms of bi-log transformation
  obj[1] = out$obj; C[1, ] = out$c
  feasibility[1] = all(C[1,] <= 0)
  prog[1] = ifelse(feasibility[1], obj[1], Inf)
  
  ## remainder of starting run
  for(k in 2:start) { ## now that problem is vectorized we can probably remove for
    out = blackbox(X[k,], ...)
    obj[k] = out$obj; C[k,] = out$c
    feasibility[k] = all(C[k,] <= 0)
    prog[k] = ifelse(feasibility[k] && obj[k] < prog[k-1], obj[k], prog[k-1]) 
  }
  
  ## the best feasible objective value (BOFV) obtained so far
  m2 = prog[start]
  since = 0
  ## best solution so far
  if(is.finite(m2)){            # if at least one feasible solution was found
    xbest = X[which.min(prog),] 
  }else{                        # if none of the solutions are feasible
    xbest = B[,2]
  }
  # xbest_unit = as.vector(normalize(xbest, B))
  
  ## initialize objective surrogate
  fgpi = newGPsep(X_unit, obj, d = dg_start[1], g = dg_start[2], dK = TRUE)
  df = mleGPsep(fgpi, param = "d", tmin = dlim[1], tmax = dlim[2], verb = verb-1)$d
  df[df<=dlim[1]] = dlim[1] * 1.001
  df[df>=dlim[2]] = dlim[2] * 0.999
  
  ## initializing constraint surrogates
  Cgpi = rep(NA, nc)
  dc = matrix(NA, nrow=nc, ncol=dim)
  for (j in 1:nc) {
    Cgpi[j] = newGPsep(X_unit, C[,j], d=dg_start[1], g=dg_start[2], dK=TRUE)
  }
  
  ## printing initial design
  if(verb > 0) {
    cat("The initial design: ")
    cat("]; xbest=[", paste(signif(xbest,4), collapse=" "), sep="")
    cat("]; ybest (prog=", m2, ", since=", since, ")\n", sep="")
  }
  
  ## iterating over the black box evaluations
  for(k in seq(from = start, to = end-1, by = nprl)) {
    ## rebuild surrogates periodically under new normalized responses
    if(k > start && (ceiling((k-start)/nprl) %% urate == 0)) {
      ## objective surrogate
      deleteGPsep(fgpi)
      fgpi = newGPsep(X_unit, obj, d=df, g=dg_start[2], dK=TRUE)
      df = mleGPsep(fgpi, param = "d", tmin = dlim[1], tmax = dlim[2], verb=verb-1)$d
      df[df<=dlim[1]] = dlim[1] * 1.001
      df[df>=dlim[2]] = dlim[2] * 0.999
      
      ## constraint surrogates 
      for(j in 1:nc) {
        dc[j,] = mleGPsep(Cgpi[j], param = "d", tmin = dlim[1], tmax = dlim[2], verb=verb-1)$d
        dc[j, dc[j,]<=dlim[1]] = dlim[1] * 1.001
        dc[j, dc[j,]>=dlim[2]] = dlim[2] * 0.999
        deleteGPsep(Cgpi[j])
        Cgpi[j] = newGPsep(X_unit, C[,j], d=dc[j,], g=dg_start[2], dK=TRUE)
      }
    }
    
    
    ## Maximize the PEIC AF via the PSO algorithm.
    m2_old = m2
    for (prl in 1:nprl) {
      out_AF = psoptim(
        rep(NA,dim), fn = AF_PEIC,
        fgpi=fgpi, Cgpi=Cgpi, fmin=m2_old,
        df=df, point_update=tail(X_unit, prl-1),
        lower = 0, upper = 1,
        control = list(maxit = 100,   # generations
                       s = 100,       # swarm size
                       fnscale = -1)) # for maximization
      
      ## calculate next point
      xnext_unit = matrix(out_AF$par, nrow = 1)
      X_unit = rbind(X_unit, xnext_unit)
      xnext = unnormalize(xnext_unit, B)
      X = rbind(X, xnext)
      
      ## new run
      out = blackbox(xnext, ...)
      fnext = out$obj; obj = c(obj, fnext)
      C = rbind(C, out$c)
      
      ## check if the BFOV has changed
      feasibility = c(feasibility, all(out$c <= 0))
      since = since + 1
      if(tail(feasibility,1) && fnext < tail(prog,1)) {
        m2 = fnext; since = 0
      } # otherwise m2 unchanged; should be the same as prog[k-1]
      prog = c(prog, m2)
      
      ## best solution so far
      if(is.finite(m2)){            # if at least one feasible solution was found
        xbest = X[which.min(prog),] 
      }else{                        # if none of the solutions are feasible
        xbest = B[,2]
      }
      # xbest_unit = as.vector(normalize(xbest, B))
      
      ## progress meter
      if(verb > 0) {
        cat("k=", k+prl, " ", sep="")
        cat("PEIC=", paste(signif(out_AF$value,4), collapse=" "))
        cat("; xnext ([", paste(signif(xnext,4), collapse=" "), 
            "], feasibility=", feasibility[k], ")\n", sep="")
        cat(" xbest=[", paste(signif(xbest,4), collapse=" "), sep="")
        cat("]; ybest (prog=", paste(signif(m2,4), collapse=" ")) 
        cat(", since=", since, ")\n", sep="")
      }
    }

    
    ## plot progress
    if (plotprog) {
      par(ps=16, mfrow=c(1,2))
      ## progress
      if(is.finite(m2)){
        plot(prog, type="l", lwd=1.6, 
             xlim=c(start, k+nprl),
             ylim=range(prog[max(which(is.finite(prog))[1], start):(k+nprl)]), 
             xlab="n", ylab="BFOV", main="progress")
      }else{
        plot(prog, type="l", lwd=1.6, 
             xlim=c(start, k+nprl),
             ylim=range(obj[start:(k+nprl)]), 
             xlab="n", ylab="BFOV", main="progress")
      }
      
      # the updating points per iteration in parallel
      plot(tail(X[,1:2], nprl), 
           xlim = c(B[1,1], B[1,2]+0.1*abs(B[1,2])), ylim = B[2,], 
           pch = 18, col="red",
           xlab = "x1", ylab = "x2", main="updating points")
      text(tail(X[,1:2], nprl), labels = seq_len(nprl), pos=4)
      points(X[1:k,1:2], pch = 16, cex = 0.6)
      par(mfrow=c(1,1))
    }
    
    ## update GP fits
    updateGPsep(fgpi, tail(X_unit, nprl), tail(obj,nprl), verb = 0)
    for(j in 1:nc){
      updateGPsep(Cgpi[j], tail(X_unit, nprl), tail(C[,j],nprl), verb = 0)
    }
  }
  
  ## delete GP surrogates
  deleteGPsep(fgpi)
  for(j in 1:nc) deleteGPsep(Cgpi[j])
  
  return(list(prog = prog, xbest = xbest, obj = obj, C=C, X = X))
}
