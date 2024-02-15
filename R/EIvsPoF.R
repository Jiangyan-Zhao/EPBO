#' @title EI versus PoF
#' 
#' @description Black-box optimization under inequality constraints via a bi-objective optimization.
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
#' @param ab prior parameters; see \code{\link{darg}} describing the prior used on the
#' lengthscale parameter during emulation(s) for the constraints
#' @param dlim 2-vector giving bounds for the lengthscale parameter(s) under MLE/MAP inference,
#' thereby augmenting the prior specification in \code{ab}
#' @param plotPareto \code{logical} indicating if the Pareto plots should be made after each inner iteration;
#' the plots show two panels tracking the Pareto front, and the Pareto set
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
#' @references Parr, J. M., A. J. Keane, A. I. Forrester, and C. M. Holden (2012). Infill sampling criteria
#' for surrogate-based optimization with constraint handling. \emph{Engineering Optimization} 44(10), 1147–1166.
#' 
#' @keywords optimize
#' @keywords design
#' 
#' @import laGP
#' @import tgp
#' @import mco
#' @importFrom stats dnorm
#' @importFrom stats kmeans
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
#' EIvsPoF = optim.EIvsPoF(MTP, B, verb = 0)
#' # progress, best feasible value of the objective over the trials
#' EIvsPoF$prog
#' # the recommended solution
#' EIvsPoF$xbest


optim.EIvsPoF = function(
    blackbox, B, nprl=3, start=10, end=100, 
    Xstart=NULL, urate=ceiling(10/nprl),
    dg_start=c(1e-2*sqrt(nrow(B)), 1e-6), 
    dlim=c(1e-4, 1)*sqrt(nrow(B)), 
    plotPareto=FALSE, verb=2, ...)
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
  obj[1] = out$obj; C[1, ] = bilog(out$c)
  feasibility[1] = all(C[1,] <= 0)
  prog[1] = ifelse(feasibility[1], obj[1], Inf)

  ## remainder of starting run
  for(k in 2:start) { ## now that problem is vectorized we can probably remove for
    out = blackbox(X[k,], ...)
    obj[k] = out$obj; C[k,] = bilog(out$c)
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
    Cgpi[j] = newGPsep(X_unit, C[,j], d=dg_start[1], g=dg_start[2], dK=TRUE)
    dc[j,] = mleGPsep(Cgpi[j], param = "d", tmin = dlim[1], tmax = dlim[2], ab = ab, verb=verb-1)$d
    deleteGPsep(Cgpi[j])
    dc[j, dc[j,]<dlim[1]] = 10*dlim[1]
    dc[j, dc[j,]>dlim[2]] = dlim[2]/10
    Cgpi[j] = newGPsep(X_unit, C[,j], d=dc[j,], g=dg_start[2], dK=TRUE)
    dc[j,] = mleGPsep(Cgpi[j], param = "d",  tmin = dlim[1], tmax = dlim[2], ab = ab, verb=verb-1)$d
  }
  
  ## printing initial design
  if(verb > 0) {
    cat("The initial design: ")
    cat("ab=[", paste(signif(ab,3), collapse=", "), sep="")
    cat("]; xbest=[", paste(signif(xbest,3), collapse=" "), sep="")
    cat("]; ybest (prog=", m2, ", since=", since, ")\n", sep="")
  }
  
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
        Cgpi[j] = newGPsep(X_unit, C[,j], d=dc[j,], g=dg_start[2], dK=TRUE)
        dc[j,] = mleGPsep(Cgpi[j], param = "d",  tmin = dlim[1], tmax = dlim[2], ab = ab, verb=verb-1)$d
      }
    }
    
    
    ## Approximate the Pareto front and Pareto set via the NSGA-II algorithm.
    AF_Pareto = nsga2(
      fn=AF_EIvsPoF, idim=dim, odim=2,
      fgpi=fgpi, fmean=fmean, fsd=fsd, Cgpi=Cgpi, fmin=m2,
      generations=100, popsize=100*nprl,
      cprob=0.9, cdist=20, mprob=0.1, mdist=20,
      lower.bounds=rep(0, dim), upper.bounds=rep(1, dim))
    AF_PF = AF_Pareto$value # Pareto front
    AF_PS = AF_Pareto$par   # Pareto set
    
    
    ## Perform k-means clustering on the AF's Pareto front
    # Here, a small disturbance is added to PF to avoid PS aggregation, 
    # resulting in clustering failure.
    AF_PF_dis = AF_PF + rnorm(length(AF_PF), mean = 0, sd = 0.001)
    AF_cl = kmeans(x = AF_PF_dis, centers = nprl, 
                   iter.max = 100, nstart = 5*nprl)
    for (cl in 1:nprl) {
      ## a single member with highest EIC is selected from each cluster
      cl_idx = which(AF_cl$cluster == cl)
      cl_PS = matrix(AF_PS[cl_idx, ], nrow = length(cl_idx))
      cl_PF = matrix(AF_PF[cl_idx, ], nrow = length(cl_idx))
      # cl_EIC = rowSums(-cl_PF)
      cl_EIC = exp(-cl_PF[,1]) * (-cl_PF[,2])
      max_EIC_idx = which.max(cl_EIC)
      
      ## calculate next point
      xnext_unit = cl_PS[max_EIC_idx, ]
      X_unit = rbind(X_unit, xnext_unit)
      xnext = unnormalize(xnext_unit, B)
      X = rbind(X, xnext)
      
      ## new run
      out = blackbox(xnext, ...)
      fnext = out$obj; obj = c(obj, fnext)
      C = rbind(C, bilog(out$c))
      
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
        cat("k=", k+cl, " ", sep="")
        cat("; xnext ([", paste(signif(xnext,3), collapse=" "), 
            "], feasibility=", feasibility[k], ")\n", sep="")
        cat(" xbest=[", paste(signif(xbest,3), collapse=" "), sep="")
        cat("]; ybest (prog=", paste(signif(m2,3), collapse=" ")) 
        cat(", since=", since, ")\n", sep="")
      }
    }
    
    ## plot Pareto Front and Pareto set
    if (plotPareto) {
      par(ps=16, mfrow=c(1,2))
      # Pareto Front
      plot(AF_PF, ylim = c(-1, 0), 
           pch = 1 + AF_cl$cluster, col = 1 + AF_cl$cluster,
           xlab="-logEI", ylab="-PoF", main="Objective space")
      # Pareto Set
      plot(unnormalize(AF_PS[,1:2], B), 
           xlim = B[1,], ylim = B[2,], cex = 0.6,
           xlab = "x1", ylab = "x2", main="Parameter space")
      points(tail(X[,1:2], nprl), 
             col = 1 + (1:nprl), pch = 18, cex = 1.6)
      par(mfrow=c(1,1))
    }
    
    ## update GP fits
    updateGPsep(fgpi, tail(X_unit, nprl), (tail(obj,nprl)-fmean)/fsd, verb = 0)
    df = mleGPsep(fgpi, param = "d", tmin = dlim[1], tmax = dlim[2], ab = ab, verb = 0)$d
    for(j in 1:nc){
      updateGPsep(Cgpi[j], tail(X_unit, nprl), tail(C[,j],nprl), verb = 0)
      dc[j,] = mleGPsep(Cgpi[j], param = "d",  tmin = dlim[1], tmax = dlim[2], ab = ab, verb = 0)$d
    }
  }
  
  ## delete GP surrogates
  deleteGPsep(fgpi)
  for(j in 1:nc) deleteGPsep(Cgpi[j])
  
  return(list(prog = prog, xbest = xbest, obj = obj, C=C, X = X))
}
