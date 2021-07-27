# Internal function
#
# @param dt data.frame of time variables
# @param X matrix or NULL, covariates
# @param changePts numeric change point(s) of log hazard ratio
# @param constantVE TRUE: constant VE; FALSE: waning VE
# @param timePts numeric upper bound of intervals
# @param rightCens TRUE: standard Cox; FALSE: proposed method
# @param plots T/F, whether to plot curves of VE_a and VE_h or not
# @param threshold numeric, threshold of the convergence criterion
# @param maxit integer, maximum number of iterations

#' @include postProcess.R 
#' @importFrom stats pnorm

CoxReg <- function(dt,  
                   X,  
                   changePts,
                   constantVE,  
                   timePts,  
                   plots,  
                   threshold,  
                   maxit) {

  n <- nrow(x = dt)

  if (is.null(x = X)) {
    p <- 0L
    varname <- NULL
    SD <- NULL
    X <- matrix(data = 0.0, nrow = n, ncol = 1L)
  } else {
    p <- ncol(x = X)
    varname <- colnames(x = X)

    # standardize covariates X
    SD <- apply(X = X, MARGIN = 2L, FUN = sd, na.rm = TRUE)
    X <- scale(x = X, center = FALSE, scale = SD)
  }
  
  # knots (in days) of linear splines 
  if (is.null(x = changePts)) {
    knots <- as.list((4L:8L)*7L)
  } else {
    knots <- list(changePts)
  }

  if (constantVE) {
    npc <- length(x = knots[[1L]])
  } else {
    npc <- length(x = knots[[1L]]) + 1L
  }
  
  message("performing the standard Cox regression") 

  YD <- dt[,"right_time"]
  delta <- rep(x = 1L, times = n) # time to first detection YD and status delta
  YD[dt[,"right_time"]<0L] <- dt[dt[,"right_time"]<0L,"left_time"]
  delta[dt[,"right_time"]<0L] <- 0L
  tau <- max(YD) + 1L
    
  # sort the data by YD
  sorted_index <- order(YD)
  YD <- YD[sorted_index]
  delta <- delta[sorted_index]
  dt <- dt[sorted_index,]
  X <- X[sorted_index,,drop=FALSE]
    
  args <- list("beta"= rep(x = 0.0, times = p),
               "gamma" = rep(x = 0.0, times = npc), 
               "time" = YD,
               "delta" = delta,
               "X" = X,
               "W" = dt[,"entry_time"],
               "S" = dt[,"vaccination_time"], 
               "constantVE" = constantVE, 
               "threshold" = threshold, 
               "maxit" = maxit)

  if (p == 0L) {
    funcCox <- "Cox_noX"
  } else {
    funcCox <- "Cox"
  }

  # fit the standard Cox model

  if (length(x = knots) == 1L) {

    args[[ "knots" ]] <- knots[[ 1L ]]

    fit.cox <- tryCatch(expr = do.call(what = funcCox, args = args),
                        error = function(e){message(e$message); return( NULL )})

    if (is.null(x = fit.cox)) stop("calculation aborted", call. = FALSE)

    if (p != 0L) {
      theta <- c(args[[ "beta" ]], args[[ "gamma" ]])
    } else {
      theta <- args[[ "gamma" ]]
    }

    ll0 <- fit.cox[[ 1L ]]
    covar <- fit.cox[[ 2L ]]
    finalknot <- 1L

    if (!anyNA(x = theta)) {
      message("Number of subjects: ", n) 
      message("Partial log-likelihood at final estimates: ", ll0)
    } else {
      stop("NA values were produced in the standard Cox regression")
    }

  } else {

    curll <- -Inf
    finalknot <- NA

    for (i in 1L:length(x = knots)) {


      args[[ "beta" ]] <- rep(x = 0.0, times = p)
      args[[ "gamma" ]] <- rep(x = 0.0, times = npc)

      args[[ "knots" ]] <- knots[[ i ]]

      fit.cox <- tryCatch(expr = do.call(what = funcCox, args = args),
                          error = function(e){message(e$message); return( NULL )})

      if (is.null(x = fit.cox)) next

      if (p != 0L) {
        temptheta <- c(args[[ "beta" ]], args[[ "gamma" ]])
      } else {
        temptheta <- args[[ "gamma" ]]
      }

      ll <- fit.cox[[ 1L ]]
      
      if (!anyNA(x = temptheta) & ll>curll) {
        finalknot <- i
        curll <- ll
        theta <- temptheta
        covar <- fit.cox[[ 2L ]]
      }
    }
    
    if (!is.na(x = finalknot)) {
      message("Day ", knots[[ finalknot ]], " (week ", knots[[ finalknot ]]/7, 
              ") was selected as the change point by AIC")
      message("Number of subjects: ", n) 
      message("Partial log-likelihood at final estimates: ", curll)
      
    } else {
      stop("all candidate change points produced NA values. ",
           "No change point was selected by AIC", call. = FALSE)
    }
  }

  knots <- knots[[ finalknot ]]

  res <- .postProcess(theta = theta, 
                      covMat = covar, 
                      plots = plots, 
                      SD = SD, 
                      varname = varname, 
                      constantVE = constantVE, 
                      tau = tau,
                      knots = knots,
                      timePts = timePts)

  res[[ "knots" ]] <- knots
  res[[ "gamma" ]] <- theta[{p+1L}:{p+npc}]
  res[[ "covgamma" ]] <- as.matrix(covar[{p+1L}:{p+npc}, {p+1L}:{p+npc}])
  res[[ "tau" ]] <- tau

  
  return( res )
}
