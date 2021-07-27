# Internal function
#
# @param dt data.frame containing all time variables
# @param X matrix or NULL, covariates
# @param changePts numeric change point(s) of log hazard ratio
# @param constantVE TRUE: constant VE; FALSE: waning VE
# @param timePts numeric upper bound of intervals
# @param rightCens TRUE: standard Cox; FALSE: proposed method
# @param plots T/F, whether to plot curves of VE_a and VE_h or not
# @param threshold numeric, threshold of the convergence criterion
# @param maxit integer, maximum number of iterations

#' @include postProcess.R 
#' @importFrom stats pnorm sd

EMmeth <- function(dt,  
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

    # these simplify later logic
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
    npc <- length(x = knots[[ 1L ]])
  } else {
    npc <- length(x = knots[[ 1L ]]) + 1L
  }
  

  # apply the propose method for interval censored data
  # check if the data need to be modified due to no right-censored subjects
  if (max(dt[,"left_time"]) < max(dt[,"right_time"])) {
    dt[dt[,"right_time"] > max(dt[,"left_time"]), "right_time"] <- -1L
    message("all subjects whose right-interval time is greater than ",
            "the biggest left-interval time made to be right-censored.")
  }
    
  message("performing nonparametric maximum likelihood") 

  Rstar <- dt[,"right_time"]
  Rstar[dt[,"right_time"] < 0L] <- dt[dt[,"right_time"] < 0L, "left_time"]
    
  # sort the data by Rstar
  sorted_index <- order(Rstar)
  dt <- dt[sorted_index,]
  X <- X[sorted_index,,drop=FALSE]
  Rstar <- Rstar[sorted_index]

  times <- sort(x = unique(x = c(dt[,"left_time"], Rstar)))
  if (times[1L] == 0L) times <- times[-1L]

  m <- length(x = times)
  tau <- times[m] + 1L

  # initialize EM algorithm
  tind <- Init(dt[,"left_time"], Rstar, dt[,"entry_time"], times)

  resetW <- function() {
    w1 <- matrix(data = 0.0, nrow = n, ncol = m)
    tst <- dt[,"left_time"] == dt[,"right_time"]
    w1[cbind(which(tst),tind[tst,2]+1)] <- 1.0
    return( w1 )
  }

  args <- list("w" = resetW(), 
               "beta"= rep(x = 0.0, times = p),
               "gamma" = rep(x = 0.0, times = npc), 
               "lambda" = rep(x = 1.0/m, times = m), 
               "L" = dt[,"left_time"], 
               "R" = dt[,"right_time"], 
               "Rstar" = Rstar, 
               "t" = times, 
               "tind" = tind, 
               "X" = X,
               "S" = dt[,"vaccination_time"], 
               "constantVE" = constantVE, 
               "threshold" = threshold, 
               "maxit" = maxit)

  if (p == 0L) {
    funcEM <- "EM_noX"
    funcLL <- "LL_noX"
    funcPLL <- "PLL_noX"
  } else {
    funcEM <- "EM"
    funcLL <- "LL"
    funcPLL <- "PLL"
  }

  if (length(x = knots) == 1L) {

    args[[ "knots" ]] <- knots[[ 1L ]]

    tryCatch(expr = do.call(what = funcEM, args = args),
             error = function(e){stop(e$message, call. = FALSE)})

    if (p != 0L) {
      theta <- c(args[[ "beta" ]], args[[ "gamma" ]])
    } else {
      theta <- args[[ "gamma" ]]
    }

    lambda <- args[[ "lambda" ]]

    ll0 <- tryCatch(expr = do.call(what = funcLL, args = args),
                    error = function(e){stop(e$message, call. = FALSE)})

    finalknot <- 1L

    message("Number of subjects: ", n) 
    message("Number of unique time points: ", m)
    message("Log-likelihood at final estimates: ", sum(ll0))

  } else {

    curll <- -Inf
    finalknot <- NA

    for (i in 1L:length(x = knots)) {

      args[[ "beta" ]] <- rep(x = 0.0, times = p)
      args[[ "gamma" ]] <- rep(x = 0.0, times = npc)
      args[[ "lambda" ]] <- rep(x = 1.0/m, times = m)
      args[[ "w" ]] <- resetW()

      args[[ "knots" ]] <- knots[[ i ]]

      tst <- tryCatch(expr = do.call(what = funcEM, args = args),
                      error = function(e){message(e$message); return( NA )})

      if (!is.null(x = tst)) next

      if (p != 0L) {
        temptheta <- c(args[[ "beta" ]], args[[ "gamma" ]])
      } else {
        temptheta <- args[[ "gamma" ]]
      }

      if (!anyNA(x = temptheta)) {
        ll <- tryCatch(expr = do.call(what = funcLL, args = args),
                       error = function(e){stop(e$message, call. = FALSE)})
        templl = sum(ll)

        if(templl > curll) {
          finalknot <- i
          curll <- templl
          theta <- temptheta
          templambda <- args[[ "lambda" ]]
          ll0 <- ll
        }
      }
      
    }
    
    if (!is.na(x = finalknot)) {
      message("Day ", knots[[ finalknot ]], " (week ", knots[[ finalknot ]]/7, 
              ") was selected as the change point by AIC")
      message("Number of subjects: ", n) 
      message("Partial log-likelihood at final estimates: ", curll)
      lambda <- templambda
    } else {
      stop("all candidate change points produced NA values. ",
           "No change point was selected by AIC", call. = FALSE)
    }
  }

  args[[ "knots" ]] <- knots[[ finalknot ]]
  args[[ "lambda" ]] <- lambda
      
  # use the first-order derivative to estimate covariance
  hn <- 5/sqrt(n)
  pll <- matrix(data = 0.0, nrow = n, ncol = p+npc)
      
  for (j in 1L:(p+npc)) {
    theta1 <- theta
    theta1[j] <- theta1[j] + hn

    if (p != 0L) args[[ "beta" ]] <- theta1[1L:p]
    args[[ "gamma" ]] <- theta1[{p+1L}:{p+npc}]

    pll[,j] <- tryCatch(expr = do.call(what = funcPLL, args = args),
                       error = function(e){stop(e$message, call. = FALSE)})

  }

  numdif <- {pll - matrix(data = ll0, nrow = n, ncol = p+npc, byrow = FALSE)}/hn
  covar_inv <- crossprod(x = numdif)
  if (is.finite(x = determinant(x = covar_inv)$modulus)) {
    covar <- solve(a = covar_inv)
  } else {
    stop("covariance matrix estimator is singular", call. = FALSE)
  }

  res <- .postProcess(theta = theta, 
                      covMat = covar, 
                      plots = plots, 
                      SD = SD, 
                      varname = varname, 
                      constantVE = constantVE, 
                      tau = tau,
                      knots = knots[[ finalknot ]],
                      timePts = timePts)

  res[[ "knots" ]] <- knots[[ finalknot ]]
  res[[ "gamma" ]] <- theta[{p+1L}:{p+npc}]
  res[[ "covgamma" ]] <- as.matrix(covar[{p+1L}:{p+npc}, {p+1L}:{p+npc}])
  res[[ "tau" ]] <- tau
  
  return( res )
}
