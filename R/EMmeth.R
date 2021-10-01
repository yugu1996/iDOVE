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
    SD <- numeric(length = 0)
    X <- matrix(data = 0.0, nrow = n, ncol = 1L)
  } else {
    p <- ncol(x = X)

    varname <- colnames(x = X)

    # standardize covariates X
    SD <- apply(X = X, MARGIN = 2L, FUN = sd, na.rm = TRUE)
    X <- scale(x = X, center = FALSE, scale = SD)
  }

  npc <- length(x = changePts[[ 1L ]]) + 1L - constantVE
  
  message("performing nonparametric maximum likelihood") 

  # Rstar resets cases without right boundary to use left boundary
  Rstar <- dt[,"right_time"]
  Rstar[dt[,"right_time"] < 0L] <- dt[dt[,"right_time"] < 0L, "left_time"]
    
  # sort the data by Rstar
  sorted_index <- order(Rstar)
  dt <- dt[sorted_index,]
  X <- X[sorted_index,,drop=FALSE]
  Rstar <- Rstar[sorted_index]

  # all unique interval times
  times <- sort(x = unique(x = c(dt[,"left_time"], Rstar)))
  if (times[1L] == 0L) times <- times[-1L]

  m <- length(x = times)
  tau <- times[m]

  # initialize EM algorith
  # Returns {n x 3} matrix, containing the element of t
  # nearest the entry, left interval, and right interval times
  tind <- Init(dt[,"left_time"], Rstar, dt[,"entry_time"], times)

  holdw <- .resetW(m = m, 
                   left_time = dt[,"left_time"], 
                   right_time = dt[,"right_time"],
                   tind = tind)

  args <- list("w" = holdw, 
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
    # functions and arguments when no covariates are in model
    funcEM <- "EM_noX"
    funcLL <- "LL_noX"
    funcPLL <- "PLL_noX"
    EMargs <- c("w", "gamma", "lambda", "L", "R",
                "Rstar", "t", "tind", "S", "knots",
                "constantVE", "threshold", "maxit")
    LLargs <- c("gamma", "lambda", "L", "R", "t",
                "tind", "S", "knots", "constantVE")
    PLLargs <- c("w", "gamma", "lambda_init", "L", "R",
                 "Rstar", "t", "tind", "S", "knots",
                 "constantVE", "threshold", "maxit")
  } else {
    # functions and arguments when covariates are in model
    funcEM <- "EM"
    funcLL <- "LL"
    funcPLL <- "PLL"
    EMargs <- c("w", "beta", "gamma", "lambda", "L", 
                "R", "Rstar", "t", "tind", "X", 
                "S", "knots", "constantVE", "threshold", "maxit")
    LLargs <- c("beta", "gamma", "lambda", "L", "R", 
                "t", "tind", "X", "S", "knots", 
                "constantVE")
    PLLargs <- c("w", "beta", "gamma", "lambda_init", "L", 
                 "R", "Rstar", "t", "tind", "X", 
                 "S", "knots", "constantVE", "threshold", "maxit")
  }

  if (length(x = changePts) == 1L) {

    res <- .noaic(changePts = changePts, 
                  args = args, 
                  p = p, 
                  funcEM = funcEM, 
                  EMargs = EMargs, 
                  funcLL = funcLL, 
                  LLargs = LLargs,
                  n = n,
                  m = m)

  } else {

    res <- .aicSearch(changePts = changePts, 
                      args = args, 
                      p = p, 
                      npc = npc, 
                      m = m, 
                      holdw = holdw, 
                      funcEM = funcEM, 
                      EMargs = EMargs,
                      funcLL = funcLL, 
                      LLargs = LLargs,
                      n = n)


  }

  args[[ "knots" ]] <- changePts[[ res$finalknot ]]
  args[[ "lambda_init" ]] <- res$lambda

  covar <- .covStep(n = n, p = p, npc = npc, 
                    theta = res$theta, 
                    args = args, 
                    funcPLL = funcPLL, 
                    PLLargs = PLLargs, 
                    ll0 = res$ll0)
      
  post <- .postProcess(theta = res$theta, 
                       covMat = covar, 
                       plots = plots, 
                       SD = SD, 
                       varname = varname, 
                       constantVE = constantVE, 
                       tau = tau,
                       changePts = changePts[[ res$finalknot ]],
                       timePts = timePts)

  post$changePts <- changePts[[ res$finalknot ]]

  return( post )
}

.noaic <- function(changePts, args, p, funcEM, EMargs, funcLL, LLargs, n, m) {

  args[[ "knots" ]] <- changePts[[ 1L ]]

  tryCatch(expr = do.call(what = funcEM, args = args[ EMargs ]),
           error = function(e){ stop(e$message, call. = FALSE) })

  if (p != 0L) {
    theta <- c(args[[ "beta" ]], args[[ "gamma" ]])
  } else {
    theta <- args[[ "gamma" ]]
  }

  ll0 <- tryCatch(expr = do.call(what = funcLL, args = args[ LLargs ]),
                  error = function(e){ stop(e$message, call. = FALSE) })

  message("Number of subjects: ", n) 
  message("Number of unique time points: ", m)
  message("Log-likelihood at final estimates: ", 
          ifelse(test = abs(x = sum(ll0)) > 1.0, 
                 yes = round(x = sum(ll0), digits = 2L),
                 no = round(x = sum(ll0), digits = 6L)))

  return( list("finalknot" = 1L, 
               "lambda" = args[[ "lambda" ]],
               "theta" = theta,
               "ll0" = ll0) )
}

.aicSearch <- function(changePts, args, p, npc, m, holdw, funcEM, EMargs,
                       funcLL, LLargs, n) {

  curll <- -Inf
  finalknot <- NA

  for (i in 1L:length(x = changePts)) {

    message("evaluating change point: ", changePts[[i]])

    args[[ "beta" ]] <- rep(x = 0.0, times = p)
    args[[ "gamma" ]] <- rep(x = 0.0, times = npc)
    args[[ "lambda" ]] <- rep(x = 1.0/m, times = m)
    args[[ "w" ]] <- holdw

    args[[ "knots" ]] <- changePts[[ i ]]

    tst <- tryCatch(expr = do.call(what = funcEM, args = args[ EMargs ]),
                    error = function(e){ message(e$message); return( NA ) })

    if (!is.null(x = tst)) next

    if (p != 0L) {
      temptheta <- c(args[[ "beta" ]], args[[ "gamma" ]])
    } else {
      temptheta <- args[[ "gamma" ]]
    }

    if (!anyNA(x = temptheta)) {

      ll <- tryCatch(expr = do.call(what = funcLL, args = args[ LLargs ]),
                     error = function(e){ stop(e$message, call. = FALSE) })

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
    message("Day ", changePts[[ finalknot ]], 
            " (week ", changePts[[ finalknot ]]/7, 
            ") was selected as the change point by AIC")
    message("Number of subjects: ", n) 
    message("Number of unique time points: ", m)
    message("Log-likelihood at final estimates: ", 
            ifelse(test = abs(x = sum(curll)) > 1.0, 
                   yes = round(x = sum(curll), digits = 2L),
                   no = round(x = sum(curll), digits = 6L)))
  } else {
    stop("all candidate change points produced NA values; ",
         "no change point was selected by AIC", call. = FALSE)
  }

  return( list("finalknot" = finalknot, 
               "lambda" = templambda,
               "theta" = theta,
               "ll0" = ll0) )
}

# initializes w matrix to be one at t_m nearest right_time for individuals
# with exactly observed events and zero otherwise.
.resetW <- function(m, left_time, right_time, tind) {

  w1 <- matrix(data = 0.0, nrow = length(x = left_time), ncol = m)

  tst <- left_time == right_time

  w1[cbind(which(tst),tind[tst,2]+1)] <- 1.0

  return( w1 )
}

.covStep <- function(n, p, npc, theta, args, funcPLL, PLLargs, ll0) {

  # use the first-order derivative to estimate covariance
  hn <- 5/sqrt(n)
  pll <- matrix(data = 0.0, nrow = n, ncol = p+npc)
      
  for (j in 1L:(p+npc)) {

    theta1 <- theta
    theta1[j] <- theta1[j] + hn

    if (p != 0L) args[[ "beta" ]] <- theta1[1L:p]
    args[[ "gamma" ]] <- theta1[{p+1L}:{p+npc}]

    pll[,j] <- tryCatch(expr = do.call(what = funcPLL, args = args[ PLLargs ]),
                       error = function(e){ stop(e$message, call. = FALSE) })

  }

  numdif <- {pll - matrix(data = ll0, nrow = n, ncol = p+npc, byrow = FALSE)}/hn
  covar_inv <- crossprod(x = numdif)

  if (is.finite(x = determinant(x = covar_inv)$modulus)) {
    covar <- solve(a = covar_inv)
  } else {
    stop("covariance matrix estimator is singular", call. = FALSE)
  }

  return( covar )
}
