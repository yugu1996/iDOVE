#' Durability of Vaccine Efficacy Against Asymptomatic SARS-CoV-2 Infection
#' 
#' Assesses potentially time-varying vaccine efficacy (VE) against SARS-CoV-2 
#' infection under staggered enrollment and time-varying community transmission, 
#' allowing crossover of placebo volunteers to the vaccine arm. The infection 
#' time data can be either interval- or right-censored, with the latter being 
#' a special case of the former.
#' 
#' The information required for an analysis is 
#'   \describe{
#'     \item{Entry Time:}{The time when the participant enters the trial 
#'        in whole units days.}
#'     \item{Left Interval Time:}{The last examination time when the test is 
#'        negative in whole units days.}
#'     \item{Right Interval Time:}{The first examination time when the test 
#'        is positive in whole units days. If the participant does not 
#'        test positive during the trial, use NA or Inf.}
#'     \item{Vaccination Time:}{The time when vaccination takes place in 
#'        whole units days. If the participant is not vaccinated during 
#'        the trial, use NA or Inf.}
#'     \item{Covariates:}{Baseline covariates (e.g., priority group, age, 
#'        ethnicity).}
#'    }
#'    
#' The covariates can include categorical variables, for which
#' all other categories are compared to the first category.
#'
#' Note that all the time variables are measured from the start of the 
#' clinical trial and are specified in whole units of days. Though they
#' need not be provided as integer, all non-NA and finite values must be able 
#' to be cast as integers without loss of information. For each 
#' individual, the entry_time, left_time, and right_time satisfy
#' entry_time \eqn{\le} left_time \eqn{\le} right_time. 
#' 
#' The general structure of the formula input is
#'   \preformatted{
#'   intCens(entry_time, left_time, right_time, vaccination_time) ~ covariates 
#'   }
#' 
#' The response variable contains all of the time information. It must be
#' specified through function 'intCens()'. Specifically, 
#' \preformatted{intCens(entry_time, left_time, right_time, vaccination_time)}
#' If entry_time > left_time, or left_time > right_time, the case will be 
#' removed from the analysis and a message will be generated.  
#' 
#' When right_time - left_time \eqn{\le} 2 for all individuals whose infection 
#' times are truly interval-censored (i.e., right_time is finite and less than
#' the end of the trial), idove() assumes that examinations are completed
#' daily or every two days and performs a standard Cox regression, regardless
#' of the value specified for input rightCens. 
#' 
#' The log hazard ratio with respect to vaccination, denoted by \eqn{\eta(t)}, 
#' is approximated by linear splines. Specifically, in \eqn{K}-piece linear 
#' splines,
#' \deqn{\eta(t) = \gamma_1t+\gamma_2(t-t_1)_+ + \gamma_3(t-t_2)_+ +\dots + \gamma_K(t-t_{K-1})_+,}
#' where \eqn{t_1, \dots, t_{K-1}} are the \eqn{K-1} pre-specified change points 
#' and \eqn{\gamma_1, \dots, \gamma_K} are the \eqn{K} spline parameters. 
#' The first \eqn{K-1} spline parameters are always estimated from the data, 
#' whereas the treatment of the last parameter \eqn{\gamma_K} depends on the 
#' value of input constantVE.
#' If constantVE = TRUE, the slope of the last piece is assumed to be zero, 
#' that is, \eqn{\gamma_K = -\sum_{k=1}^{K-1}\gamma_k}. Thus, \eqn{\gamma_K} is 
#' a redundant parameter and is not estimated. If constantVE = FALSE, the slope 
#' of the last piece can be nonzero and \eqn{\gamma_K} is also estimated from 
#' the data. 
#' 
#' @rdname idove
#' @name idove 
#' 
#' @references Lin, D-Y, Gu, Y., Zeng, D., Janes, H. E., and Gilbert, P. B. (2021). 
#'   Evaluating Vaccine Efficacy Against SARS-CoV-2 Infection. Submitted.
#' 
#' @param formula A formula object, with the response (all the time information)
#'   on the left hand side of a '~' operator and the covariates on the right. 
#'   The response must be specified through the intCens() function. 
#'   See ?intCens and Details for further information. 
#'   
#' @param data A data.frame object. The data.frame in which to interpret the
#'   variable names in formula. Must contain the entry time, the left interval,
#'   time, the right interval time, the vaccination time, and the covariates. 
#'   See Details.
#'   
#' @param constantVE A logical object. If TRUE, VE is assumed to be constant. 
#'   If FALSE (default), VE is assumed to be waning. Note that constant versus 
#'   waning VE pertains only to the period after ramping VE.
#'   
#' @param rightCens A logical object. If TRUE, the standard Cox model will be 
#'   fitted, by treating the time of first positive test as a potentially 
#'   right-censored event time and performing maximum partial likelihood 
#'   estimation. If FALSE (default), the interval-censored event time method 
#'   will be applied. When examinations are completed daily or every two days,
#'   this input is ignored and a standard Cox regression is performed. 
#'   See Details.  
#' 
#' @param plots A logical object. If TRUE (default), plots of the estimated 
#'   VE in reducing attack rate, the estimated VE in reducing the hazard rate, 
#'   and their 95\% confidence intervals will be automatically generated. 
#'   If FALSE, plots will not be generated.
#'   
#' @param changePts An integer vector object or NULL. The potential change
#'   points (in days) of the piece-wise log-linear hazard ratio. See Details 
#'   for further information. If NULL, AIC will be used to select one change 
#'   point from \{28, 35, 42, 49, 56\} (weeks 4, 5, 6, 7, 8).
#'   
#' @param timePts An integer vector object or NULL. The endpoints (in days) 
#'   of the time periods for which the VE in reducing the attack rate are to
#'   be estimated. The estimated VE in reducing the hazard rate at 
#'   these endpoints are also returned. If NULL, a default sequence 
#'   \eqn{t_1, 2t_1, 3t_1, \dots} will be used, where \eqn{t_1} be the first 
#'   change point. The sequence ends at the maximum of the left and
#'   and right ends of the time intervals from all participants. This input is 
#'   ignored when constantVE = TRUE.  
#'   
#' @param tol A numeric scalar. The convergence threshold for the EM or 
#'   Newton-Raphson algorithm. The default value is 0.0001. 
#'   
#' @param maxit A positive integer. The maximum number of iterations for the 
#'   EM or Newton-Raphson algorithm. The default value is 2000. 
#'   
#' @returns An S3 object of class iDOVE containing a list with elements
#'
#'   \item{call}{The unevaluated call.}
#'
#'   \item{covariates}{A matrix containing the estimated (log) hazard ratio of
#'     each covariate, together with the estimated standard error, the 95\%
#'     confidence intervals, and the two-sided p-value for testing no covariates
#'     effect. NA if only an intercept is given as the right hand side in 
#'     input formula.}
#'     
#'   \item{vaccine}{A list containing one or three elements, depending on the 
#'     value of constantVE. 
#'     If constantVE = TRUE, the only element is named 'VE' and is a vector 
#'     containing the estimate of constant VE, its standard error estimate, 
#'     and the 95\% confidence interval.
#'     If constantVE = FALSE, three matrices are returned. The first matrix 
#'     named 'VE_a' contains the estimates of the VE in reducing the attack 
#'     rate at all time points given in timePts, together with the 95\%
#'     confidence intervals. The second matrix named 'VE_h' contains the 
#'     estimates of the VE in reducing the hazard rate at timePts. 
#'     The third matrix named 'VE_period' contains the estimates of VE in
#'     reducing the attack rate over successive time periods according to
#'     timePts, together with the 95\% confidence intervals.}
#'     
#'   Objects of class iDOVE have additional attributes, knots, tau, gamma, and
#'     covgamma, which are included for plotting capabilities
#'
#' @export
#' @import methods
#'
#' @useDynLib iDOVE
#' 
#' @include CoxReg.R EMmeth.R 
#' @import Rcpp
#' @importFrom graphics axis legend lines plot plot.new
#' @importFrom stats model.response update.formula complete.cases terms
#' 
#' @examples
#' data(idoveData)
#'
#' set.seed(1234)
#' smp <- sample(1L:nrow(x = idoveData), size = 500L)
#' 
#' # Fit the model with default settings
#' idove(formula = intCens(entry.time, left.time, right.time, vaccine.time) ~ 1, 
#'       data = idoveData[smp,])
#' 
#' # Specify Week 4 as the change point
#' # Assume a potentially waning VE after 4 weeks
#' # Estimate VE_a over 0-4, 4-16, 16-28, 28-40 weeks
#' idove(formula = intCens(entry.time, left.time, right.time, vaccine.time) ~ 1, 
#'       data = idoveData[smp,],
#'       changePts = 4*7,
#'       timePts = c(4, 16, 28, 40)*7)
#'       
#' # Specify multiple change points at Weeks 4 and 8
#' # Assume a constant VE after 8 weeks
#' idove(formula = intCens(entry.time, left.time, right.time, vaccine.time) ~ 1, 
#'       data = idoveData[smp,],
#'       changePts = c(4, 8)*7,
#'       constantVE = TRUE)

idove <- function(formula,  
                  data, 
                  constantVE = FALSE,  
                  rightCens = FALSE,  
                  plots = TRUE, 
                  changePts = NULL,  
                  timePts = NULL, 
                  tol = 0.0001,  
                  maxit = 2000) {

  # matched call to include with returned object
  cl <- match.call()
  
  # formula and data must be provided

  if (missing(x = formula)) {
    stop("a formula argument must be provided", call. = FALSE)
  }
  
  if (missing(x = data)) {
    stop("a data argument must be provided", call. = FALSE)
  }
  
  # add intercept from model if provided
  # this was added 6/17/21 to ensure that factors are handled properly
  if (attr(x = stats::terms(x = formula), which = "intercept") == 0L) {
    formula = update.formula(old = formula, new = .~. +1)
  }
  
  # reset options to allow for keeping NA values
  opt <- options()
  options(na.action = 'na.pass')
  on.exit(options(opt))
  
  # try to obtain the model.frame
  mf <- tryCatch(expr = stats::model.frame(formula = formula, data = data),
                 error = function(e) {
                   message("unable to obtain model.frame")
                   stop(e$message, call. = FALSE)
                 })
  
  # extract covariates (X is a matrix)
  X <- suppressMessages(stats::model.matrix(object = mf, data = data))
  # 6/17/21 remove intercept
  int <- attr(x = X, which = "assign") != 0L
  X <- X[,int, drop = FALSE]
  
  if (ncol(x = X) == 0L) X <- NULL
  
  # extract time variables
  dt <- suppressMessages(stats::model.response(data = mf))
  
  if(!all.equal(target = colnames(x = dt), 
                current = c("entry_time", 
                            "left_time", 
                            "right_time", 
                            "vaccination_time"))) {
    stop("the LHS of formula did not contain an appropriate intCens() object",
         call. = FALSE)
  }
  
  # remove any cases that have NA in dt or X
  use <- stats::complete.cases(cbind(X, dt))
  
  if (all(!use, na.rm = TRUE)) {
    stop("input checks result in all NA -- verify inputs", call. = FALSE)
  }

  if (any(!use, na.rm = TRUE)) {
    dt <- dt[use,]
    if (!is.null(x = X)) X <- X[use,]
  } 
  
  if (sum(!use, na.rm = TRUE) > 0L) {
    message(sum(!use), 
            " cases removed from the analysis due to NA values")
  }
  
  # verify changePts
  if (!is.null(x = changePts)) {

    if (!is.numeric(x = changePts)) {
      stop("changePts must be a numeric vector or NULL", call. = FALSE)
    }
    
    if (length(x = changePts) == 0L) {
      message("changePts has zero length; default values will be used")
      changePts <- NULL
    } else if (any(changePts < 0.0)) {
      stop("all changePts must be positive", call. = FALSE)
    }

    if (!is.integer(x = changePts)) {
      tmp <- as.integer(x = round(x = changePts, digits = 0L))
      if (!isTRUE(x = all.equal(target = changePts, current = tmp))) {
        stop("all changePts must be integer (days)", call. = FALSE)
      }
      changePts <- tmp
    }
  } else {
    message("changePts not given; using AIC to select from {28, 35, 42, 49, 56}")
  }
  
  # verify timePts
  if (!is.null(x = timePts)) {
    if (!is.numeric(x = timePts)) {
      stop("timePts must be a numeric vector or NULL", call. = FALSE)
    }
    
    if (length(x = timePts) == 0L) {
      message("timePts has zero length; default values will be used")
      timePts <- NULL
    } else if (any(timePts < 0.0)) {
      stop("all timePts must be positive", call. = FALSE)
    }
    if (!is.integer(x = timePts)) {
      tmp <- as.integer(x = round(x = timePts, digits = 0L))
      if (!isTRUE(x = all.equal(target = timePts, current = tmp))) {
        stop("all timePts must be integer (days)", call. = FALSE)
      }
      timePts <- tmp
    }
  } else if (!constantVE) {
    message("timePts not given; default values will be used")
  }
  
  # create separate time variables (E,L,R,S)

  # right_time is set negative if it is provided as Inf in input
  # do not consider those in this test
  ind <- dt[,"right_time"] >= 0L
  is.rightCens <- prod(dt[ind,"right_time"]-dt[ind,"left_time"] <= 2L) | rightCens
  
  # apply the right-censored version of the proposed method (standard Cox)
  if (is.rightCens && !rightCens) {
    message("gap between two successive examinations is <= 2 days")
    rightCens <- TRUE
  }

  ### main analysis
  
  if (!rightCens) {

    res <- EMmeth(dt = dt,  
                  X = X,  
                  changePts = changePts,  
                  constantVE = constantVE,  
                  timePts = timePts,  
                  plots = plots,  
                  threshold = tol,  
                  maxit = maxit)
  } else {
  
    res <- CoxReg(dt = dt,  
                  X = X,  
                  changePts = changePts,  
                  constantVE = constantVE,  
                  timePts = timePts,  
                  plots = plots,  
                  threshold = tol,  
                  maxit = maxit)
  }

  res[["call"]] <- cl

  class(x = res) <- "iDOVE"
  attr(x = res, which = "knots") <- res[[ "knots" ]]
  attr(x = res, which = "gamma") <- res[[ "gamma" ]]
  attr(x = res, which = "covgamma") <- res[[ "covgamma" ]]
  attr(x = res, which = "tau") <- res[[ "tau" ]]

  res[[ "knots" ]] <- NULL
  res[[ "gamma" ]] <- NULL
  res[[ "covgamma" ]] <- NULL
  res[[ "tau" ]] <- NULL
  
  return( res )
}
