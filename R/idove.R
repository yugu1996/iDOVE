#' Durability of Vaccine Efficacy Against Asymptomatic SARS-CoV-2 Infection
#' 
#' Assesses potentially time-varying vaccine efficacy (VE) against SARS-CoV-2 
#' infection under staggered enrollment and time-varying community transmission, 
#' allowing crossover of placebo volunteers to the vaccine arm. The infection 
#' time data are interval-censored, and the log hazard ratio is assumed to be 
#' a piece-wise linear function of time. 
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
#' A model without covariates is also allowed.
#'
#' Note that all of the time variables are measured from the start of the 
#' clinical trial and are specified in whole units of days. Though they
#' need not be provided as integer, all non-NA and finite values must be able 
#' to be cast as integers without loss of information. For each 
#' individual, the entry_time and left_time should satisfy
#' entry_time \eqn{\le} left_time. For each individual that tests positive, 
#' entry_time \eqn{\le} left_time \eqn{\le} right_time. For each individual 
#' that is vaccinated, entry_time \eqn{\le} vaccination_time.
#' 
#' The general structure of the formula input is
#'   \preformatted{
#'   intCens(entry_time, left_time, right_time, vaccination_time) ~ covariates 
#'   }
#' 
#' The left-hand side contains all of the time information. It must be
#' specified through function 'intCens()'. Specifically, 
#' \preformatted{intCens(entry_time, left_time, right_time, vaccination_time)}
#' If entry_time > left_time, or left_time > right_time, the case will be 
#' removed from the analysis and a message will be generated. 
#' 
#' The special case of right-censored data is implemented by dove2() in the 
#'   DOVE package available through CRAN.
#' 
#' @rdname idove
#' @name idove 
#' 
#' @references Lin, D-Y, Gu, Y., Zeng, D., Janes, H. E., and Gilbert, P. B. (2021). 
#'   Evaluating Vaccine Efficacy Against SARS-CoV-2 Infection.
#'   Clinical Infectious Diseases, ciab630, https://doi.org/10.1093/cid/ciab630
#' 
#' @param formula A formula object, with all of the time variables
#'   on the left hand side of a '~' operator and the covariates on the right. 
#'   The time variables must be specified through the intCens() function. 
#'   See ?intCens and Details for further information. 
#'   
#' @param data A data.frame object. The data.frame in which to interpret the
#'   variable names in formula. Must contain the entry time, the left interval
#'   time, the right interval time, the vaccination time, and any covariates. 
#'   See Details.
#'   
#' @param constantVE A logical object. If FALSE (default), VE is assumed to be 
#'   potentially waning after the last change point; 
#'   otherwise it is assumed to be constant after the last change point.
#'   
#' @param plots A logical object. If TRUE (default), plots of the estimated 
#'   VE in reducing attack rate, the estimated VE in reducing the hazard rate, 
#'   and their 95\% confidence intervals will be automatically generated. 
#'   If FALSE, plots will not be generated.
#'   
#' @param changePts An integer vector object or NULL. The potential change
#'   points (in days) of the piece-wise log-linear hazard ratio. See Details 
#'   for further information. 
#'   If NULL, the Akaike information criterion (AIC) will be used to select 
#'   one change point from \{28, 35, 42, 49, 56\} (weeks 4, 5, 6, 7, 8).
#'   
#' @param timePts An integer vector object or NULL. The endpoints (in days) 
#'  of the time periods for which the VE in reducing the attack rate are to
#'  be estimated. The estimated VE in reducing the hazard rate at 
#'  these endpoints are also returned. If NULL, a default sequence 
#'  \eqn{t_1, 2t_1, 3t_1, \dots} will be used, where \eqn{t_1} is the first 
#'  change point. The sequence ends at the maximum of the left and
#'  and right ends of the time intervals from all participants. This input is 
#'  ignored when constantVE = TRUE.  
#'   
#' @param tol A numeric scalar object. The convergence threshold for the EM algorithm. 
#'   The default value is 0.0001. 
#'   
#' @param maxit A positive integer object. The maximum number of iterations for the 
#'   EM algorithm. The default value is 2000. 
#'   
#' @returns An S3 object of class iDOVE containing a list with elements
#'
#'   \item{call}{The unevaluated call.}
#'
#'   \item{changePts}{The changePts of the analysis.}
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
#'    named 'VE_a' contains the estimates of the VE in reducing the attack 
#'    rate at all time points given in timePts, together with the 95\%
#'    confidence intervals. The second matrix named 'VE_h' contains the 
#'    estimates of the VE in reducing the hazard rate at timePts. 
#'    The third matrix named 'VE_period' contains the estimates of VE in
#'    reducing the attack rate over successive time periods according to
#'    timePts, together with the 95\% confidence intervals.}
#'
#' @export
#' @import methods
#'
#' @useDynLib iDOVE
#' 
#' @include EMmeth.R verifyInputs.R
#' @import Rcpp
#' 
#' @examples
#' data(idoveData)
#'
#' set.seed(1234)
#' smp <- sample(1L:nrow(x = idoveData), size = 250L)
#' 
#' # NOTE: This sample size is chosen for example only -- larger data sets
#' # should be used.
#' # See the vignette for a full analysis of the idoveData dataset
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

# 8/11/21: remove CoxReg and additional attributes knots, gamma, covgamma, and tau.
idove <- function(formula,  
                  data, 
                  constantVE = FALSE,  
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

  inputs <- .verifyInputs(formula = formula,  
                          data = data, 
                          constantVE = constantVE,  
                          plots = plots, 
                          changePts = changePts,  
                          timePts = timePts, 
                          tol = tol,  
                          maxit = maxit)

  # changePts (in days) of linear splines 

  tau <- max(inputs$dt[,"left_time"], inputs$dt[,"right_time"])


  if (!is.null(x = inputs$changePts)) {
    def <- inputs$changePts
    use <- def < tau
    changePts <- list(def[use])
    if (length(def[use]) > 0L) {
      message("changePts: {", paste(def[use],collapse=","), "}")
    } else {
      message("all changePts > tau; using default setting")
      changePts <- NULL
    }
  }

  if (is.null(x = inputs$changePts)) {
    def <- (4L:8L)*7L
    use <- def < tau
    changePts <- as.list(x = def[use])
    message("changePts not given; using AIC to select from {",
            paste(def[use],collapse=","), "}")
  }

  ### main analysis
  
  res <- EMmeth(dt = inputs$dt,  
                X = inputs$X,  
                changePts = changePts,  
                constantVE = inputs$constantVE,  
                timePts = inputs$timePts,  
                plots = inputs$plots,  
                threshold = inputs$threshold,  
                maxit = inputs$maxit)
  
  res[["call"]] <- cl

  class(x = res) <- "iDOVE"
  
  return( res )
}
