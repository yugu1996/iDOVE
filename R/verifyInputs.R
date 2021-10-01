.verifyInputs <- function(formula,  
                          data, 
                          constantVE,  
                          plots, 
                          changePts,  
                          timePts, 
                          tol,  
                          maxit) {

  if (!is.logical(x = plots)) stop("plots must be logical", call. = FALSE)
  if (!is.numeric(x = tol)) stop("tol must be numeric", call. = FALSE)
  if (!is.numeric(x = maxit)) stop("maxit must be integer", call. = FALSE)
  maxit <- as.integer(x = round(x = maxit, digits = 0L))

  data <- .verifyData(formula = formula, data = data)

  changePts <- .verifyChangePts(changePts = changePts)

  timePts <- .verifyTimePts(timePts = timePts, constantVE = constantVE)

  return( list("dt" = data$data,  
               "X" = data$X,  
               "changePts" = changePts,  
               "constantVE" = constantVE,  
               "timePts" = timePts,  
               "plots" = plots,  
               "threshold" = tol,  
               "maxit" = maxit) )
  
}


#' @importFrom stats terms update.formula 
#' @importFrom stats model.frame model.matrix model.response
.verifyData <- function(formula, data) {

  # reset options to allow for keeping NA values
  opt <- options()
  options(na.action = 'na.pass')
  on.exit(options(opt))
  
  # add intercept to model if not provided
  # this was added 6/17/21 to ensure that factors are handled properly
  if (attr(x = stats::terms(x = formula), which = "intercept") == 0L) {
    formula = stats::update.formula(old = formula, new = .~. +1)
  }

  # try to obtain the model.frame
  mf <- tryCatch(expr = stats::model.frame(formula = formula, data = data),
                 error = function(e) {
                   message("unable to obtain model.frame")
                   stop(e$message, call. = FALSE)
                 })

  # extract covariates (X is a matrix)
  X <- suppressMessages(expr = stats::model.matrix(object = mf, data = data))
  # 6/17/21 keep only non-intercept
  int <- attr(x = X, which = "assign") != 0L
  X <- X[,int, drop = FALSE]

  if (ncol(x = X) == 0L) X <- NULL
  
  # extract time variables
  dt <- suppressMessages(expr = stats::model.response(data = mf))

  if (!is(object = dt, class2 = "intCens")) {
    stop("the LHS of formula did not contain an appropriate intCens() object",
         call. = FALSE)
  }

  # remove any cases that have NA in the entry times or X
  if (is.null(x = X)) {
    use <- !is.na(x = dt[,"entry_time"])
  } else {
    use <- stats::complete.cases(X) & !is.na(x = dt[,"entry_time"])
  }
  
  if (sum(use) == 0L) {
    stop("input checks result in all NA -- verify inputs", call. = FALSE)
  }

  dt <- dt[use,,drop=FALSE]
  if (!is.null(x = X)) X <- X[use,,drop=FALSE]
  
  if (sum(use) != length(x = use)) {
    message(sum(!use), 
            " cases removed from the analysis due to NA values")
  }

  rownames(x = X) <- NULL
  rownames(x = dt) <- NULL

  # apply the propose method for interval censored data
  # check if the data need to be modified due to no right-censored subjects
  if (max(dt[,"left_time"]) < max(dt[,"right_time"])) {
    tst <- dt[,"right_time"] > max(dt[,"left_time"])
    dt[tst, "right_time"] <- -1L
    message("all subjects whose right-interval time is greater than ",
            "the biggest left-interval time made to be right-censored.")
  }
    
  return( list("X" = X, "data" = dt) )
}

.verifyChangePts <- function(changePts) {

  if (is.null(x = changePts)) {
    return( NULL )
  }

  if (!is.numeric(x = changePts)) {
    stop("changePts must be a numeric vector or NULL", call. = FALSE)
  }
    
  if (length(x = changePts) == 0L) {
    message("changePts has zero length; using AIC to select from {28, 35, 42, 49, 56} days")
    return( NULL )
  } else if (any(changePts <= 0.0)) {
    stop("all changePts must be positive", call. = FALSE)
  }

  if (!is.integer(x = changePts)) {
    tmp <- as.integer(x = round(x = changePts, digits = 0L))
    if (!isTRUE(x = all.equal(target = changePts, current = tmp))) {
      stop("all changePts must be integer (days)", call. = FALSE)
    }
    changePts <- tmp
  }

  return( changePts )
}

.verifyTimePts <- function(timePts, constantVE) {

  if (!is.logical(x = constantVE) || is.na(x = constantVE)) {
    stop("constantVE must be logical", call. = FALSE)
  }

  if (constantVE) {
    message("constantVE selected")
    return( NULL )
  }

  if (is.null(x = timePts)) {
    message("timePts not given; default values will be used")
    return( NULL )
  }

  if (!is.numeric(x = timePts)) {
    stop("timePts must be a numeric vector or NULL", call. = FALSE)
  }
    
  if (length(x = timePts) == 0L) {
    message("timePts has zero length; default values will be used")
    return( NULL )
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

  return( timePts )
}

