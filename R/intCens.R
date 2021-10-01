#' Specify Time Variables
#'
#' This function is used in the model statement of idove() to specify
#'  the entry time, left interval time, right interval time, and vaccination 
#'  time.
#'
#' Times must obey the following relationships:
#'  (i) For all participants, entry_time <= left_time;
#'  (ii) For all participants that tested positive during the trial,
#'  entry_time <= left_time <= right_time; and 
#'  (iii) For all participants that received vaccination,
#'  entry_time <= vaccination_time. If a case is found to violate one or more
#'  of these relationships, its entry_time is set to NA.
#'
#' @param entry_time The variable for the time when 
#'   the participant enters the trial. Entry times must be integer
#'   (or be able to be cast as integer without loss of information),
#'   non-negative, and complete.
#'   
#' @param left_time The variable for the last examination time when
#'   the test is negative. Left interval times must be integer
#'   (or be able to be cast as integer without loss of information),
#'   non-negative, and complete.
#'   
#' @param right_time The variable for the first examination time when
#'   the test is positive. Right interval times must be integer
#'   (or be able to be cast as integer without loss of information),
#'   non-negative, and NA or Inf should be used if the participant never 
#'   tested positive during the trial.
#'
#' @param vaccination_time The variable for the time when 
#'   vaccination takes place. Vaccination times must be integer
#'   (or be able to be cast as integer without loss of information),
#'   non-negative, and NA or Inf should be used if the participant was not 
#'   vaccinated during the trial.
#'
#' @returns This function is intended to be used only in the model statement 
#'  of idove(). The result, a matrix, is used internally.
#'    
#' @name intCens
#' @rdname intCens
#' @export
#
# Note to future developer: Function sets Infs to NAs in right_time and
# vaccination_time. Sets entry_time to NA if E <= L <= R is not satisfied.
# 8/11/21: vaccination time can be NA or Inf if a participant is not vaccinated 
# during the clinical trial (lines 112-115). 
# 8/16/21: simplification of code; adjusted message and stop statements

intCens = function(entry_time, left_time, right_time, vaccination_time) {
  
  if (missing(x = entry_time) || missing(x = left_time) ||
      missing(x = right_time) || missing(x = vaccination_time)) {
    stop("intCens inputs must include entry_time, left_time, right_time, and, ", 
         "vaccination_time", call. = FALSE)
  }

  ### entry time
  
  # must be provided as a numeric vector.
  # must be integer-like
  # must be complete
  # must be non-negative
  
  entry_time <- .basicTests_noNA(x = entry_time, name = "entry_time")

  ### left interval time
  
  # must be provided as a numeric vector.
  # must be integer-like
  # must be complete
  # must be non-negative
  
  left_time <- .basicTests_noNA(x = left_time, name = "left_time")

  ### right interval time
  
  # must be provided as a numeric vector.
  # must be integer-like
  # can be incomplete (NA or Inf, Inf replaced with NA)
  # must be non-negative
  
  right_time <- .basicTests_NAorInf(x = right_time, name = "right_time")
  right_time[is.na(x = right_time)] <- -1L

  ### time of vaccination
  
  # must be provided as a numeric vector.
  # must be integer-like
  # can be incomplete (NA or Inf)
  # must be non-negative
  
  vaccination_time <- .basicTests_NAorInf(x = vaccination_time, 
                                          name = "vaccination_time")
  vaccination_time[is.na(x = vaccination_time)] <- 999999L
  
  ### check if all inputs are of the same length
  if (length(x = entry_time) != length(x = left_time) || 
      length(x = entry_time) != length(x = right_time) ||
      length(x = entry_time) != length(x = vaccination_time)) {
    stop("intCens inputs must be of same length", call. = FALSE)
  }
  
  ### check if entry_time <= left_time <= right_time and
  ### entry_time <= vaccination_time
  
  tst <- entry_time <= left_time

  eventStatus <- right_time >= 0L
  if (sum(eventStatus) == 0L) {
    stop("no events in the data", call. = FALSE)
  }

  tst[eventStatus] <- tst[eventStatus] & 
                      {left_time[eventStatus] <= right_time[eventStatus]}

  vacStatus <- vaccination_time < 999999L

  if (sum(vacStatus) == 0L) {
    stop("no vaccinated participants", call. = FALSE)
  }

  tst[vacStatus] <- tst[vacStatus] & 
                    {entry_time[vacStatus] <= vaccination_time[vacStatus]}

  violate <- length(x = entry_time) - sum(tst)
 
  if (violate > 0L) {
    message(violate, 
            ifelse(test = violate > 1L, 
                   yes = " cases do not ", 
                   no = " case does not "),
            "satisfy required time relationships; ",
            ifelse(test = violate > 1L, 
                   yes = "cases removed", 
                   no = "case removed"))
    entry_time[!tst] <- NA
  }

  if (sum(eventStatus & vacStatus) == 0L) {
    stop("no vaccinated participants experienced an event; ", 
         "vaccine efficacy cannot be estimated", call. = FALSE)
  }

  if (all(left_time[eventStatus] == right_time[eventStatus])) {
    warning("all event times are either exactly observed or right-censored; ",
            "the user is recommended to use dove2() instead for more stable ",
            "variance estimation")
  }

  dm <- cbind(entry_time, left_time, right_time, vaccination_time)
  
  cname <- c("entry_time", "left_time", "right_time", "vaccination_time")
  dimnames(x = dm) <- list(NULL, cname)

  class(x = dm) <- c("intCens", class(x = dm))

  return( dm )
}

# test of input vector that must satisfy:
#   must be provided as a numeric vector.
#   must be integer-like
#   must be complete -- no NA or Inf
#   must be non-negative
.basicTests_noNA <- function(x, name) {

  if (anyNA(x = x)) {
    stop(name, " must be complete", call. = FALSE)
  }

  if (any(is.infinite(x = x))) {
    stop(name, " contains infinities", call. = FALSE)
  }
  
  if (is.numeric(x = x)) {
    if (!is.integer(x = x)) {
      tmp <- as.integer(x = round(x = x, digits = 0L))
      if (!isTRUE(x = all.equal(target = x, current = tmp))) {
        stop(name, " is not integer (days)", call. = FALSE)
      }
      x <- tmp
    }
  } else {
    stop(name, " is not integer (days)", call. = FALSE)
  }

  if (any(x < 0L)) {
    stop(name, " must be non-negative", call. = FALSE)
  }

  return( x )
}


# test of input vector that must satisfy:
#   must be provided as a numeric vector.
#   must be integer-like
#   can be incomplete (NA or Inf w/ Inf replaced by NA)
#   must be non-negative
.basicTests_NAorInf <- function(x, name) {

  x[is.infinite(x = x)] <- NA

  if (all(is.na(x = x))) return( x )

  if (is.numeric(x = x)) {
    if (!is.integer(x = x)) {
      tmp <- as.integer(x = round(x = x, digits = 0L))
      if (!isTRUE(x = all.equal(target = x, current = tmp))) {
        stop(name, " is not integer (days)", call. = FALSE)
      }
      x <- tmp
    }
  } else {
    stop(name, " is not integer (days)", call. = FALSE)
  }

  if (any(x < 0L, na.rm = TRUE)) {
    stop(name, " must be non-negative", call. = FALSE)
  }

  return( x )
}
