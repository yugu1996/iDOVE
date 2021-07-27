#' Specify Time Variables
#'
#' This function is used in the model statement of idove() to specify
#'  the entry time, left interval time, right interval time, and vaccination time.
#'
#' @param entry_time The variable for the time when 
#'   the participant enters the trial. 
#'   
#' @param left_time The variable for the last examination time when
#'   the test is negative
#'   
#' @param right_time The variable for the first examination time when
#'   the test is positive; Inf or NA if the participant is never tested
#'   positive during the follow-up
#'
#' @param vaccination_time The variable for the time when 
#'   vaccination takes place. 
#'
#' @returns This function is intended to be used only in the model statement 
#'  of idove(). The result, a matrix, is used internally.
#'    
#' @name intCens
#' @rdname intCens
#' @export

intCens = function(entry_time, left_time, right_time, vaccination_time) {
  
  ### entry time
  
  # must be provided as a numeric vector. 
  
  if (missing(x = entry_time)) {
    stop("must provide a entry_time argument", 
         call. = FALSE)
  }
  
  if (!is.numeric(x = entry_time)) {
    stop ("entry_time is not numeric", call. = FALSE)
  }

  if (!is.integer(x = entry_time)) {
    tmp <- as.integer(x = round(x = entry_time, digits = 0L))
    if (!isTRUE(x = all.equal(target = entry_time, current = tmp))) {
      stop("all times must be integer (days)", call. = FALSE)
    }
    entry_time <- tmp
  }
  
  ### left interval time
  
  # must be provided as a numeric vector. 
  
  if (missing(x = left_time)) {
    stop("must provide a left_time argument", 
         call. = FALSE)
  }
  
  if (!is.numeric(x = left_time)) {
    stop ("left_time is not numeric", call. = FALSE)
  }
  
  if (!is.integer(x = left_time)) {
    tmp <- as.integer(x = round(x = left_time, digits = 0L))
    if (!isTRUE(x = all.equal(target = left_time, current = tmp))) {
      stop("all times must be integer (days)", call. = FALSE)
    }
    left_time <- tmp
  }

  ### right interval time
  
  # must be provided as a numeric vector. 
  
  if (missing(x = right_time)) {
    stop("must provide a right_time argument", 
         call. = FALSE)
  }
  
  if (!is.numeric(x = right_time)) {
    stop ("right_time is not numeric", call. = FALSE)
  }
  
  # right times that are provided as infinite are set to -1 to allow
  # for integer conversion
  id <- is.infinite(x = right_time) | is.na(x = right_time)
  right_time[id] <- -1

  if (!is.integer(x = right_time)) {
    tmp <- as.integer(x = round(x = right_time, digits = 0L))
    if (!isTRUE(x = all.equal(target = right_time, current = tmp))) {
      stop("all times must be integer (days)", call. = FALSE)
    }
    right_time <- tmp
  }
  ### time of vaccination
  
  # must be provided as a numeric vector. 
  
  if (missing(x = vaccination_time)) {
    stop("must provide a vaccination_time argument", 
         call. = FALSE)
  }
  
  if (!is.numeric(x = vaccination_time)) {
    stop ("vaccination_time is not numeric", call. = FALSE)
  }
  
  if (!is.integer(x = vaccination_time)) {
    tmp <- as.integer(x = round(x = vaccination_time, digits = 0L))
    if (!isTRUE(x = all.equal(target = vaccination_time, current = tmp))) {
      stop("all times must be integer (days)", call. = FALSE)
    }
    vaccination_time <- tmp
  }

  ### check if all inputs are of the same length
  
  if (length(x = entry_time) != length(x = left_time)) {
    stop("all inputs must be of same length", call. = FALSE)
  }
  
  if (length(x = entry_time) != length(x = right_time)) {
    stop("all inputs must be of same length", call. = FALSE)
  }
  
  if (length(x = entry_time) != length(x = vaccination_time)) {
    stop("all inputs must be of same length", call. = FALSE)
  }
  
  ### check if entry_time <= left_time <= right_time
  
  tst <- {{entry_time > left_time} | {left_time > right_time}} & right_time >= 0
  
  if (any(tst, na.rm = TRUE)) {
    message(sum(tst, na.rm = TRUE), 
            " don't satisfy entry_time <= left_time <= right_time; cases set to NA")
    entry_time[tst] <- NA
  }
  
  dm <- cbind(entry_time, left_time, right_time, vaccination_time)
  
  cname <- c("entry_time", "left_time", "right_time", "vaccination_time")
  dimnames(x = dm) <- list(NULL, cname)
  
  return( dm )
}
