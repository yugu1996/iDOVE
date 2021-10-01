# Internal functions for calculating VE_a, VE_h and VE_a over a time period
# as well as their SE and 95% CI
# t and knots are in months
# 6/30/2021 extracted from iDOVE v1.2 without modification
# @param a numeric scalar
# @param a numeric vector of length n_gamma - 1
# @param a numeric vector
# @param a numeric matrix of dim {n_gamma x n_gamma}
.r <- function(time, knots, gamma, covgamma) {

  if (!is.numeric(x = time) || !is.numeric(x = knots) ||
      !is.numeric(x = gamma) || !is.numeric(x = covgamma)) {
    stop("inappropriate input type to r(); contact maintainer", call. = FALSE)
  }

  if (!is.vector(x = time) ||
      !is.vector(x = knots) ||
      !is.vector(x = gamma) || 
      !is.matrix(x = covgamma)) {
    stop("inappropriate input to r(); contact maintainer", call. = FALSE)
  }

  if (is.na(x = time) || any(is.na(x = knots)) ||
      any(is.na(x = gamma)) || any(is.na(x = covgamma))) {
    stop("NAs found in r() inputs; contact maintainer", call. = FALSE)
  }

  if ({length(x = time) != 1L} ||
      {{length(x = knots) + 1L} != length(x = gamma)} ||
      {length(x = gamma) != ncol(x = covgamma)} ||
      {length(x = gamma) != nrow(x = covgamma)}) {
    stop("inappropriate input dimensions for r(); contact maintainer", 
         call. = FALSE)
  }

  # calculate r(t) and its SE
  B <- pmax(time-c(0.0,knots),0.0)
  estr <- B %*% gamma
  sdr <- sqrt(B %*% covgamma %*% B)

  return( c(estr, sdr) )
}

# @param time a numeric scalar
# @param lower a numeric vector of length n_gamma
# @param upper a numeric vector of length n_gamma
# @param gamma a numeric vector
.partial <- function(time, lower, upper, gamma) {

  if (!is.numeric(x = time) || !is.numeric(x = lower) ||
      !is.numeric(x = gamma) || !is.numeric(x = upper)) {
    stop("inappropriate input type to partial; contact maintainer", call. = FALSE)
  }

  if (is.na(x = time) || any(is.na(x = lower)) ||
      any(is.na(x = gamma)) || any(is.na(x = upper))) {
    stop("NAs found in partial inputs; contact maintainer", call. = FALSE)
  }

  if (!is.vector(x = time) ||
      !is.vector(x = lower) ||
      !is.vector(x = upper) || 
      !is.vector(x = gamma)) {
    stop("inappropriate input to partial; contact maintainer", call. = FALSE)
  }

  if (length(x = time) != 1L || 
      {length(x = upper) != length(x = lower)} ||
      {length(x = gamma) != length(x = lower)}) {
    stop("inappropriate input dimensions for partial; contact maintainer", 
         call. = FALSE)
  }

  if (!any(time > lower)) {
    return( list("Vt" = 0.0, "pVt" = rep(x = 0.0, times = length(x = gamma))) )
  }

  # calculate V(t) and V'(t)
  cumgamma <- cumsum(x = gamma)
  cumgammay <- cumsum(x = gamma*lower)

  ebig <- exp(x = cumgamma*pmin(time,upper))
  esmall <- exp(x = cumgamma*lower)
  ediff <- ebig - esmall

  coef <- {1.0/cumgamma}*exp(x = -cumgammay)
  coeft <- {1.0/cumgamma^2}*ediff

  Vt <- coef[time > {lower+1e-8}] %*% ediff[time > {lower+1e-8}]

  npc <- length(x = gamma)
  der <- rep(x = 0.0, times = npc)
  te <- pmin(time, upper)*ebig - lower*esmall
  for (k in 1L:npc) {
    if (time > {lower[k]+1e-8}) {
      tediff <- te - lower[k]*ediff
      coef2 <- {1.0/cumgamma}*tediff - coeft
      tst <- {time > lower} & {lower > {lower[k]-1e-8}}
      der[k] <- coef2[tst] %*% exp(x = -cumgammay[tst])
    }
  }
  return( list("Vt" = Vt, "pVt" = der) )
}

.V <- function(time, lower, upper, gamma, covgamma) {

  if (!is.matrix(x = covgamma)) {
    stop("inappropriate input to V(); contact maintainer", call. = FALSE)
  }

  # calculate V(t) and its SE
  temp <- .partial(time = time, lower = lower, upper = upper, gamma = gamma)
  estV <- temp[[ 1L ]]
  der <- temp[[ 2L ]]

  if ({length(x = der) != ncol(x = covgamma)} ||
      {length(x = der) != nrow(x = covgamma)}) {
    stop("inappropriate input dimensions for V(); contact maintainer", 
         call. = FALSE)
  }

  sdV <- sqrt(x = der %*% covgamma %*% der)

  return( c(estV, sdV) )
}
