# Internal functions for calculating VE_a, VE_h and VE_a over a time period
# as well as their SE and 95% CI
# t and knots are in months

r <- function(time, knots, gamma, covgamma) {
  # calculate r(t) and its SE
  B <- pmax(time-c(0.0,knots),0.0)
  estr <- B %*% gamma
  sdr <- sqrt(B %*% covgamma %*% B)
  return( c(estr, sdr) )
}

partial <- function(time, lower, upper, gamma) {

  # calculate V(t) and V'(t)
  cumgamma <- cumsum(x = gamma)
  cumgammay <- cumsum(x = gamma*lower)

  ebig <- exp(x = cumgamma*pmin(time,upper))
  esmall <- exp(x = cumgamma*lower)
  ediff <- ebig - esmall

  coef <- {1.0/cumgamma}*exp(x = -cumgammay)
  coeft <- {1.0/cumgamma^2}*ediff

  Vt <- coef[time > lower] %*% ediff[time > lower]

  npc <- length(x = gamma)
  der <- rep(x = 0.0, times = npc)
  te <- pmin(time, upper)*ebig - lower*esmall
  for (k in 1L:npc) {
    if (time > lower[k]) {
      tediff <- te - lower[k]*{ebig-esmall}
      coef2 <- {1.0/cumgamma}*tediff - coeft
      tst <- time>lower & lower>=lower[k]
      der[k] <- coef2[tst] %*% exp(x = -cumgammay[tst])
    }
  }
  return( list("Vt" = Vt, "pVt" = der) )
}

V <- function(time, lower, upper, gamma, covgamma) {
  # calculate V(t) and its SE
  temp <- partial(time = time, lower = lower, upper = upper, gamma = gamma)
  estV <- temp[[ 1L ]]
  der <- temp[[ 2L ]]
  sdV <- sqrt(x = der %*% covgamma %*% der)
  return( c(estV, sdV) )
}
