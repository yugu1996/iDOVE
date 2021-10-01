#' @importFrom stats pnorm
#' @include VEcal.R VEplot.R
# 8/11/21: set timept to 0,1,...tau, compute daily VE estimates
# and generate the plots based on VE_a and VE_h (lines 85-118, 161-166)

.postProcess <- function(theta, 
                         covMat, 
                         plots, 
                         SD, 
                         varname, 
                         constantVE, 
                         tau,
                         changePts,
                         timePts) {

  nTheta <- length(x = theta)

  nBeta <- length(x = SD)

  nGam <- nTheta - nBeta

  beta.output <- .postBeta(theta = theta, 
                           SD = SD,  
                           covMat = covMat,  
                           varname = varname)

  gamma <- theta[{nBeta+1L}:nTheta]
  covgamma <- as.matrix(x = covMat[{nBeta+1L}:nTheta, {nBeta+1L}:nTheta])

  if (constantVE) {

    VE.output <- .postGammaConstant(knots = changePts, 
                                    gamma = gamma, 
                                    covgamma = covgamma)
   
  } else {
    
    VE.output <- .postGammaNotConstant(timePts = timePts, 
                                       knots = changePts,  
                                       tau = tau,  
                                       gamma = gamma,  
                                       covgamma = covgamma,
                                       plots = plots)
    
  }

  return( list("covariates" = beta.output,
               "vaccine" = VE.output) )

}

.postBeta <- function(theta, SD, covMat, varname) {

  nBeta <- length(x = SD)

  if (nBeta == 0L) return( NA )

  beta <- theta[1L:nBeta]/SD
  covbeta <- covMat[1L:nBeta, 1L:nBeta, drop = FALSE]

  # scale back beta and its SE
  sebeta <- sqrt(x = diag(x = covbeta))/SD
  zbeta <- beta/sebeta

  beta.output <- cbind("coef" = beta, 
                       "se(coef)" = sebeta,
                       "z" = zbeta,
                       "Pr(>|z|)" = 2.0*pnorm(q = abs(x = zbeta), lower.tail = FALSE),
                       "exp(coef)" = exp(x = beta),
                       "lower .95" = exp(x = beta-1.96*sebeta),
                       "upper .95" = exp(x = beta+1.96*sebeta))
 
  rownames(x = beta.output) <- varname

  return( beta.output )  

}

.postGammaConstant <- function(knots, gamma, covgamma) {

  npc <- length(x = knots)

  B <- knots[npc] - c(0, knots[-npc])
  B <- B*0.0329

  HRest <- exp(x = sum(B * gamma))
  HRsdest <- sqrt(x = HRest^2 * {B %*% covgamma %*% B})

  VE.out <- c("VE" = 1.0 - HRest,
              "se" = HRsdest,
              "lower .95" = 1.0 - HRest*exp(x = 1.96*HRsdest/HRest),
              "upper .95" = 1.0 - HRest*exp(x = -1.96*HRsdest/HRest))
    
  VE.output = list("VE" = VE.out)
 
  return( VE.output )   
}

.postGammaNotConstant <- function(timePts, knots, tau, gamma, covgamma, plots) {

  if (is.null(x = timePts)) {
    timePts <- seq(from = knots[1L], to = tau, by = knots[1L])
  } else {
    if (any(timePts > tau)) {
      message("timePts > tau have been removed")
      timePts <- timePts[timePts < tau]
    }
  }
    
  # express knots in months
  knots <- knots*0.0329

  timept <- {0L:tau}*0.0329

  VE_h <- .VEh(timept = timept, 
               knots = knots, 
               gamma = gamma, 
               covgamma = covgamma)

  VE_a <- .VEa(timept = timept, 
               knots = knots, 
               gamma = gamma, 
               covgamma = covgamma)

  VE_a[1L,] <- VE_h[1L,]

  VE_ave <- .VEave(timePts = timePts, 
                   knots = knots, 
                   gamma = gamma, 
                   covgamma = covgamma)

  if (plots) .VEplot(VE_a = VE_a, VE_h = VE_h)
    
  return( list("VE_a" = VE_a,
               "VE_h" = VE_h,
               "VE_period" = VE_ave) )
}

.VEa <- function(timept, knots, gamma, covgamma) {

  # express knots in months
  lower <- c(0.0, knots)
  upper <- c(knots, Inf)
    
  npt <- length(x = timept)

  VE_a <- matrix(data = NA, nrow = npt, ncol = 5L)
  colnames(x = VE_a) <- c("time", "VE_a", "se", "lower .95", "upper .95")

  Vt <- sapply(X = timept, 
               FUN = .V,
               lower = lower,
               upper = upper,
               gamma = gamma,
               covgamma = covgamma)

  VE_a[,1L] <- timept/0.0329
  VE_a[,2L] <- 1.0 - Vt[1L,]/timept
  VE_a[,3L] <- Vt[2L,]/timept
  VC <- Vt[1L,]/timept
  VE_a[,4L] <- 1.0 - VC*exp(x =  1.96*VE_a[,3L]/VC)
  VE_a[,5L] <- 1.0 - VC*exp(x = -1.96*VE_a[,3L]/VC)

  return( VE_a )

}

.VEh <- function(timept, knots, gamma, covgamma) {

  npt <- length(x = timept)

  # VE_h and VE_a: Est | SE | Lower | Upper
  VE_h <- matrix(data = NA, nrow = npt, ncol = 5L)
  colnames(x = VE_h) <- c("time", "VE_h", "se", "lower .95", "upper .95")

  rt <- sapply(X = timept, 
               FUN = .r,
               knots = knots,
               gamma = gamma,
               covgamma = covgamma)

  VE_h[,1L] <- timept/0.0329
  VE_h[,2L] <- 1.0 - exp(x = rt[1L,])
  VE_h[,3L] <- exp(x = rt[1L,])*rt[2L,]
  VE_h[,4L] <- 1.0 - exp(x = rt[1L,]+1.96*rt[2L,])
  VE_h[,5L] <- 1.0 - exp(x = rt[1L,]-1.96*rt[2L,])

  return( VE_h )

}

.VEave <- function(timePts, knots, gamma, covgamma) {

  lower <- c(0.0, knots)
  upper <- c(knots, Inf)

  nGam <- length(x = gamma)

  # VE_a over time periods
  # 8/12/21: reset intt
  # intt <- timept
  intt <- sort(x = unique(x = c(0.0,timePts)))*0.0329
  nint <- length(x = intt)
  dint <- diff(x = intt)

  VE_ave <- matrix(data = NA, nrow = nint-1L, ncol = 6L)
  VE_ave[,1L] <- intt[-nint]/0.0329
  VE_ave[,2L] <- intt[-1L]/0.0329
  colnames(x = VE_ave) <- c("left", "right", "VE_a", "se", 
                            "lower .95", "upper .95")

  Vt <- sapply(X = intt, 
               FUN = .V,
               lower = lower,
               upper = upper,
               gamma = gamma,
               covgamma = covgamma)

  der <- sapply(X = intt, 
                FUN = function(x, lower, upper, gamma, covgamma) {
                        .partial(time = x, 
                                 lower = lower,  
                                 upper = upper,
                                 gamma = gamma)[[ 2L ]]
                     },
                lower = lower,
                upper = upper,
                gamma = gamma)

  VE_ave[,3L] <- 1.0 - diff(x = Vt[1L,])/dint
  tempvar <- Vt[2L,-nint]^2 + Vt[2L,-1L]^2 - 
             2.0 * diag(x = t(x = der[,-nint,drop=FALSE]) %*% 
                            covgamma %*% der[,-1L,drop=FALSE])
  VE_ave[,4L] <- sqrt(x = tempvar)/dint
  VC_ave <- 1.0 - VE_ave[,3L]
  VE_ave[,5L] <- 1.0 - VC_ave*exp(x = 1.96*VE_ave[,4L]/VC_ave)
  VE_ave[,6L] <- 1.0 - VC_ave*exp(x = -1.96*VE_ave[,4L]/VC_ave)
  
  return( VE_ave )
}
