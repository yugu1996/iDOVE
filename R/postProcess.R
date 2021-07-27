#' @importFrom stats pnorm
#' @include VEcal.R VEplot.R
.postProcess <- function(theta, 
                         covMat, 
                         plots, 
                         SD, 
                         varname, 
                         constantVE, 
                         tau,
                         knots,
                         timePts) {

  nTheta <- length(x = theta)

  if (is.null(x = SD)) {
    nBeta <- 0L
  } else {
    nBeta <- length(x = SD)
  }

  nGam <- nTheta - nBeta

  if (nBeta > 0L) {
  
    beta <- theta[1L:nBeta]/SD
    covbeta <- as.matrix(covMat[1L:nBeta, 1L:nBeta])

    # scale back beta and its SE
    sebeta <- sqrt(x = diag(x = covbeta))/SD
    zbeta <- beta/sebeta
  
    beta.output <- cbind(beta, 
                         sebeta,
                         zbeta,
                         2*pnorm(q = abs(x = zbeta), lower.tail = FALSE),
                         exp(x = beta),
                         exp(x = beta-1.96*sebeta),
                         exp(x = beta+1.96*sebeta))
 
    colnames(x = beta.output) <- c("coef",
                                   "se(coef)", "z", "Pr(>|z|)",
                                   "exp(coef)",
                                   "lower .95", "upper .95")
   
    rownames(x = beta.output) <- varname
  } else {
    beta.output <- NA
  }
  
  gamma <- theta[{nBeta+1L}:nTheta]
  covgamma <- as.matrix(covMat[{nBeta+1L}:nTheta, {nBeta+1L}:nTheta])

  if (constantVE) {
    npc <- length(x = knots)

    B <- knots[npc] - c(0, knots[-npc])
    B <- B*0.0329

    HRest <- exp(x = sum(B * gamma))
    HRcov <- HRest^2 * {B %*% covgamma %*% B}
    HRsdest <- sqrt(x = HRcov)
    VE.lower <- 1.0 - HRest*exp(x = 1.96*HRsdest/HRest)
    VE.upper <- 1.0 - HRest*exp(x = -1.96*HRsdest/HRest)
    VE.out <- c("VE" = 1.0 - HRest,
                "se" = HRsdest,
                "lower .95" = VE.lower,
                "upper .95" = VE.upper)
    
    VE.output = list("VE" = VE.out)
    
  } else {
    
    if (is.null(x = timePts)) {
      timePts = seq(from = knots[1L], to = tau, by = knots[1L])
    }
    
    # express knots in months
    knots <- knots*0.0329
    lower <- c(0, knots)
    upper <- c(knots, Inf)
    
    timept <- sort(x = unique(x = c(0.0,timePts)))*0.0329
    npt <- length(x = timept)

    # VE_h and VE_a: Est | SE | Lower | Upper
    VE_h <- matrix(NA, nrow = npt-1L, ncol = 4L)
    VE_a <- matrix(NA, nrow = npt-1L, ncol = 4L)
    for (i in 2L:npt) {
      rt <- r(time = timept[i], 
              knots = knots, 
              gamma = gamma, 
              covgamma = covgamma)
      VE_h[i-1L,1L:2L] <- c(1.0 - exp(x = rt[1L]), exp(x = rt[1L])*rt[2L])
      VE_h[i-1L,3L] <- 1.0 - exp(x = rt[1L]+1.96*rt[2L])
      VE_h[i-1L,4L] <- 1.0 - exp(x = rt[1L]-1.96*rt[2L])

      Vt <- V(time = timept[i], 
              lower = lower,  
              upper = upper,  
              gamma = gamma,  
              covgamma = covgamma)
      VE_a[i-1L,1L:2L] <- c(1.0 - Vt[1L]/timept[i], Vt[2L]/timept[i])
      VC <- 1.0 - VE_a[i-1L,1L]
      VE_a[i-1L,3L] <- 1.0 - VC*exp(x = 1.96*VE_a[i-1L,2L]/VC)
      VE_a[i-1L,4L] <- 1.0 - VC*exp(x = -1.96*VE_a[i-1L,2L]/VC)
    }
    VE_h <- cbind(timept[-1L]/0.0329, VE_h)
    VE_a <- cbind(timept[-1L]/0.0329, VE_a)
    
    colnames(x = VE_a) <- c("time", "VE_a", "se", "lower .95", "upper .95")
    colnames(x = VE_h) <- c("time", "VE_h", "se", "lower .95", "upper .95")

    # VE_a over time periods
    intt <- timept
    nint <- length(x = intt) - 1L
    intl <- intt[1L:nint]
    intu <- intt[2L:{nint+1L}]
    dint <- intu - intl
    VE_ave <- matrix(data = NA, nrow = nint, ncol = 4L)

    if (!anyNA(gamma) & !anyNA(covgamma)) {

      VV <- matrix(data = 0.0, nrow = nint + 1L, ncol = 2L)
      der <- NULL
      for (i in 1L:{nint+1L}) {
        VV[i,] <- V(time = intt[i], 
                    lower = lower, 
                    upper = upper, 
                    gamma = gamma, 
                    covgamma = covgamma)

        der <- cbind(der, 
                     partial(time = intt[i], 
                             lower = lower, 
                             upper = upper, 
                             gamma = gamma)[[ 2L ]])
      }

      for (i in 1L:nint) {
        VE_ave[i,1L] <- 1.0 - {VV[i+1L,1L]-VV[i,1L]}/dint[i]
        tempvar <- VV[i,2L]^2 + VV[i+1L,2L]^2 - 
                   2.0*{der[,i] %*% covgamma %*% der[,i+1L]}
        VE_ave[i,2L] <- sqrt(x = tempvar)/dint[i]
        VC_ave <- 1.0 - VE_ave[i,1L]
        VE_ave[i,3L] <- 1.0 - VC_ave*exp(x = 1.96*VE_ave[i,2L]/VC_ave)
        VE_ave[i,4L] <- 1.0 - VC_ave*exp(x = -1.96*VE_ave[i,2L]/VC_ave)
      }
    }

    VE_ave <- cbind(intl/0.0329, intu/0.0329, VE_ave)
    colnames(x = VE_ave) <- c("left", "right", "VE_a", "se", 
                              "lower .95", "upper .95")
    
    if (plots) {
      .VEplot(knots = knots/0.0329, 
              tau = tau, 
              gamma = gamma, 
              covgamma = covgamma)
    }
    
    VE.output <- list("VE_a" = VE_a,
                      "VE_h" = VE_h,
                      "VE_period" = VE_ave)
    
  }

  return( list("covariates" = beta.output,
               "vaccine" = VE.output) )

}
