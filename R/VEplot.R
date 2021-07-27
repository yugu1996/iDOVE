# Internal function: plot the curves of VE_a and VE_h 
# knots and tau are provided in days
.VEplot <- function(knots, tau, gamma, covgamma) {

  knots <- knots*0.0329
  lower <- c(0, knots)
  upper <- c(knots, Inf)

  timept = seq(0, tau*0.0329, 0.01)
  npt = length(timept)

  VE_h <- matrix(NA, nrow = npt, ncol = 4L)
  VE_a <- matrix(NA, nrow = npt, ncol = 4L)
  for (i in 1L:npt) {
    rt <- r(time = timept[i], 
            knots = knots, 
            gamma = gamma, 
            covgamma = covgamma)
    VE_h[i,1L:2L] <- c(1.0 - exp(x = rt[1L]), exp(x = rt[1L])*rt[2L])
    VE_h[i,3L] <- 1.0 - exp(x = rt[1L]+1.96*rt[2L])
    VE_h[i,4L] <- 1.0 - exp(x = rt[1L]-1.96*rt[2L])

    Vt <- V(time = timept[i], 
            lower = lower,  
            upper = upper,  
            gamma = gamma,  
            covgamma = covgamma)
    VE_a[i,1L:2L] <- c(1.0 - Vt[1L]/timept[i], Vt[2L]/timept[i])
    VC <- 1.0 - VE_a[i,1L]
    VE_a[i,3L] <- 1.0 - VC*exp(x = 1.96*VE_a[i,2L]/VC)
    VE_a[i,4L] <- 1.0 - VC*exp(x = -1.96*VE_a[i,2L]/VC)
  }
  VE_h <- cbind(timept/0.0329, VE_h)
  VE_a <- cbind(timept/0.0329, VE_a)
  VE_a[1L,] <- VE_h[1L,]
  
  ymin.a <- min(VE_a[,c(2L,4L,5L)]*100.0)
  ymin.h <- min(VE_h[,c(2L,4L,5L)]*100.0)
  ymax.a <- max(VE_a[,c(2L,4L,5L)]*100.0)
  ymax.h <- max(VE_h[,c(2L,4L,5L)]*100.0)
  ymin.a <- ymin.a - ymin.a %% 10L - 10.0
  ymin.h <- ymin.h - ymin.h %% 10L - 10.0
  ymax.a <- min(ymax.a - ymax.a %% 10L + 20.0, 100.0)
  ymax.h <- min(ymax.h - ymax.h %% 10L + 20.0, 100.0)
  
  graphics::par(mfrow=c(1L,2L))
  
  graphics::plot(VE_a[,1L], VE_a[,2L]*100,
                 xlab = "Time since vaccination (in days)", 
                 ylab = "Vaccine efficiency in reducing attack rate (%)", 
                 type = 'l',
                 yaxt="n", main="", ylim = c(ymin.a, ymax.a))
  graphics::lines(VE_a[,1L], VE_a[,4L]*100, col = 3L)
  graphics::lines(VE_a[,1L], VE_a[,5L]*100, col = 3L)
  graphics::axis(side = 2L,
                 at = seq(ymin.a, ymax.a, 20),
                 labels = paste0(seq(ymin.a, ymax.a, 20), "%"), 
                 cex = 0.8)
  graphics::legend(x = "bottom",
                   legend = c("VE_a", "95% CI"), 
                   lty = c(1L,1L),
                   col = c(1L,3L), 
                   bg = "gray95",
                   cex = 0.75)
  
  graphics::plot(VE_h[,1L], VE_h[,2L]*100,
                 xlab = "Time since vaccination (in days)", 
                 ylab = "Vaccine efficiency in reducing hazard rate (%)", 
                 type = 'l',
                 yaxt="n", main="", ylim = c(ymin.h, ymax.h))
  graphics::lines(VE_h[,1L], VE_h[,4L]*100, col = 3L)
  graphics::lines(VE_h[,1L], VE_h[,5L]*100, col = 3L)
  graphics::axis(side = 2L,
                 at = seq(ymin.h, ymax.h, 20),
                 labels = paste0(seq(ymin.h, ymax.h, 20), "%"), 
                 cex = 0.8)
  graphics::legend(x = "bottom",
                   legend = c("VE_h", "95% CI"), 
                   lty = c(1L,1L),
                   col = c(1L,3L), 
                   bg = "gray95",
                   cex = 0.75)
  
  return()
}
