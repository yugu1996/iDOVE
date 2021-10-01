# Internal function: plot the curves of VE_a and VE_h 
# knots and tau are provided in days
# @param VE_a is an {npt x 5} matrix with colums
#   "time", "VE_a", "se", "lower .95", "upper .95"
# @param VE_h is an {npt x 5} matrix with colums
#   "time", "VE_h", "se", "lower .95", "upper .95"
#
#' @importFrom grDevices dev.new
# 6/30/2021 introduced with new dove2() function in v1.7
# 8/12/21: change the input to be VE_a and VE_h
# 8/20/21: broke out .graph function
.VEplot <- function(VE_a, VE_h) {

  .graph(x = VE_a, 
         ylab = "Vaccine Efficiency in Reducing Attack Rate (%)",
         leg = "VE_a")
  
  grDevices::dev.new()

  .graph(x = VE_h, 
         ylab = "Vaccine Efficiency in Reducing Hazard Rate (%)",
         leg = "VE_h")
  
  return()
}

#' @importFrom graphics par plot lines axis legend 
.graph <- function(x, ylab, leg) {

  if (!is.matrix(x = x)) {
    stop("inappropriate input type to VEplot_type2", call. = FALSE)
  }

  if (ncol(x = x) != 5L) {
    stop("inappropriate dimension for input to VEplot_type2", call. = FALSE)
  }

  if (!(all(x[,4L] >= x[,2L]) || all(x[,4L] <= x[,2L])) ||
      !(all(x[,5L] >= x[,2L]) || all(x[,5L] <= x[,2L]))) {
    stop("confidence intervals do not have expected behaviors",
         call. = FALSE)
  }

  if (any(x[,2L] > 1.0)) {
    stop("vaccine efficacy > 100", call. = FALSE)
  }

  opts20 <- seq(-200, 200, 20)
  opts10 <- seq(-200, 200, 10)

  ymin <- floor(x = floor(x = min(x[,c(2L,4L,5L)]*100.0)) * 1.1)
  tst <- which(x = opts20 > ymin)[1L]
  if (tst > 1L) {
    ymin20 <- opts20[which(x = opts20 > ymin)[1L] - 1L]
    ymin10 <- opts10[which(x = opts10 > ymin)[1L] - 1L]
  } else {
    ymin20 <- ymin
    ymin10 <- ymin
  }

  ymax <- ceiling(x = ceiling(x = max(x[,c(2L,4L,5L)]*100.0)) * 1.1)
  ymax20 <- min(opts20[which(x = opts20 > ymax)[1L]], 100.0)
  ymax10 <- min(opts10[which(x = opts10 > ymax)[1L]], 100.0)

  if (ymin20 == ymin10 && ymax20 == ymax10) {
    ymin <- ymin20
    ymax <- ymax20
    dy <- 20L
  } else {
    ymin <- ymin10
    ymax <- ymax10
    dy <- 10L
  }

  graphics::plot(x = x[,1L], y = x[,2L]*100,
                 xlab = "Time Since Vaccination (in days)", 
                 ylab = ylab, 
                 type = 'l',
                 yaxt="n", main="", ylim = c(ymin, ymax))
  graphics::lines(x = x[,1L], y = x[,4L]*100, col = 3L)
  graphics::lines(x = x[,1L], y = x[,5L]*100, col = 3L)
  graphics::axis(side = 2L,
                 at = seq(ymin, ymax, by = dy),
                 labels = paste0(seq(ymin, ymax, by = dy), "%"), 
                 cex = 0.8)
  graphics::legend(x = "bottom",
                   legend = c(leg, "95% CI"), 
                   lty = c(1L,1L),
                   col = c(1L,3L), 
                   bg = "gray95",
                   cex = 0.75)

  return()
}
