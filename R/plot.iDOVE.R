#' Plot Estimated Vaccine Efficacy
#'
#' Generates plots of the estimated vaccine efficacy in reducing attack rate, 
#'   the estimated vaccine efficacy in reducing the hazard rate, 
#'   and their 95\% confidence intervals.
#'
#' @param x An iDOVE object. The value object returned by idove().
#'
#' @param ... ignored
#' 
#' @method plot iDOVE
#'
#' @returns No return value, called to produce graphical elements.
#'
#' @examples
#'
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
#' result <- idove(formula = intCens(entry.time, left.time, right.time, vaccine.time) ~ 1, 
#'                 data = idoveData[smp,])
#' 
#' plot(x = result)
#' 
#' @name plot
#' @export
#' @importFrom graphics plot
#'
plot.iDOVE <- function(x, ...) {
  
  # 8/11/21: replace the original arguments of knots, gamma, covgamma, and tau
  # to VE_a and VE_h.

  .VEplot(VE_a = x$vaccine$VE_a, VE_h = x$vaccine$VE_h)

  return( NULL )

}
