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
#' smp <- sample(1L:nrow(x = idoveData), size = 500L)
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

  .VEplot(knots = attr(x = x, which = "knots"),
          tau = attr(x = x, which = "tau"),
          gamma = attr(x = x, which = "gamma"),
          covgamma = attr(x = x, which = "covgamma"))

}
