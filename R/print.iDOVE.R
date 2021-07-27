#' Print the Primary Results of an idove() Analysis
#'
#' Print the primary results of an idove() analysis.
#'
#' @param x An iDOVE object. The value object returned by a call to idove()
#'
#' @param ... ignored
#'
#' @export
#' @name print
#' @method print iDOVE
#'
#' @returns No return value, called to display key results.
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
#' print(x = result)
#' 
print.iDOVE <- function(x, ...) {

  attr(x = x, which = "knots") <- NULL
  attr(x = x, which = "gamma") <- NULL
  attr(x = x, which = "covgamma") <- NULL
  attr(x = x, which = "tau") <- NULL

  x <- unclass(x = x)

  print(x = x)

}
