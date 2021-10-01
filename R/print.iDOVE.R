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
#' print(x = result)
#' 
# 8/11/21: removed lines regarding additional attributes

print.iDOVE <- function(x, ...) {

  print(x = unclass(x = x))

  return( NULL )

}
