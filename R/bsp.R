#' Define a Beta-Stacy Process
#'
#' @param support The support where the BSP is defined.
#' @param centeringMeasure The mean or centering measure of the BSP
#' @param precision The precision of the BSP
#'
#' @return An object representing a Beta-Stacy Process.
#' @export
#'
#' @examples
#' bsp(c(1,5),c(0.2,0.4),c(2,1))
#'
#'
bsp <- function(support, centeringMeasure, precision) {

  # Add checks to ensure the input is valid
  # support is a numeric vector (sorted)
  # centeringMeasure is a vector of the same length as support all numbers must be in [0,1]
  #     and nondecreasing
  #     An alternate input is a CDF function that will be evaluated at the
  #     support points
  # precision is a nonnegative vector of numeric values and should be the same length as support
  #     if only one number is supplied it is replicated n times
  #    should these be nonincreasing???????

  structure(list(support=support, centeringMeasure=centeringMeasure,
                 precision=precision), class="betaStacyProcess")
}
