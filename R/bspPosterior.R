#' Update a BSP prior with data and return a BSP posterior
#'
#' @param bspPrior The BSP object that defines the BSP prior
#' @param data An nx2 matrix of data the first column are the actual data
#'             the second column contains the censoring variable
#'             0 - if right censored 1 -
#' @param precision The precision of the BSP
#'
#' @return An object representing a Beta-Stacy Process.
#' @export
#'
#' @examples
#' bsp(c(1,5),c(0.2,0.4),c(2,1))
#'
#'
bspPosterior <- function(bspPriorObject, data) {

  # Add checks to ensure the input is valid
  # check bspPriorObject is a bsp object
  # check data has observation times and censoring info (if not then censoring info is assumed??)


  # need a lot more here

  return(bspPriorObject)

}
