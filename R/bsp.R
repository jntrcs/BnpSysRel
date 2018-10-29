#' Define a Beta-Stacy Process
#'
#' @param support The support where the BSP is defined.
#' @param centeringMeasure The mean or centering measure of the BSP
#' @param precision A function that can be evaluated at any value of the support to obtain the precision, or a scalar to indicate constant precision
#' @note If you want precision to change as a stepwise function, it's recommended you construct the function using ifelse()
#'
#' @return An object representing a Beta-Stacy Process.
#' @export
#'
#' @examples
#' bsp(c(1,5),c(0.2,0.4),c(2,1))
#'
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

  if (any(!support==sort(support))) stop("support should be an increasing series of time points")
  if (is.function(centeringMeasure)){
    centeringMeasure<-centeringMeasure(support)
  }
  if (length(centeringMeasure)!=length(support))stop("centeringMeasure and support length differ")
  if(any(centeringMeasure>1 |centeringMeasure<0))stop("All centeringMeasure points must be between 0 and 1")
  #print(precision)
  if (!is.function(precision))
  if(is.numeric(precision) & length(precision)==1){
    precisionNum<-precision
    precision<-function(x) rep(precisionNum, length(x))
  }
  else stop("Precision must be single scalar or function")

  if (length(precision(support))!=length(support))stop("Precision function must return vector of precisions same length as t")
  #if (length(precision)!=length(support))stop("precision and support length differ")

  #centeringMeasure<-c(0, centeringMeasure)
  #support<-c(0, support)

  evaluate_center<-function(t){
    indices<-sapply(t, FUN=function(x)sum(support<=x))
    cMeas<-centeringMeasure[indices]
    indices[indices!=0]<-cMeas
    indices
  }

  structure(list(support=support, centeringMeasure=centeringMeasure,
                 evaluate_alpha= precision,
                 evaluate_g=evaluate_center), class="betaStacyProcess")
}

