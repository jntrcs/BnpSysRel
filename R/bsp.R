#' Define a Beta-Stacy Process
#'
#' @param support The support where the BSP is defined.
#' @param centeringMeasure The mean or centering measure of the BSP at each point on the support
#' or a function that can be evaluated at each point
#' @param precision The precision can be one of: a scalar to indicate constant precision,
#'  a vector of equal length to the support indicating precision at each time,
#'  a function that will be evaluated at each time of the support to obtain the precision
#'
#' @return An object representing a Beta-Stacy Process.
#' @export
#'
#' @examples
#' bsp(support=c(1,5), centeringMeasure=c(0.2,0.4), precision=c(2,1))
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
  if (is.function(precision)){
    precision<-precision(support)
    if (length(precision)!=length(support))stop("Precision function must evaluate to same length as support")

  }
  if (length(precision)==1){
    precision <- rep(precision, length(support))
  }

  precisionAfter=precision

  #if (length(precision)!=length(support))stop("precision and support length differ")

  #centeringMeasure<-c(0, centeringMeasure)
  #support<-c(0, support)
  #includePoint=rep(c(T,F), 100)[1:length(support)]
  makeBSP(support, centeringMeasure, precision, precisionAfter)

}

makeBSP<-function(support, centeringMeasure, precision, precisionAfter){
  structure(list(support=support, centeringMeasure=centeringMeasure,
                 precisionAt=precision, precisionAfter=precisionAfter,
                 class="betaStacyProcess"))
}

#' Evaluate centering measure at specific times
#'
#' @param bsp A BSP object
#' @param times A list of times that the centering measure should be determined for
#'
#' @return A vector of same length as times with the centering measure for each time
#' @export
#'
#'@note The centering measure is considered constant
#' @examples
#' evaluate_centering_measure(bsp(c(1,2), c(.2,.6), 1), c(.5,1.5,2.5))

evaluate_centering_measures<-function(bsp, times){
  support<-bsp$support
  centeringMeasure<-bsp$centeringMeasure
  indices<-sapply(times, FUN=function(x)sum(support<=x))
  cMeas<-centeringMeasure[indices]
  indices[indices!=0]<-cMeas
  indices
}

#' Evaluate precision at specific times
#'
#' @param bsp A BSP object
#' @param times A list of times that the centering measure should be determined for
#'
#' @return A vector of same length as times with the centering measure for each time
#' @export
#'
#' @examples
#' evaluate_centering_measure(bsp(c(1,2), c(.2,.6), 1), c(.5,1.5,2.5))

evaluate_precision<-function(bsp, times){
  support<-bsp$support
  precisionAt<-bsp$precisionAt
  precisionAfter<-bsp$precisionAfter
  jumpDown<-bsp$jumpDown
  sapply(times, FUN=function(time){
    if (time<min(support)){warning("Precision not specfied for time < min(support), assumed to be 0")
      return(0.0001)}else{
      index<-max(which(time>=support))
      return(ifelse(time==support[index], precisionAt[index], precisionAfter[index]))
    }
  })
}

