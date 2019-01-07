#' Define a Beta-Stacy Process
#'
#' @param support A numeric vector indicating the points at which the support for the BSP is defined
#' @param centeringMeasure The mean or centering measure of the BSP at each point on the support
#' or a function that can be evaluated at each point
#' @param precision A constant indicating the precision of the BSP (Currently, the precision
#' cannot vary over the support).
#'
#' @return An object representing a Beta-Stacy Process.
#' @export
#'
#' @details
#' The most common use case for this function will be defining BSPs to use as priors.
#' The precision reflects the uncertainty in the prior. In general, a one unit increase
#' in the precision is equivalent to one additional observation.
#'
#' @examples
#' bsp(support=c(1,5), centeringMeasure=c(0.2,0.4), precision=2)
#'
#'
#'
bsp <- function(support, centeringMeasure, precision, calculateMoments=FALSE) {

  # Add checks to ensure the input is valid
  # support is a numeric vector (sorted)
  # centeringMeasure is a vector of the same length as support all numbers must be in [0,1]
  #     and nondecreasing
  #     An alternate input is a CDF function that will be evaluated at the
  #     support points
  # precision is a nonnegative vector of numeric values and should be the same length as support
  #     if only one number is supplied it is replicated n times
  #    should these be nonincreasing???????
  if (is.function(centeringMeasure)){
    centeringMeasure<-centeringMeasure(support)
  }

  if (is.unsorted(support)) stop("support should be an increasing series of time points")

  if (length(centeringMeasure)!=length(support))stop("centeringMeasure and support length differ")
  if(any(centeringMeasure>1 |centeringMeasure<0))stop("All centeringMeasure points must be between 0 and 1")

  if (min(support)!=0){
    support<-c(0, support)
    centeringMeasure<-c(0, centeringMeasure)
  }else if (centeringMeasure[1]!=0)stop("Centering measure cannot be greater than 0 when t=0")
  if (!is.numeric(precision)| length(precision)!=1)
    stop("Precision must be a single number")
  precision<-rep(precision, length(support))

  #if (is.function(precision)){
   # precision<-precision(support)
    #if (length(precision)!=length(support))stop("Precision function must evaluate to same length as support")

  #}
  #if (length(precision)==1){
   # precision <- rep(precision, length(support))
  #}

  #if (length(precision)!=length(support))stop("precision and support length differ")

  bsp =makeBSP(support, centeringMeasure, precision)
  if (calculateMoments)
    bsp=E1E2(bsp)
  return(bsp)
}

#'@export
makeBSP<-function(support, centeringMeasure, precision){
  structure(list(support=support, centeringMeasure=centeringMeasure,
                 precision=precision), class="betaStacyProcess")
}

#' Evaluate centering measure at specific times
#'
#' @param bsp A BSP object
#' @param times A list of times that the centering measure should be determined for
#'
#' @return A vector of same length as times with the centering measure for each time
#' @export
#'
#'@note For times in between jumps on the support, the centering measure is
#'considered equal to its value after the last jump
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
#' evaluate_precision(bsp(c(1,2), c(.2,.6), 1), c(.5,1.5,2.5))

evaluate_precision<-function(bsp, times){
  support<-bsp$support
  support[support==0]=-.1
  precision<-bsp$precision


  # sapply(times, FUN=function(time){
  #   if (time<min(support)){warning("Precision not specfied for time < min(support), assumed to be 0")
  #     return(0.0001)}else{
  #     index<-max(which(time>=support))
  #     return(ifelse(time==support[index], precisionAt[index], precisionAfter[index]))
  #   }
  # })
}

