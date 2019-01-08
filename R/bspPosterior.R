#' Update a BSP prior with data and return a BSP posterior
#'
#' @param bspPrior The BSP object that defines the BSP prior
#' @param data An nx2 matrix of data the first column are the actual data
#'             the second column contains the censoring variable
#'             0 - if right censored 1 - if fully observed
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

  ###CHECKS
  if (!is.matrix(data))stop("data must be an nx2 matrix")
  #Should we have a "check valid BSP" function?

  ####ARRANGING DATA
  data<-data[order(data[,1]),, drop=F]
  censor<-data[,2]
  data<-data[,1]
  prior<-bspPriorObject

  #####COMPUTATION
  prior_ts<-prior$support
  prior_gs<-prior$centeringMeasure
  prior_precs<-prior$precision
  if (any(data!=sort(data)))stop("Data must be sorted (increasing)")
  m<-function(data, t) sapply(t, FUN=function(x) sum(data>=x))
  j<-function(data, censor, t) sapply(t, FUN=function(x) sum(censor*(data==x)))
  prior_g<-function(t)sapply(t, FUN=function(t)c(0,prior_gs)[sum(prior_ts<=t)+1])
  allJumps<-sort(unique(c(prior_ts, data)))
  precAtJumps<-evaluate_precision(prior, allJumps)
  GAtJumps<-prior_g(allJumps)
  GBeforeJumps<-c(0, GAtJumps[-length(GAtJumps)])
  allJs<-j(data, censor, allJumps)
  allMs<-m(data, allJumps)

  cumulativeProduct = cumprod(1-
                                ((precAtJumps)*(GAtJumps-GBeforeJumps)+allJs)/
                                (precAtJumps*(1-GBeforeJumps)+allMs)
  )




  centeringMeasure = 1-cumulativeProduct
  isna=F
  alpha<- (precAtJumps*(1-GAtJumps)+allMs-allJs)/(1-centeringMeasure)
  if(any(is.na(alpha))){
    isna=T
    first.na=which(is.na(alpha))[1]
    newPoint<-mean(allJumps[(first.na-1):first.na])

    new.alpha<-(evaluate_precision(prior, newPoint)*(1-prior_g(newPoint))+
                  m(data, newPoint)-
                  j(data, censor, newPoint))/
      (1-centeringMeasure[first.na-1])
  }
  lastTime=allJumps[length(allJumps)]+1
  lastAlpha<-(evaluate_precision(prior, lastTime)*(1-prior_gs[length(prior_gs)]))/#+m(data, lastTime)-j(data, censor, lastTime))/
    (1-centeringMeasure[length(centeringMeasure)])
  #measuresAfter<-calc(support+.00001)
  #precisionAfter<-measuresAfter$Precisions
  alphaAfter<-c(alpha[-1], lastAlpha)
  if(isna) alphaAfter[first.na-1]<-new.alpha
  return(makeBSP(allJumps, centeringMeasure, alphaAfter))


}


#'Creates a matrix where each column represents a sampled CDF from the provided BSP
#'
#'@param bsp Any bsp object
#'@param reps Number of samples to draw, default 10000
#'
#'@return A matrix. Each column can be plotted against the BSP's support to view a sampled CDF from the BSP
#'@export
#'
bsp_sampling<-function(bsp, reps=10000){
  draws<-matrix(0,ncol = reps, nrow=length(bsp$support))
  firstTime<-bsp$support[1]
  alpha1=evaluate_precision(bsp, firstTime)*bsp$centeringMeasure[1]
  beta1=evaluate_precision(bsp, firstTime)*(1-bsp$centeringMeasure[1])
  draws[1,]<-rbeta(reps, alpha1, beta1)
  for (i in 2:length(bsp$support)){
    prec<-evaluate_precision(bsp, bsp$support[i])
    alphak<-prec*(bsp$centeringMeasure[i]-bsp$centeringMeasure[i-1])
    betak=prec*(1-(bsp$centeringMeasure[i]))
    epsk<-1-draws[i-1,]
    #if (any(epsk>1 | epsk<0))print(i)
    draws[i,]<-draws[i-1,]+rbeta(reps, alphak, betak)*epsk
  }
  draws
}
