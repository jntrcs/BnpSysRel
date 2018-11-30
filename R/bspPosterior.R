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


  #####COMPUTATION
  make_calculator<-function(prior_alpha, prior_ts, prior_gs, data, censor){
    if (any(data!=sort(data)))stop("Data must be sorted (increasing)")
    m<-function(data, t) sapply(t, FUN=function(x) sum(data>=x))
    j<-function(data, censor, t) sapply(t, FUN=function(x) sum(censor*(data==x)))
    prior_g<-function(t)sapply(t, FUN=function(t)c(0,prior_gs)[sum(prior_ts<=t)+1])
    single_t<-function(t){
      previous<-unique(data[data<=t & censor])
      previous<-sort(unique(c(previous, prior_ts[prior_ts<=t])))
      if(length(previous)==0) {
        previous=0
      }
      prior_at_prev<-prior_alpha(previous)
      g = 1-prod(1-(prior_at_prev*(prior_g(previous)-prior_g(previous-.00001)) +j(data, censor, previous))/
                   (prior_at_prev*(1-prior_g(previous-.00001))+m(data, previous)))
      alpha = (prior_alpha(t)*(1-prior_g(t))+m(data, t)-j(data, censor, t))/(1-g)

      return(list(CenteringMeasure=g, Alpha=alpha))
    }
    function(ts){
      a=sapply(ts, FUN=single_t)
      return(list(centeringMeasures=as.numeric(a[1,]), Precisions=as.numeric(a[2,]), support=ts))
    }
  }
  calc<-make_calculator(function(x)evaluate_precision(bspPriorObject, x), bspPriorObject$support,
                        bspPriorObject$centeringMeasure, data[,1], data[,2])


  support = unique(sort(c(bspPriorObject$support, data[,1])))
  newMeasures<-calc(support)
  centeringMeasure = newMeasures$centeringMeasures
  alpha<-newMeasures$Precisions
  precisionAfter<-alpha
  measuresAfter<-calc(support+.00001)
  precisionAfter<-measuresAfter$Precisions
  return(makeBSP(support, centeringMeasure, alpha,precisionAfter))

}


#'Creates a matrix where each column represents a sampled CDF from the provided BSP
#'
#'@param BSP Any BSP object
#'@param reps Number of samples to draw, default 10000
#'
#'@return A matrix. Each column can be plotted against the BSP's support to view a sampled CDF from the BSP
#'@export
#'
bsp_sampling<-function(BSP, reps=10000){
  draws<-matrix(0,ncol = reps, nrow=length(BSP$support))
  firstTime<-BSP$support[1]
  alpha1=BSP$evaluate_alpha(firstTime)*BSP$centeringMeasure[1]
  beta1=BSP$evaluate_alpha(firstTime)*(1-BSP$centeringMeasure[1])
  draws[1,]<-rbeta(reps, alpha1, beta1)
  for (i in 2:length(BSP$support)){
    prec<-BSP$evaluate_alpha(BSP$support[i])
    alphak<-prec*(BSP$centeringMeasure[i]-BSP$centeringMeasure[i-1])
    betak=prec*(1-(BSP$centeringMeasure[i]))
    epsk<-1-draws[i-1,]
    #if (any(epsk>1 | epsk<0))print(i)
    draws[i,]<-draws[i-1,]+rbeta(reps, alphak, betak)*epsk
  }
  draws
}
