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
  times<-data[,1]
  C = data[,2] #C as used in the paper

  get_base<-function(BSP, t){
    a=which(t>=BSP$support)
    if(length(a)==0)return(0)
    return(BSP$centeringMeasure[max(a)])
  }

  get_precision<-function(BSP, t){
    BSP$precision[max(which(t>=BSP$support))]
  }

  ts<-sort(unique(c(bspPriorObject$support,times)))
  M<-function(t){
    sum(times>=t)
  }
  J<-function(t){
    if (t==0)return(0) #What should this be?
    return(sum(as.numeric(times*C==t)))
  }
  runningProd=1
  centeringMeasure=rep(0, length(ts))
  precision<-rep(0, length(ts))
  for (i in 1:length(ts)){
    t<-ts[i]
    alpha_T<-get_precision(bspPriorObject,t)
    tMinusOne = ifelse(i==1, 0, ts[i-1])
    centeringMeasure_T=get_base(bspPriorObject, t)
    centeringMeasureTMinusOne=get_base(bspPriorObject, tMinusOne)
    M_T<-M(t)
    J_T<-J(t)
    runningProd = runningProd* (1-(alpha_T*(centeringMeasure_T- centeringMeasureTMinusOne)+J_T)/
                                  (alpha_T*(1-centeringMeasureTMinusOne)+M_T))
    centeringMeasure[i]<-1-runningProd
    precision[i]<- (alpha_T*(1-centeringMeasure_T) + M_T -J_T)/(1-centeringMeasure[i])
  }
  return(bsp(support = ts, centeringMeasure = centeringMeasure, precision = precision))

}
