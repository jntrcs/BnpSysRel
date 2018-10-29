#' Prints a Beta-Stacy process object as a data frame
#'
#'@param x The bsp object to be printed
#' @return none
#' @export
#'
#' @examples
#' print(bsp(c(1,5),c(0.2,0.4),c(2,1)))
#'
#'
print.betaStacyProcess<-function(x){
  res=data.frame(lapply(x[1:2], FUN=as.numeric))
  res$alpha<-x$evaluate_alpha(res$support)
  res
}

#' Plots a Beta-Stacy process over its support
#'
#' @param x The bsp object to be plotted
#' @param withConfInt Currently ignored. Future: If true, will calculate and plot approximate 95\% confidence intervals
#' @return A ggplot object
#' @export
#'
#' @examples
#' plot(bsp(c(1,5),c(0.2,0.4),c(2,1)))
#'
#'
plot.betaStacyProcess<- function(x, withConfInt=FALSE) {

  data=data.frame(times=x$support, centeringMeasure=x$centeringMeasure)
  ggplot2::ggplot(data, ggplot2::aes(x=times, y=centeringMeasure))+ggplot2::geom_step()+ggplot2::geom_point()+
    ggplot2::xlab("Time (t)")+ggplot2::ylab("Centering Measure G(t)") #add some precision stuff
}
