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
  res=data.frame(lapply(x[1:3], FUN=as.numeric))
  print(res)
}

#' Plots a Beta-Stacy process over its support
#'
#' @param x The bsp object to be plotted
#' @param withConfInt If true, will calculate and plot approximate 95\% confidence intervals, default FALSE
#' @param withPaths If true, will plot 100 sampled CDFs in gray (must also set withConfInt=T), default FALSE
#' @param conf.level Confident level for credible interval, default 0.95
#' @return A ggplot object
#' @export
#'
#' @examples
#' plot(bsp(c(1,5),c(0.2,0.4),c(2,1)))
#'
#'
plot.betaStacyProcess<- function(x, withConfInt=FALSE,
                                 withPaths=FALSE, conf.level=.95) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"pkg\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (withPaths &!withConfInt)warning("withPaths ignored if withConfInt set to false")
  if (conf.level<=0 | conf.level>=1)stop("Confidence must be in (0, 1)")

  data=data.frame(times=x$support[-1], centeringMeasure=x$centeringMeasure[-1])
  p=ggplot2::ggplot(data,  ggplot2::aes(x=times, y=centeringMeasure))+
    ggplot2::geom_line(size=1.5)+ggplot2::geom_point()+
    ggplot2::xlab("Time (t)")+ggplot2::ylab("Centering Measure G(t)")
  #add some precision stuff
  if(withConfInt){
    samples=bspSampling(x, reps=1000)
    cred_int=data.frame(lower=apply(samples, 1, quantile, (1-conf.level)/2, na.rm=T),
    upper =apply(samples, 1, quantile, 1-(1-conf.level)/2, na.rm=T))
    p=p+ggplot2::geom_line(data=cred_int,size=1.5,
            ggplot2::aes(x=x$support[-1],y=lower, color="Credible Interval"))+
      ggplot2::geom_line(data=cred_int,size=1.5,
            ggplot2::aes(x=x$support[-1],y=upper, color="Credible Interval"))
    if (withPaths){
      paths=data.frame(samples)
      p= p+ ggplot2::geom_line(data=paths, ggplot2::aes(y=X1, x=x$support[-1]),
                        inherit.aes = F,
                        alpha=.2)
      for (i in 2:100){
        p=p+ggplot2::geom_line(data=NULL,
                               ggplot2::aes_(y=paths[[i]], x=x$support[-1]),
                              inherit.aes = F,
                              alpha=.2, show.legend = F)
      }
    }
  }
  plot(p)
  return(invisible(p))
}


#' Calculates the expected time to failure for various quantiles
#'
#'@param bsp The bsp object to calculate from
#'@param probs numeric vector of probabilities
#' @return a named numeric vector of times associated with each quantile
#' @export
#'
#'@details This function calculates the expected failure time for a given quantile from
#'the centering measure of the BSP. Linear interpolation is to choose time between two discrete points
#'
#' @examples
#' bsp=bsp(1:3, c(.25,.5,.75), 1)
#' quantile(bsp, c(.333, .666))

quantile.betaStacyProcess<-function(bsp, probs){
  indices<-sapply(probs, function(p)sum(bsp$centeringMeasure<=p))
  #qs<-sapply(probs, function(p)bsp$support[sum(c(bsp$centeringMeasure,1)<=p)])
  p1=bsp$centeringMeasure[indices]
  p2= c(bsp$centeringMeasure,1)[indices+1]

  weights <- (probs-p1)/(p2-p1)
  qs=(1-weights)*bsp$support[indices] + weights*c(bsp$support,1)[indices+1]

  names(qs)<-paste0(probs*100,"%")
  qs


}

summary.betaStacyProcess<-function(x){
  print(paste("Median time to failure:",
              round(x$support[which.max(which(x$centeringMeasure<=.5))], 3)))
}


