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
#' @param withConfInt Currently ignored. Future: If true, will calculate and plot approximate 95\% confidence intervals
#' @return A ggplot object
#' @export
#'
#' @examples
#' plot(bsp(c(1,5),c(0.2,0.4),c(2,1)))
#'
#'
plot.betaStacyProcess<- function(x, withConfInt=FALSE,
                                 withPaths=FALSE, conf.level=.95) {

  data=data.frame(times=x$support[-1], centeringMeasure=x$centeringMeasure[-1])
  p=ggplot2::ggplot(data,  ggplot2::aes(x=times, y=centeringMeasure))+
    ggplot2::geom_step(size=1.5)+ggplot2::geom_point()+
    ggplot2::xlab("Time (t)")+ggplot2::ylab("Centering Measure G(t)")
  #add some precision stuff
  if(withConfInt){
    samples=bspSampling(x, reps=1000)
    cred_int=data.frame(lower=apply(samples, 1, quantile, (1-conf.level)/2),
    upper =apply(samples, 1, quantile, 1-(1-conf.level)/2))
    p=p+ggplot2::geom_line(data=cred_int,size=1.5,
            ggplot2::aes(x=x$support[-1],y=lower, color="Credible Interval"))+
      ggplot2::geom_line(data=cred_int,size=1.5,
            ggplot2::aes(x=x$support[-1],y=upper, color="Credible Interval"))
    if (withPaths){
      paths=data.frame(samples)
      p= p+ ggplot2::geom_line(data=paths, aes(y=X1, x=x$support[-1]),
                        inherit.aes = F,
                        alpha=.2)
      for (i in 2:100){
        p=p+ggplot2::geom_line(data=NULL,
                               ggplot2::aes_(y=paths[[i]], x=x$support[-1]),
                              inherit.aes = F,
                              alpha=.2, show.legend = F)
        #plot(p)
      }
    }
  }
  return(p)
}

summary.betaStacyProcess<-function(x){
  print(paste("Median time to failure:",
              round(x$support[which.max(which(x$centeringMeasure<=.5))], 3)))
}

summary.bspPosteriorList<-function(posteriors){
  print(names(posteriors))
  print("This should print some more descriptive statistics")
}
