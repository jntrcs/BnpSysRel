##Single Update demo


library(BnpSysRel)
library(gridExtra)
library(ggplot2)

set.seed(14)

gg_curve<-function (f, lower, upper, color="black")
{
  data = data.frame(x = seq(lower, upper, length.out = 200))
  data$y = f(data$x)
  return(geom_line(data = data, aes(x = x, y = y), color=color))
}

dDIFF <-  function(x, shape1, shape2=shape1, scale1=1, scale2=scale1) {
  res  <- x
  for (xx in seq(along=x)) {
    res[xx]  <-  integrate(function(y) dweibull(y+x[xx], shape1, scale1) * dweibull(y, shape2, scale2),  lower=0.0,  upper=+Inf)$value
  }
  return(res)
}



b=function(scale, prob){
  a= function(x)dDIFF(x, shape1=4, shape2=4, scale1=15, scale2=scale)
  abs(integrate(a, 0, 100)$value-prob)
}

scales=list(`30` = optimize(b, interval = c(0, 30), prob=.7)$minimum,
     `60` = optimize(b, interval = c(0, 30), prob=.40)$minimum,
     `90` = optimize(b, interval = c(0, 30), prob=.10)$minimum)


plots<-list()
scenarios=expand.grid(n=c(10,100,1000), Uncensored=c(30, 60, 90))

for (s in 1:nrow(scenarios)){
  n= scenarios$n[s]
  uncens =scenarios$Uncensored[s]

  startEndValve<-qweibull(c(.7,.99), 4,15)
  times<- rweibull(n, 4,15)
  cenTimes<-rweibull(n, 4, scales[[paste(uncens)]])

  cens = cenTimes< times


  times[cens]<-cenTimes[cens]
  data<-cbind(times, as.numeric(!cens))
  post=bspPosterior(createUninformativePrior(data), data)
  plots[[length(plots)+1]]=
    plot(post, withConfInt = T, withLegend = F)+
    gg_curve(f=function(x)pweibull(x,4,15), min(post$support), max(post$support), color="green")+
    ggtitle(paste0("n = ", n, " ", uncens, "% Uncensored"))+
    theme(plot.title = element_text(size=8))
}





a=arrangeGrob(grobs=plots, ncol=3, nrow=3,  widths=c(4,4,4))
plot(a)







