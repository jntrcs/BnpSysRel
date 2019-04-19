#Sim3ResultsDisplay
library(xtable)
library(BnpSysRel)
load("./Simulation/Simulation3.RData")

set.seed(14)
valves<-replicate(10,rweibull(10000, 4,20))
generators<-replicate(10, rchisq(10000, 3))
plumbing<-apply(valves, 1, min)
power<-apply(generators, 1, max)
system<-pmin(plumbing, power)

ns<-c(10,50)
censoring<-c(T, F) #100% of the data and 80% of the data will be uncensored

evaluation_times<-list()
evaluation_times$valve<-qweibull(c(.2,.5,.8), 4, 20)
evaluation_times$generator<-qchisq(c(.2,.5,.8),3)
evaluation_times$plumbing<-quantile(plumbing, c(.2,.5,.8))
evaluation_times$power<-quantile(power, c(.2,.5,.8))
evaluation_times$system<-quantile(system, c(.2,.5,.8))

coverageArray<-array(data=0, dim=c(5, 3, 2, 2),
                     dimnames = list(Parts=names(evaluation_times), Tails=c("Left", "Middle", "Right"),
                                     N=c(10, 50), Censored=c("FALSE", "TRUE")))

coverage2Array<-array(data=0, dim=c(5, 3, 2, 2),
                      dimnames = list(Parts=names(evaluation_times), Tails=c("Left", "Middle", "Right"),
                                      N=c(10, 50), Censored=c("FALSE", "TRUE")))

biasArray<-array(data=0, dim=c(5, 3, 2, 2),
                 dimnames = list(Parts=names(evaluation_times), Tails=c("Left", "Middle", "Right"),
                                 N=c(10, 50), Censored=c("FALSE", "TRUE")))
isBiasedArray<-array(data=0, dim=c(5, 3, 2, 2),
                     dimnames = list(Parts=names(evaluation_times), Tails=c("Left", "Middle", "Right"),
                                     N=c(10, 50), Censored=c("FALSE", "TRUE")))
#coverageArray

extractor<-function(results, measurement, n, cens, part="all", time="all"){
  n<-as.character(n)
  cens<-as.character(cens)
  part_names<-c(paste0("valve", 1:10), paste0("generator",1:10), "plumbing", "power", "system")
  time_names<-c("first", "middle", "last")
  if (part=="all") columns=1:23 else columns = which(gsub("^|\\d+$", "", part)==gsub("^|\\d+$", "", part_names))
  if (time=="all")rows=1:3 else rows= which(time==time_names)
  sapply(results[[n]][[cens]], FUN=function(x)x[[measurement]][rows, columns])
}


#par(ask=T)
for (n in ns){
  for (c in censoring){
    for (p in names(evaluation_times)){
      for (h in c("first", "middle", "last")){
        n<-as.character(n)
        c<-as.character(c)

        cov1=extractor(results, "Coverage", n, c, part=p, time=h )
        cov2=extractor(results, "Coverage2", n, c, part=p, time=h )
        bias=extractor(results, "Bias", n, c, part=p, time=h )

        cov2[is.na(cov2)]<-0
        if (length(bias)!=100)print(n)
        h=c("Left", "Middle", "Right")[which(h==c("first", "middle", "last"))]
        #print(paste( p, h, n,c))
        coverageArray[p, h, n, c]<-mean(cov1)
        coverage2Array[p, h, n, c]<-mean(cov2)
        biasArray[p, h, n, c]<-mean(bias)
        mcInt<-mean(bias)+qnorm(c(.005, .995))*sqrt(var(bias)/length(bias))
        isBiasedArray[p, h, n, c]<-!( mcInt[1]<0 & mcInt[2]>0)
      }
    }

  }
}

isBiasedArray

apply(coverageArray, 4, function(x)apply(x, 3, xtable))

apply(coverageArray, c(3,4), xtable)

a=extractor(results, 'Bias',50, T, 'plumbing', 'last')
mean(a)+qnorm(c(.025,.975))*sqrt(var(a)/length(a))
hist(a)

require(reshape2)
require(ggplot2)
meltedCoverage=melt(coverageArray)
ggplot(meltedCoverage, aes(y=value, x=N, color=Parts, linetype=Censored))+ facet_grid(. ~ Tails)+
  geom_line()+geom_point()+
  geom_hline(aes(yintercept = .95))+
  ylab("Coverage")+ggtitle("Coverage for 95% Credible Intervals")


meltedBias=melt(biasArray)
meltedBias=meltedBias[meltedBias$Parts%in% c("power", "plumbing", "system"),]
meltedIsBiased=melt(isBiasedArray)
meltedIsBiased=meltedIsBiased[meltedIsBiased$Parts%in% c("power", "plumbing", "system"),]

meltedBias$Significant= meltedIsBiased$value
ggplot(meltedBias, aes(y=value, x=N, color=Parts,
                       linetype=Censored))+ facet_grid(. ~ Tails)+
  geom_line()+geom_point(size=2,aes(shape=factor(Significant)))+
  ylab("Bias")+geom_hline(yintercept = 0)+ggtitle("Average Bias")+
  guides(shape=guide_legend(title="MC 95% CI"))+scale_shape_discrete(labels=c("Contains 0", "Does not"))#+
  #ggsave("Simulation/Sim3Bias.pdf")


# ##Appending on the bias results from justParts (fixing the error in Sim2Results.Rdata)
# load("JustPartResults.RData")
# biasArrayJP<-array(data=0, dim=c(2, 3, length(ns), length(censoring)),
#                    dimnames = list(Parts=c("valve","generator"), Tails=c("Left", "Middle", "Right"),
#                                    N=c(10, 50, 250), Censored=c("FALSE", "TRUE")))
# isBiasedArrayJP<-array(data=0, dim=c(2, 3, length(ns), length(censoring)),
#                        dimnames = list(Parts=c("valve","generator"), Tails=c("Left", "Middle", "Right"),
#                                        N=c(10, 50, 250), Censored=c("FALSE", "TRUE")))
#
# extractor<-function(results, measurement, n, cens, part="all", time="all"){
#   n<-as.character(n)
#   cens<-as.character(cens)
#   part_names<-c("valve","generator")
#   time_names<-c("first", "middle", "last")
#   if (part=="all") columns=1:5 else columns = which(gsub("^|\\d+$", "", part)==gsub("^|\\d+$", "", part_names))
#   if (time=="all")rows=1:3 else rows= which(time==time_names)
#   sapply(results[[n]][[cens]], FUN=function(x)x[[measurement]][rows, columns])
# }

