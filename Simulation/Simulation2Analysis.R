#Simulation2Analysis.R
library(xtable)
library(BnpSysRel)
load("./Simulation/Simulation2Results.RData")

coverageArray<-array(data=0, dim=c(5, 3, 3, 2),
                     dimnames = list(Parts=c("valve", "generator", "plumbing",
                                             "power", "system"),
                                     Tails=c("Left", "Middle", "Right"),
                                     N=c(10, 50, 250), Censored=c("FALSE", "TRUE")))



biasArray<-array(data=0, dim=c(5, 3, 3, 2),
                 dimnames = list(Parts=c("valve", "generator", "plumbing",
                                         "power", "system"),
                                 Tails=c("Left", "Middle", "Right"),
                                 N=c(10, 50, 250), Censored=c("FALSE", "TRUE")))
isBiasedArray<-array(data=0, dim=c(5,3,3,2),
                     dimnames = list(Parts=c("valve", "generator", "plumbing",
                                             "power", "system"),
                                     Tails=c("Left", "Middle", "Right"),
                                     N=c(10, 50, 250), Censored=c("FALSE", "TRUE")))

extractor<-function(results, measurement, n, cens, part="all", time="all"){
  n<-as.character(n)
  cens<-as.character(cens)
  part_names<-c(paste0("valve", 1:3), paste0("generator",1:2), "plumbing", "power", "system")
  time_names<-c("first", "middle", "last")
  if (part=="all") columns=1:5 else columns = which(gsub("^|\\d+$", "", part)==gsub("^|\\d+$", "", part_names))
  if (time=="all")rows=1:3 else rows= which(time==time_names)
  sapply(results[[n]][[cens]], FUN=function(x)x[[measurement]][rows, columns])
}

ns<-c(10,50,250)
censoring=c(T, F)
#par(ask=T)
for (n in ns){
  for (c in censoring){
    for (p in c("valve", "generator", "plumbing",
                "power", "system")){
      for (h in c("first", "middle", "last")){
        n<-as.character(n)
        c<-as.character(c)

        cov1=extractor(results, "Coverage", n, c, part=p, time=h )
        bias=extractor(results, "Bias", n, c, part=p, time=h )

        h=c("Left", "Middle", "Right")[which(h==c("first", "middle", "last"))]
        #print(paste( p, h, n,c))
        coverageArray[p, h, n, c]<-mean(cov1)
        biasArray[p, h, n, c]<-mean(bias)
        mcInt<-mean(bias)+qnorm(c(.005, .995))*sqrt(var(c(bias))/length(bias))
        isBiasedArray[p, h, n, c]<-!( mcInt[1]<0 & mcInt[2]>0)
      }
    }

  }
}


library(reshape2)
library(ggplot2)
meltedCoverage=melt(coverageArray)
meltedCoverage$Significant = abs(meltedCoverage$value-.95)>qnorm(.975)*sqrt(.95*.05/3000)
ggplot(meltedCoverage, aes(y=value, x=N, color=Parts, linetype=Censored))+ facet_grid(. ~ Tails)+
  geom_line()+#geom_point(aes(shape=Significant))+
  geom_hline(aes(yintercept = .95))+
  scale_x_continuous(breaks=c(10,50, 250),trans="log2")+
  ylab("Coverage")+ggtitle("Coverage for 95% Credible Intervals")#+
  #ggsave("Simulation/Simulation2Coverage.pdf", width=5, height=6, units="in")


meltedBias=melt(biasArray)
meltedIsBiased=melt(isBiasedArray)

meltedBias$Significant= meltedIsBiased$value
ggplot(meltedBias, aes(y=value, x=N, color=Parts,
                       linetype=Censored))+ facet_grid(. ~ Tails)+
  geom_line()+geom_point(size=2,aes(shape=factor(Significant)))+
  ylab("Bias")+geom_hline(yintercept = 0)+ggtitle("Average Bias")+
  scale_x_continuous(breaks=c(10,50, 250),trans="log2")+
  guides(shape=guide_legend(title="MC 95% CI"))+
  scale_shape_discrete(labels=c("Contains 0", "Does not"))#+
  #ggsave("Simulation/Simulation2Bias.pdf", width=5, height=6, units="in")


