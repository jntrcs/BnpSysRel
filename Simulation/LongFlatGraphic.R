
library(BnpSysRel)

set.seed(104)
data=rweibull(25, 4, 20)
censored<-1-rbinom(25, 1, prob=data/sum(data))
censored[data>24 & data<28]=0

a=bspPosterior(createUninformativePrior(cbind(data,censored)), cbind(data, censored))
plot(a)+geom_point(aes(x=a$support[-1], y=a$centeringMeasure[-1]), size=3.5, shape=19)#+
 # ggsave("Simulation/LongFlat2.jpg")
