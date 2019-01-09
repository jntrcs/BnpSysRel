##Simulation questions:

#1. Are the posterior estimates asymptotically consistent in complex systems with
#both parallel and series components

#steps:
# a. Create reliability diagrams for the system we want to test
# b. Define each components true DGD (data generating distribution, eg. rnorm(5,4))
# c. Define priors for each component
# d. Generate data using DGD
# e. Calculate posterior
# f. Compare estimates of F(t) in the middle and tails
# g. Track variance of estimates for various levels of data.

require(BnpSysRel)

dataList<-list()

#valve ~ Weibull(4, 10)
curve(dweibull(x, 4, 20),  xlim=c(0,30))
set.seed(100)
valve_data<-rweibull(30, 4, 20)
censored<-valve_data>25
valve_data[censored]<-25
dataList$valve<-cbind(valve_data, as.numeric(!censored))

#Generator~chisq(3)
curve(dchisq(x, 3), xlim=c(0,30))
gen_data<-rchisq(15, 3)
g_censored<-gen_data>7
gen_data[g_censored]<-7
dataList$generator<-cbind(gen_data, as.numeric(!g_censored))




posteriors=estimateSystemReliability("Simulation/System1.txt", priorList =list(),
                          dataList = dataList)

valves<-matrix(rweibull(5000, 4, 20), ncol=5, byrow=T)
plumbing<-apply(valves, 1, min)
generators<-matrix(rchisq(2000, 3), ncol=2, byrow=T)
power<-apply(generators, 1, max)

true_system<-apply(cbind(plumbing, power),1, min)
hist(power)
mean(true_system==plumbing)

gg_curve<-function(f, lower, upper){
  data=data.frame(x=seq(lower, upper, length.out=200))
  data$y=f(data$x)
  return(geom_line(data=data, aes(x=x, y=y)))
}

plot(posteriors$valve)+gg_curve(function(x)pweibull(x, 4,20), 12,25)

plot(posteriors$generator)+gg_curve(function(x)pchisq(x, 3), 0,8)


