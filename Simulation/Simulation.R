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

##Calculate the truth
valves<-replicate(5,rweibull(100000, 4,20))
generators<-replicate(2, rchisq(100000, 3))
plumbing<-apply(valves, 1, min)
power<-apply(generators, 1, max)
system<-pmin(plumbing, power)

ns<-c(50)
censoring<-c(T) #100% of the data and 80% of the data will be uncensored

evaluation_times<-list()
evaluation_times$valve<-c(10, 20, 24)
evaluation_times$generator<-c(1, 4, 6)
evaluation_times$plumbing<-c(7, 13,16)
evaluation_times$power<-c(2,4,7)
evaluation_times$system<-c(2, 4, 7)

allTimes<-sort(unique(unlist(evaluation_times)))

valve_prior<-bsp(evaluation_times$valve, pweibull(evaluation_times$valve, 4,20), .2)
generator_prior<-bsp(evaluation_times$generator, pchisq(evaluation_times$generator,3), .2)

priorList=list(valve=valve_prior, generator= generator_prior)
#plot(generator_prior, withConfInt = T, withPaths=T)+gg_ecdf(generators[,1])

true_failure_rates<-list()
true_failure_rates$valve<-pweibull(evaluation_times$valve, 4, 20)
true_failure_rates$generator<-pchisq(evaluation_times$generator, 3)
true_failure_rates$plumbing<-sapply(evaluation_times$plumbing, FUN=function(x)mean(plumbing<=x))
true_failure_rates$power<-sapply(evaluation_times$power, FUN=function(x)mean(power<=x))
true_failure_rates$system<-sapply(evaluation_times$system, FUN=function(x)mean(system<=x))

results<-list()
whichNeg<-rep(0, length(true_failure_rates))
nNeg<-rep(0,3)
names(whichNeg)<-names(true_failure_rates)
names(nNeg)<-c('10', '50', '500')

set.seed(50)

for (n in ns){
  nc=as.character(n)
  results[[nc]]<-list()
  for (c in censoring){
    cc<-as.character(c)
    results[[nc]][[cc]]<-list()
    if (c){
      startEndValve<-qweibull(c(.7,.99), 4,20)
      startEndGen<-qchisq(c(.7,.99), 3)
      startEndPlumbing<-quantile(plumbing,c(.7,.99))
      startEndPower<-quantile(power,c(.7,.99))
      startEndSystem<-quantile(system,c(.7,.99))
    }
    # cutoff<-qweibull(c, 4, 20)
    # g_cutoff<-qchisq(c, 3)
    # plumb_cutoff<-quantile(plumbing, c)
    # power_cutoff<-quantile(power, c)
    # sys_cutoff<-quantile(system, c)
    startOn=1
    for(i in startOn:10){
    results[[nc]][[cc]][[i]]<-list()
    dataList<-list()

    if (c){
      cutoff<-runif(n, startEndValve[1], startEndValve[2])
      g_cutoff<-runif(n, startEndGen[1], startEndGen[2])
      plumb_cutoff<-runif(n, startEndPlumbing[1], startEndPlumbing[2])
      power_cutoff<-runif(n, startEndPower[1], startEndPower[2])
      sys_cutoff<-runif(n, startEndSystem[1], startEndSystem[2])

    }else{
      cutoff<-g_cutoff<-plumb_cutoff<-power_cutoff<-sys_cutoff<-rep(Inf, n)
    }
    valve_data<-rweibull(n, 4, 20)
    censored<-valve_data>cutoff
    valve_data[censored]<-cutoff[censored]
    dataList$valve<-cbind(valve_data, as.numeric(!censored))

    gen_data<-rchisq(n, 3)
    g_censored<-gen_data>g_cutoff
    gen_data[g_censored]<-g_cutoff[g_censored]
    dataList$generator<-cbind(gen_data, as.numeric(!g_censored))

    plumb_data<-sample(plumbing, n)
    p_cens<-plumb_data>plumb_cutoff
    plumb_data[p_cens]<-plumb_cutoff[p_cens]
    dataList$plumbing<-cbind(plumb_data, as.numeric(!p_cens))

    power_data<-sample(power, n)
    pow_cens<-power_data>power_cutoff
    power_data[pow_cens]<-power_cutoff[pow_cens]
    dataList$power<-cbind(power_data, as.numeric(!pow_cens))

    sys_data<-sample(system, n)
    sys_cens<-sys_data>sys_cutoff
    sys_data[sys_cens]<-sys_cutoff[sys_cens]
    dataList$system<-cbind(sys_data, as.numeric(!sys_cens))

    posteriors=estimateSystemReliability(file="Simulation/System1.txt",
                                         priorList =priorList,
                                         dataList = dataList)

    results[[nc]][[cc]][[i]]$Bias<-sapply(names(evaluation_times), function(part){
      evaluate_centering_measures(posteriors[[part]], evaluation_times[[part]])-
                                          true_failure_rates[[part]]})

    results[[nc]][[cc]][[i]]$Coverage<-matrix(0, nrow=3, ncol=length(posteriors))
    results[[nc]][[cc]][[i]]$Coverage2<-matrix(0, nrow=3, ncol=length(posteriors))

    for (part in names(posteriors)){
      if (any(posteriors[[part]]$precision<0, na.rm =T)){
        whichNeg[names(whichNeg)==part]<-whichNeg[names(whichNeg)==part]+1
        nNeg[names(nNeg)==n]<-nNeg[names(nNeg)==n]+1
        #startOn=i+1
      #  fsdf
      }

      confint = bspConfint(posteriors[[part]], evaluation_times[[part]])
      results[[nc]][[cc]][[i]]$Coverage[,which(part==names(posteriors))] <-
        as.numeric(confint[1,]<= true_failure_rates[[part]] &
                     true_failure_rates[[part]]<= confint[2,])

      #Just using the beta to get confidence intervals.
      confint2 = bspConfint2(posteriors[[part]], evaluation_times[[part]])
      results[[nc]][[cc]][[i]]$Coverage2[,which(part==names(posteriors))] <-
        as.numeric(confint2[1,]<= true_failure_rates[[part]] &
                     true_failure_rates[[part]]<= confint2[2,])
      #if (all(results[[nc]][[cc]][[i]]$Coverage[,which(part==names(posteriors))]==0))sdkfljsd
    }
    #if (mean(results[[nc]][[cc]][[i]]$Coverage)<.5)askldfj


    }
    #save(results, whichNeg, nNeg, file="Simulation/Simulation1Results.RData")

  }
}

#save(results, whichNeg, nNeg, file="Simulation/Simulation1Results.RData")

extractor<-function(results, measurement, n, cens, part="all", time="all"){
  n<-as.character(n)
  cens<-as.character(cens)
  part_names<-c("valve", "generator", "plumbing", "power", "system")
  time_names<-c("first", "middle", "last")
  if (part=="all") columns=1:5 else columns = which(part==part_names)
  if (time=="all")rows=1:3 else rows= which(time==time_names)
  sapply(results[[n]][[cens]], FUN=function(x)x[[measurement]][rows, columns])
}

point<-mean(extractor(results, "Coverage", 10, .8, part="all", time="all"))
point+c(-1,1)*qnorm(.975)*sqrt(point*(1-point)/1000)

hist(extractor(results, "Bias", 50, T, "system", "last"))
hist(extractor(newresults, "Bias", 500, .8, "valve", "last"))

#par(ask=T)
for (n in ns){
  for (c in censoring){
    for (p in names(evaluation_times)){
      for (h in c("first", "middle", "last")){
        bias=extractor(results, "Bias", n, c, part=p, time=h )
        cove=extractor(results, "Coverage", n, c, part=p, time=h )

        #print(bias)
        if(p=="valve" &h=="last")hist(bias, main=paste(n, c))
        print(paste(n, c, p, h))
        #hist(bias, main=paste("n=",n, c, p, h))

        #print(mean(bias))
        print(mean(cove))
      }
    }

  }
}
extractor(results, "Coverage", 500, .8, part="system", time="all" )
a=results
mean(extractor(results, "Coverage", 10, 1, "valve", "last"))

survival::Surv()
#valve ~ Weibull(4, 10)
curve(dweibull(x, 4, 20),  xlim=c(0,30))
set.seed(100)


#Generator~chisq(3)
curve(dchisq(x, 3), xlim=c(0,30))







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

gg_ecdf<-function(vec){
  stat_ecdf(data=data.frame(x=vec), aes(x), geom="step", inherit.aes = F)
}

plot(posteriors$valve, withConfInt = T)+gg_curve(function(x)pweibull(x, 4,20), 12,25)

plot(posteriors$generator, withConfInt = T)+gg_curve(function(x)pchisq(x, 3), 0,8)


plot(posteriors$system, withConfInt = T)+gg_ecdf(true_system)

plot(posteriors$plumbing, withConfInt = T)+gg_ecdf(plumbing)
plot(posteriors$power, withConfInt = T)+gg_ecdf(power)
