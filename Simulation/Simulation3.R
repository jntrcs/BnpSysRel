#Sim 3
##Simulation questions:

#1. Are the estimates consistent for many components in series/parallel

require(BnpSysRel)

##Calculate the truth
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

allTimes<-sort(unique(unlist(evaluation_times)))

valve_prior<-bsp(evaluation_times$valve, pweibull(evaluation_times$valve, 4,20), .2)
generator_prior<-bsp(evaluation_times$generator, pchisq(evaluation_times$generator,3), .2)

priorList=list(valve1=valve_prior,valve2=valve_prior,valve3=valve_prior,
               valve4=valve_prior,valve5=valve_prior,valve6=valve_prior,
               valve7=valve_prior,valve8=valve_prior,valve9=valve_prior,
               valve10=valve_prior,
               generator1= generator_prior, generator2= generator_prior,
               generator3= generator_prior, generator4= generator_prior,
               generator5= generator_prior, generator6= generator_prior,
               generator7= generator_prior, generator8= generator_prior,
               generator9= generator_prior, generator10= generator_prior)
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
names(nNeg)<-c('10', '50', '250')

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
    for(i in startOn:1000){
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

      for (nvalve in 1:10){
        valve_data<-rweibull(n, 4, 20)
        censored<-valve_data>cutoff
        valve_data[censored]<-cutoff[censored]
        dataList[[paste0("valve", nvalve)]]<-cbind(valve_data, as.numeric(!censored))
      }

      for (ngen in 1:10){
        gen_data<-rchisq(n, 3)
        g_censored<-gen_data>g_cutoff
        gen_data[g_censored]<-g_cutoff[g_censored]
        dataList[[paste0("generator", ngen)]]<-cbind(gen_data, as.numeric(!g_censored))
      }
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

      posteriors=estimateSystemReliability(file="Simulation/System3.txt",
                                           priorList =priorList,
                                           dataList = dataList)

      results[[nc]][[cc]][[i]]$Bias<-sapply(names(posteriors), function(part){
        part=gsub("^|\\d+$", "", part)
        evaluate_centering_measures(posteriors[[part]], evaluation_times[[part]])-
          true_failure_rates[[part]]})

      results[[nc]][[cc]][[i]]$Coverage<-matrix(0, nrow=3, ncol=length(posteriors))
      results[[nc]][[cc]][[i]]$Coverage2<-matrix(0, nrow=3, ncol=length(posteriors))

      for (part in names(posteriors)){
        #   if (any(posteriors[[part]]$precision<0, na.rm =T)){
        #     whichNeg[names(whichNeg)==part]<-whichNeg[names(whichNeg)==part]+1
        #     nNeg[names(nNeg)==n]<-nNeg[names(nNeg)==n]+1
        #     #startOn=i+1
        #   #  fsdf
        #   }
        genericPart=gsub("^|\\d+$", "", part)
        confint = bspConfint(posteriors[[part]], evaluation_times[[genericPart]])
        results[[nc]][[cc]][[i]]$Coverage[,which(part==names(posteriors))] <-
          as.numeric(confint[1,]<= true_failure_rates[[genericPart]] &
                       true_failure_rates[[genericPart]]<= confint[2,])

        #Just using the beta to get confidence intervals.
        confint2 = bspConfint2(posteriors[[part]], evaluation_times[[genericPart]])
        results[[nc]][[cc]][[i]]$Coverage2[,which(part==names(posteriors))] <-
          as.numeric(confint2[1,]<= true_failure_rates[[genericPart]] &
                       true_failure_rates[[genericPart]]<= confint2[2,])
        #if (all(results[[nc]][[cc]][[i]]$Coverage[,which(part==names(posteriors))]==0))sdkfljsd
      }
      #if (mean(results[[nc]][[cc]][[i]]$Coverage)<.5)askldfj


    }
    save(results, file="Simulation/Simulation3.RData")

  }
}

#save(results, whichNeg, nNeg, file="Simulation/Simulation1Results.RData")

temp<-posteriors$generator1
for (i in 2:10){
  a=E1E2_parallel(temp, posteriors[[paste0("generator",i)]])
  print(a)
  temp=bspFromMoments(a)
  plot(temp)
}

for (i in 1:10){
  plot(posteriors[[paste0("generator",i)]])
}
