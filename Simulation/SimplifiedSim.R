#SimplifiedSim.R

##Simulation questions:


require(BnpSysRel)
set.seed(10)

##Calculate the truth
valve<-rweibull(100000, 4,20)
generator<-rchisq(100000, 18)

seriesSystem<-pmin(valve, generator)
parSystem<-pmax(valve, generator)

ns<-c(10, 50, 500)
censoring<-c(.75, 1) #100% of the data and 75% of the data will be uncensored

evaluation_times<-list()
evaluation_times$valve<-c(13.7,  18.2,  22.5)
evaluation_times$generator<-c(12.8, 17.3, 22.8)
evaluation_times$SeriesSystem<-quantile(seriesSystem, c(.2,.5,.8))
evaluation_times$ParSystem<-quantile(parSystem, c(.2,.5,.8))


true_failure_rates<-list()
true_failure_rates$valve<-pweibull(evaluation_times$valve, 4, 20)
true_failure_rates$generator<-pchisq(evaluation_times$generator, 18)
true_failure_rates$SeriesSystem<-true_failure_rates$ParSystem<-c(.2,.5,.8)

allTimes=sort(unique(unlist(evaluation_times)))
valve_prior<-bsp(allTimes, pweibull(allTimes, 4, 20), .2)
generator_prior<-bsp(allTimes, pchisq(allTimes, 18), .2)
priorList=list(valve=valve_prior, generator= generator_prior)
#plot(generator_prior, withConfInt = T, withPaths=T)+gg_ecdf(generators[,1])




results<-list()

for (n in ns){
  nc=as.character(n)
  results[[nc]]<-list()
  for (c in censoring){
    cc<-as.character(c)
    results[[nc]][[cc]]<-list()
    cutoff<-qweibull(c, 4, 20)
    g_cutoff<-qchisq(c, 18)

    for(i in 1:200){
      results[[nc]][[cc]][[i]]<-list()
      dataList<-list()
      valve_data<-rweibull(n, 4, 20)
      censored<-valve_data>cutoff
      valve_data[censored]<-cutoff
      dataList$valve<-cbind(valve_data, as.numeric(!censored))

      gen_data<-rchisq(n, 18)
      g_censored<-gen_data>g_cutoff
      gen_data[g_censored]<-g_cutoff
      dataList$generator<-cbind(gen_data, as.numeric(!g_censored))

      posteriors=estimateSystemReliability(file="Simulation/SimplifiedSystem.txt",
                                           priorList =priorList,
                                           dataList = dataList)

      results[[nc]][[cc]][[i]]$Bias<-sapply(names(evaluation_times), function(part){
        evaluate_centering_measures(posteriors[[part]], evaluation_times[[part]])-
          true_failure_rates[[part]]})



      results[[nc]][[cc]][[i]]$Coverage<-matrix(0, nrow=3, ncol=length(posteriors))
      for (part in names(posteriors)){
        if (any(posteriors[[part]]$precision<0))askdfj
        confint = bspConfint(posteriors[[part]], evaluation_times[[part]])
        results[[nc]][[cc]][[i]]$Coverage[,which(part==names(posteriors))] <-
          as.numeric(confint[1,]<= true_failure_rates[[part]] &
                       true_failure_rates[[part]]<= confint[2,])
        #if (all(results[[nc]][[cc]][[i]]$Coverage[,which(part==names(posteriors))]==0))sdkfljsd
      }

    }
    save(results, file="SimplifiedSimResults.RData")

  }
}

#save(results, file="Simulation/SimpleSim.RData")

extractor<-function(results, measurement, n, cens, part="all", time="all"){
  n<-as.character(n)
  cens<-as.character(cens)
  part_names<-colnames(results[[n]][[cens]][[1]][["Bias"]])
  time_names<-c("first", "middle", "last")
  if (part=="all") columns=1:length(part_names) else columns = which(part==part_names)
  if (time=="all")rows=1:length(time_names) else rows= which(time==time_names)
  sapply(results[[n]][[cens]], FUN=function(x)x[[measurement]][rows, columns])
}

for (c in censoring){
for (n in ns){
  for (part in c("valve", "generator", "ParSystem", "SeriesSystem")){
    for(area in c("first", "middle", "last")){
      print(paste("Percentage of uncensored data:", c))
      print(paste("N=",n))
      print(part)
      print(area)
      print(paste("Coverage:", mean(extractor(results, "Coverage", n, 1, part=part, time=area))) )
    }
  }
}
}
