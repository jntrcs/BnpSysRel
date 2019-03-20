#SimplifiedSim.R

##Simulation questions:


require(BnpSysRel)
set.seed(10)

##Calculate the truth
valve<-rweibull(100000, 4,20)
generator<-rweibull(100000,4,20)
generator2<-rweibull(100000,4,20)

divSystem<-pmin(valve, generator)
repSystem<-pmin(generator2, generator)

ns<-c(10, 50, 500)
censoring<-c(1) #100% of the data and 75% of the data will be uncensored

evaluation_times<-list()
evaluation_times$valve<-c(13.7,  18.2,  22.5)
evaluation_times$generator<-c(13.7,  18.2,  22.5)
evaluation_times$divSystem<-quantile(divSystem, c(.2,.5,.8))
evaluation_times$repSystem<-quantile(repSystem, c(.2,.5,.8))


true_failure_rates<-list()
true_failure_rates$valve<-pweibull(evaluation_times$valve, 4, 20)
true_failure_rates$generator<-pweibull(evaluation_times$generator, 4,20)
true_failure_rates$divSystem<-true_failure_rates$repSystem<-c(.2,.5,.8)

allTimes=sort(unique(unlist(evaluation_times)))
valve_prior<-bsp(allTimes, pweibull(allTimes, 4, 20), .2)
generator_prior<-bsp(allTimes, pweibull(allTimes, 4,20), .2)
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

    for(i in 1:200){
      results[[nc]][[cc]][[i]]<-list()
      dataList<-list()
      valve_data<-rweibull(n, 4, 20)
      censored<-valve_data>cutoff
      valve_data[censored]<-cutoff
      dataList$valve<-cbind(valve_data, as.numeric(!censored))

      gen_data<-rweibull(n, 4,20)
      g_censored<-gen_data>cutoff
      gen_data[g_censored]<-cutoff
      dataList$generator<-cbind(gen_data, as.numeric(!g_censored))

      posteriors=estimateSystemReliability(file="Simulation/SimpleRepeat.txt",
                                           priorList =priorList,
                                           dataList = dataList)

      results[[nc]][[cc]][[i]]$Bias<-sapply(names(evaluation_times), function(part){
        evaluate_centering_measure(posteriors[[part]], evaluation_times[[part]])-
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
    save(results, file="Simulation/SimpleRepeat.RData")

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

hist(extractor(results, "Bias", n=10, cens=1, part="divSystem", time="first"))
hist(extractor(results, "Bias", n=10, cens=1, part="repSystem", time="first"))
hist(extractor(results, "Bias", n=10, cens=1, part="divSystem", time="middle"))
hist(extractor(results, "Bias", n=10, cens=1, part="repSystem", time="middle"))
hist(extractor(results, "Bias", n=10, cens=1, part="divSystem", time="last"))
hist(extractor(results, "Bias", n=10, cens=1, part="repSystem", time="last"))
hist(extractor(results, "Bias", n=500, cens=1, part="divSystem", time="first"))
hist(extractor(results, "Bias", n=500, cens=1, part="repSystem", time="first"))
hist(extractor(results, "Bias", n=500, cens=1, part="divSystem", time="middle"))
hist(extractor(results, "Bias", n=500, cens=1, part="repSystem", time="middle"))
hist(extractor(results, "Bias", n=500, cens=1, part="divSystem", time="last"))
hist(extractor(results, "Bias", n=500, cens=1, part="repSystem", time="last"))


mean(extractor(results, "Coverage", n=10, cens=1, part="divSystem", time="first"))
mean(extractor(results, "Coverage", n=10, cens=1, part="repSystem", time="first"))
mean(extractor(results, "Coverage", n=10, cens=1, part="divSystem", time="middle"))
mean(extractor(results, "Coverage", n=10, cens=1, part="repSystem", time="middle"))
mean(extractor(results, "Coverage", n=10, cens=1, part="divSystem", time="last"))
mean(extractor(results, "Coverage", n=10, cens=1, part="repSystem", time="last"))
mean(extractor(results, "Coverage", n=500, cens=1, part="divSystem", time="first"))
mean(extractor(results, "Coverage", n=500, cens=1, part="repSystem", time="first"))
mean(extractor(results, "Coverage", n=500, cens=1, part="divSystem", time="middle"))
mean(extractor(results, "Coverage", n=500, cens=1, part="repSystem", time="middle"))
mean(extractor(results, "Coverage", n=500, cens=1, part="divSystem", time="last"))
mean(extractor(results, "Coverage", n=500, cens=1, part="repSystem", time="last"))
