require(BnpSysRel)

inSeries<-function(draws1, draws2){
  1-(1-draws1)*(1-draws2)
}

inParallel<-function(draws1, draws2){
  draws1*draws2
}
makeComparison<-function(posterior1, posterior2, time, series=T, location="bottomleft"){
  set.seed(400)
  if(series){
    moments<-E1E2_series
    relation<-inSeries
  }else{
    moments<-E1E2_parallel
    relation<-inParallel
  }
  approx<-bspFromMoments(moments(posterior1, posterior2))
  if(any(approx$precision<0))warning("Negative precision present in approximation")
  approx_samples<-bspSampling(approx)
  approx_dist<-samplesAt(approx_samples, time)

  valveAt<-samplesAt(bspSampling(posterior1), time)
  genAt<-samplesAt(bspSampling(posterior2), time)
  exact_dist<-relation(valveAt, genAt)
  slice=data.frame(Approx=approx_dist, Actual=exact_dist)
  ggplot(slice)+geom_density(aes(x=Approx, col="Approx."))+
          geom_density(aes(x=Actual,col="Exact"), linetype=2 )+
          #ggtitle(paste("Marginal at",round(time,2), "in", ifelse(series, "Series", "Parallel")))+
          xlab("F(X)")+theme(text = element_text(size=19))
  # plot(density(approx_dist), main=paste0("Marginal at ", round(time,3), " in ",
  #                                        ifelse(series, "Series", "Parallel")),
  #      xlab="F(T)",cex=15)
  # lines(density(exact_dist), col="red")
  # legend(location, lty=c(1,1), legend=c("Approximation", "Exact"), col=c("Black","red"))
}
##Using the same set up as before:
valves<-replicate(3,rweibull(100000, 4,20))
generators<-replicate(2, rchisq(100000, 3))
plumbing<-apply(valves, 1, min)
power<-apply(generators, 1, max)
system<-pmin(plumbing, power)

evaluation_times<-list()
evaluation_times$valve<-qweibull(c(.2,.5,.8), 4, 20)
evaluation_times$generator<-qchisq(c(.2,.5,.8),3)
evaluation_times$plumbing<-quantile(plumbing, c(.2,.5,.8))
evaluation_times$power<-quantile(power, c(.2,.5,.8))
evaluation_times$system<-quantile(system, c(.2,.5,.8))

allTimes<-sort(unique(unlist(evaluation_times)))



#1
valve_prior<-bsp(evaluation_times$valve, pweibull(evaluation_times$valve, 4,20), .2)
generator_prior<-bsp(evaluation_times$generator, pchisq(evaluation_times$generator,3), .2)
valve_data<-rweibull(10, 4,20)
valve_data2<-rweibull(10, 4,20)
valve_posterior<-bspPosterior(valve_prior, cbind(valve_data,1))
generator_posterior<-bspPosterior(valve_prior, cbind(valve_data2,1))
makeComparison(valve_posterior, generator_posterior, evaluation_times$valve[1], series=T,
               "topright") +ggsave("Simulation/SeriesValves.pdf")
makeComparison(valve_posterior, generator_posterior, evaluation_times$valve[2], series=F,
               "topright")+ggsave("Simulation/ParallelValves.pdf")

#2
valve_data<-rnorm(10,29, 1)
gen_data<-rnorm(100, 30, 1)
valve_prior<-bsp(29:32, seq(.5,.99,length.out = 4), .2)
valve_posterior<-bspPosterior(valve_prior, cbind(valve_data,0))
generator_prior<-bsp(29:32, seq(.2, .96, length.out=4), .2)
generator_posterior<-bspPosterior(generator_prior, cbind(gen_data,1))
makeComparison(valve_posterior, generator_posterior, 31.5, series=T, "topleft")+
  ggsave("Simulation/BimodalBSPs.pdf")


#3
valve_prior<-bsp(evaluation_times$valve, pweibull(evaluation_times$valve, 4,20), .2)
generator_prior<-bsp(evaluation_times$generator, pchisq(evaluation_times$generator,3), .2)
valve_data<-rweibull(10, 4,20)
valve_data2<-rweibull(1000, 4,20)
valve_posterior<-bspPosterior(valve_prior, cbind(valve_data,1))
generator_posterior<-bspPosterior(valve_prior, cbind(valve_data2,1))
makeComparison(valve_posterior, generator_posterior, 20, series=T,
               "topright") +ggsave("Simulation/DisparatePrecisions.pdf")
makeComparison(valve_posterior, generator_posterior, 20, series=F,
               "topright")#+ggsave("Simulation/ParallelValves.pdf")

#4
valve_prior<-bsp(evaluation_times$valve, pweibull(evaluation_times$valve, 4,20), .2)
generator_prior<-bsp(evaluation_times$generator, pchisq(evaluation_times$generator,3), .2)
valve_data<-rweibull(100, 4,20)
valve_data2<-rweibull(100, 8,40)
valve_posterior<-bspPosterior(valve_prior, cbind(valve_data,1))
generator_posterior<-bspPosterior(valve_prior, cbind(valve_data2,1))
makeComparison(valve_posterior, generator_posterior, 25, series=T,
               "topright") #+ggsave("Simulation/SeriesValves.pdf")
makeComparison(valve_posterior, generator_posterior, 25, series=F,
               "topright")+ggsave("Simulation/OppositeTails.pdf")
#3
require(mvtnorm)
set.seed(3)
rho=0
draws<-rmvnorm(100, mean=c(30,29), sigma=matrix(c(4, rho,rho,4), byrow = T, nrow=2))
comp1.obs<-draws[1:10,1]
comp2.obs<-draws[,2]
censor.year.1<-sample(28:35, replace=T, length(comp1.obs))
censor.year.2<-sample(28:35, replace=T, length(comp2.obs))
comp1.obs<-pmin(comp1.obs, censor.year.1)
comp2.obs<-pmin(comp2.obs, censor.year.2)
censored1<-as.numeric(comp1.obs!=censor.year.1)
censored2<-as.numeric(comp2.obs!=censor.year.2)
prior1<-bsp(20,.005, .1)
prior2<-bsp(20,.005,.1)
comp1<-bspPosterior(prior1, cbind(comp1.obs, censored1))
comp2<-bspPosterior(prior2, cbind(comp2.obs, censored2))
approx=bspFromMoments(E1E2_series(comp1,comp2))
makeComparison(comp1, comp2, time=31.5, series=T, "topleft")



