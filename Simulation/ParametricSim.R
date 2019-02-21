#parmet Bayes Sim

para_cov<-rep(0, 1000)
nonpara_cov<-rep(0, 1000)
for (i in 1:1000){
  times<-rexp(30, 5)

  #Our parameter of interest in parametric form is theta
  #Set a prior on theta of Gamma(.01,.01)
  #Then our posterior on theta is gamma with parameters:
  alpha= .01+length(times)
  beta = .01+sum(times)

  #Approx posterior dist of theta
  thetas=rgamma(10000, alpha, beta)

  #Get a 95% CI on the probability of having failed by time=.4

  ##Now let's compare that to the BSP method
  prior<-bsp(c(.1,1,2), c(.05,.8, .99), precision = .01)
  posterior<-bspPosterior(prior, cbind(times,1))

  #Now let's compare the series

  ##First let's calculate the truth Series~exp(10)
  truth<-pexp(.1,10)
  series_thetas<-sample(thetas, 10000, replace=T)+sample(thetas,10000, replace=T)
  ci =quantile(pexp(.2, series_thetas),c(.025,.975))
  para_cov[i]<-as.numeric(ci[1]<truth & truth<ci[2])

  a=bspFromMoments(E1E2_series(posterior, posterior))
  ci<-bspConfint(a, .1)
  nonpara_cov[i]<-as.numeric(ci[1,]<truth & truth<ci[2,])

}

mean(para_cov)
mean(nonpara_cov)
