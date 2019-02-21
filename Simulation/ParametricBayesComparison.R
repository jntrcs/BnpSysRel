###Analogy of using a parametric Bayesian model with our method

#Suppose we have two iid components, both following an Exp(5)
#We know from theory these components in series will follow a Exp(10)

#Now suppose we are able to collect some data
#We will analyze the data and make inference using a parameteric model and our BSP set up
times<-rexp(20, 5)

#Our parameter of interest in parametric form is theta
#Set a prior on theta of Gamma(.01,.01)
#Then our posterior on theta is gamma with parameters:
alpha= .01+length(times)
beta = .01+sum(times)

#Approx posterior dist of theta
thetas=rgamma(10000, alpha, beta)


#Get a 95% CI on the probability of having failed by time=.4
quantile(pexp(.4, thetas),c(.025,.975))

##Now let's compare that to the BSP method
prior<-bsp(c(.1,1,2), c(.05,.8, .99), precision = .01)
posterior<-bspPosterior(prior, cbind(times,1))
bspConfint(posterior, .4)

#Now let's compare the series

##First let's calculate the truth Series~exp(10)
truth<-pexp(.1,10)
series_thetas<-sample(thetas, 10000, replace=T)+sample(thetas,10000, replace=T)
quantile(pexp(.2, series_thetas),c(.025,.975))

a=bspFromMoments(E1E2_series(posterior, posterior))
bspConfint(a, .1)

dataset1<-rexp(10000, thetas)
dataset2<-rexp(10000, thetas)
post_pred_series<-pmin(dataset1, dataset2)
mean(post_pred_series<=.1)
