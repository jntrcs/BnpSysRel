##Tests
require(BnpSysRel)
a=bsp(c(1:3), centeringMeasure = c(.1,.9, .98), precision = 2, calculateMoments = T)
evaluate_second_moment(a, c(0, .5,1,1.5, 3.5))
E1E2_series(a, a)
b=bspSampling(a, 100)

plot(a)
a
print(a)

evaluate_centering_measures(a, seq(.5, 3.5, by=.5))

data<-matrix(c(1, 1,
               2.2, 1,
               3, 1
              ), byrow=T, ncol=2)
a=createUninformativePrior(data)

b=bspPosterior(a, data)
plot(b)
print(b)

datac<-matrix(c(1, 1,
               2, 0,
               3, 1
), byrow=T, ncol=2)

ce=bspPosterior(a, datac)
ce
evaluate_precision(ce, c(1.9,2,2.1))



#Example 1 from paper
prior = bsp(support=c(0,2), c(0, .8), precision=.00000000)
data<-matrix(c(1, 1,
               2,1,
               3, 1
), byrow=T, ncol=2)
posterior<-bspPosterior(prior, data)
seq(0, 4.5, by=.5)
evaluate_precision(posterior, seq(0, 4.5, by=.5))
evaluate_centering_measures(posterior,seq(0, 4.5, by=.5))
E1E2(posterior)

#Example 2 from paper
prior = bsp(support=c(2,4), c(.8,.991), precision=.00000000)
data<-matrix(c(1, 1,
               2,1,
               2,0,
               3, 1
), byrow=T, ncol=2)
posterior<-bspPosterior(prior, data)
seq(0, 4.5, by=.5)
evaluate_precision(posterior, seq(0, 4.5, by=.5))
evaluate_centering_measures(posterior,seq(0, 4.5, by=.5))
E1E2(posterior)

###EXAMPLE 3 FROM OVERLEAF
prior = bsp(support=c(2), c(.8), precision=1)
data<-matrix(c(1, 1,
               2, 1,
               2, 0,
               3, 1
), byrow=T, ncol=2)
posterior<-bspPosterior(prior, data)
seq(0, 4.5, by=.5)
evaluate_precision(posterior, seq(0, 4.5, by=.5))
evaluate_centering_measures(posterior,seq(0, 4.5, by=.5))
E1E2(posterior)


##Matches what is in the paper

# #Example 4 from paper
# prior = bsp(support=c(0,2), c(0,.5), precision=c(1,2))
# data<-matrix(c(1, 1,
#                2,1,
#                2,0,
#                3, 1
# ), byrow=T, ncol=2)
# posterior<-bspPosterior2(prior, data)
# seq(1.5, 2.5, by=.5)
# evaluate_precision(posterior, seq(1.5, 2.5, by=.5))
# evaluate_centering_measures(posterior,seq(1.5, 2.5, by=.5))
# E1E2(posterior)
# #Experiment: look at the variance of F(T=3) for the cases where there is a censored
# #observation at 2, a fully observed obs at two, and nothing at two

##Example 5 from paper
prior = bsp(c(4), c(.9), precision=0.0001)
data<-matrix(c(1,1,
               2,0,
               3,0), byrow = T, ncol=2)
posterior<-bspPosterior(prior,data)
posterior
bspFromMoments(E1E2_series(posterior, posterior))

#Nothing
prior = bsp(support=c(0,3.5), c(0, ), precision=.1)
data<-matrix(c(1, 1,
               3, 1,
               3.5,1
), byrow=T, ncol=2)
posterior1<-bspPosterior(prior, data)
post_draws1<-bsp_sampling(posterior1)
var(post_draws1[4,])

#censored
prior = bsp(support=c(0,3.5), c(0, ), precision=.1)
data<-matrix(c(1, 1,
               2,0,
               2,1,
               3, 1,


), byrow=T, ncol=2)
posterior2<-bspPosterior(prior, data)
post_draws2<-bsp_sampling(posterior2)
var(post_draws2[4,])

#Fully observed
prior = bsp(support=c(0,2,4), c(0, .8,1), precision=1)
data<-matrix(c(1, 1,
               2,1,
               3, 1,
               3.5,1
), byrow=T, ncol=2)
posterior3<-bspPosterior(prior, data)
post_draws3<-bsp_sampling(posterior3)
var(post_draws3[4,])

a=(bspPosterior(posterior, datac))
for (i in 1:5){
  a= bspPosterior(a, datac)
}
a$evaluate_alpha(5)

#hwat does the Kaplan Meier say about this data?
require(survival)
dat<-Surv(c(1,2,2, 3),c(1,1,0,1))
KM0 <- survfit(dat ~ 1,  type="kaplan-meier", conf.int=F)
plot(KM0)
?plot.survfit

#prior ts are the location of the jumps in the centering measure
#prior gs are the centering measure at each g
make_calculator<-function(prior_alpha, prior_ts, prior_gs, data, censor){
  if (any(data!=sort(data)))stop("Data must be sorted (increasing)")
  single_t<-function(t){
    m<-function(data, t) sapply(t, FUN=function(x) sum(data>=x))
    j<-function(data, censor, t) sapply(t, FUN=function(x) sum(censor*(data==x)))
    previous<-unique(data[data<=t & censor])
    previous<-sort(unique(c(previous, prior_ts[prior_ts<=t])))
    prior_g<-function(t)sapply(t, FUN=function(t)c(0,prior_gs)[sum(prior_ts<=t)+1])
    if(length(previous)==0) {
      previous=0
    }
    g = 1-prod(1-(prior_alpha(previous)*(prior_g(previous)-prior_g(previous-.00001)) +j(data, censor, previous))/
                 (prior_alpha(previous)*(1-prior_g(previous-.00001))+m(data, previous)))
    alpha = (prior_alpha(t)*(1-prior_g(t))+m(data, t)-j(data, censor, t))/(1-g)

    return(list(CenteringMeasure=g, Alpha=alpha))
  }
  function(ts){
    a=sapply(ts, FUN=single_t)
    return(list(centeringMeasures=as.numeric(a[1,]), Precisions=as.numeric(a[2,]), support=ts))
  }
}

simplest_example<-make_calculator(function(x) rep(0, length(x)),
                                  prior_ts=c(1,2,3),
                                  prior_gs = c(.1,.9,.99), data=c(1,2,3),
                                  censor = c(1,1,1))
simplest_example(c(.5,2))
simplest_example(1)
simplest_example(2)
simplest_example(2.5)
simplest_example(3)


with_prior_example<-make_calculator(prior_alpha=function(x)ifelse(x<1.5, 1,1),# rep(1, length(x)),
                                  prior_ts=c(2),
                                  prior_gs = c(.9),
                                  data=c(1,2,2,3),
                                  censor = c(1,1,0,1))
with_prior_example(.2)

with_prior_example(.7)
with_prior_example(2.9)
with_prior_example(3)
with_prior_example(3.4)

with_prior_example(3.6)


##
data=sort(rlnorm(500, 5,.2))
prior<-bsp(c(0, 500), centeringMeasure=c(0, .99), precision=.1)
posterior=bspPosterior(prior, cbind(data,1))
draws=bspSamples(posterior, reps=1000)
matplot(x=posterior$support, draws, type='l')
curve(plnorm(x,5,.2), lwd=5,add=T, col="yellow")
m<-apply(draws, 1, mean)
lines(m~posterior$support, lwd=5, col="limegreen")
