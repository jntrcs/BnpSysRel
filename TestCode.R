##Tests
a=bsp(c(.5,5), centeringMeasure = c(.1,.9), precision = .001)
plot(a)
a
print(a)
a$evaluate_alpha(0:6)
a$evaluate_g(0:6)

data<-matrix(c(1, 1,
               2.2, 1,
               3, 1
              ), byrow=T, ncol=2)

b=bspPosterior(a, data)
plot(b)
plot(b$precision~b$support, ylim=c(2.8, 3.2), type='l')
b

datac<-matrix(c(1, 1,
               2, 0,
               3, 1
), byrow=T, ncol=2)

ce=bspPosterior(a, datac)
plot(ce)
plot(ce$precision~ce$support, type='l')
ce


###EXAMPLE 3 FROM OVERLEAF
prior = bsp(support=c(0,2), c(0, .8), precision=1)
data<-matrix(c(1, 1,
               2, 1,
               2, 0,
               3, 1
), byrow=T, ncol=2)
posterior<-bspPosterior(prior, data)
posterior$evaluate_alpha(seq(0, 3.5, by=.5))
posterior$evaluate_g(seq(0, 3.5, by=.5))
##Matches what is in the paper


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




