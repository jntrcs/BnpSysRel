

alpha=.85
prior=bsp(1:200/100, pexp(1:200/100), precision = alpha)
timePeriod=qexp(.5)
truth<-pexp(timePeriod)

cov<-numeric(0)
alphas<-seq(.5, 1, length.out = 100)
for (prec in alphas){
alpha=prec
inInterval<-numeric(5000)
counts<-numeric(5000)
for (i in 1:5000){
data=rexp(10)
#posterior<-bspPosterior(prior, cbind(data, 1))
#a=bspConfint(posterior, timePeriod)
F0<-evaluate_centering_measure(prior, timePeriod)
counts[i]=sum(data<=timePeriod)
a=qbeta(c(.025,.975), alpha*evaluate_centering_measure(prior, timePeriod) +
        sum(data<=timePeriod), alpha*(1-evaluate_centering_measure(prior, timePeriod))+
                            sum(data>timePeriod))
inInterval[i]<-a[1]<=truth & a[2]>=truth
}
cov<-c(cov, mean(inInterval))
}
plot(cov~alphas)

alphas<-seq(2, 0, length.out = 40)
coverage<-numeric(0)
for (alpha in alphas){
inInterval<-numeric(10000)
allSame=logical(10000)
for (i in 1:10000){
  true.p<-.5#runif(1)
  data=rbinom(1,12, true.p)
  allSame[i]= data==12 | data==0

  a=qbeta(c(.025, .975), data+alpha, 12-data+alpha)
  inInterval[i]<-a[1]<=true.p & a[2]>=true.p
  #if(allSame[i]&inInterval[i])ashj

}
coverage<-c(coverage, mean(inInterval))
}
plot(coverage~alphas)


