##Beta-Binomial Simulation

alphas<-seq(0,4, length.out=41)
ns<-7:18
true.ps<-1:19/20
coverageResults<-array(0,dim=c(length(alphas), length(ns), length(true.ps)),
                       dimnames=list(Alphas=alphas, Ns=ns, ProbSuccess=true.ps))

for (alpha in alphas){
  for (n in ns){
    for (p in true.ps){
      dataSets=rbinom(10000, n, p)
      lowers=qbeta(c(.025), alpha+dataSets, alpha+n-dataSets)
      uppers=qbeta(c(.975), alpha+dataSets, alpha+n-dataSets)
      coverageResults[as.character(alpha), as.character(n), as.character(p)]<-
        mean(lowers<p & p<uppers)
    }
  }
}
save(coverageResults, file="BetaBinomResults.RData")

require(reshape2)
meltedCov<-melt(coverageResults)
save(meltedCov, file="Simulation/ShinyData.RData")

meltedCov$NProb<-paste(meltedCo)
ggplot(meltedCov, aes(y=value, x=Alphas, color=ProbSuccess, shape=factor(Ns)))+geom_point()
