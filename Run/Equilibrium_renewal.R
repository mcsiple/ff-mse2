# Equilibrium Renewal Model

lh.pars = lh.test
# Function to calculate spawning biomass per recruit
sbprFN <- function(lh.pars, F=0) {
  with (lh.pars,{
    n.age<-rep(NA, length(ages))
    n.age[1]<-1
    for (i in 2:length(ages)){
      n.age[i]<- n.age[i-1]*exp(-F*selectivity[i-1,2]-M)
    }
    sbpr <- sum(n.age * maturity*w.at.age)
  })
}

yprFN<- function(lh.pars, F=0) {
  with (lh.pars,{
    n.age<-rep(NA, length(ages))
    y.age <- rep(NA, length(ages))
    n.age[1]<-1
    f.at.age <- F * selectivity[,2]
    for (i in 2:length(ages)){
      n.age[i]<- n.age[i-1]*exp(-f.at.age[i-1]-M)
    }
    b.age <- n.age * w.at.age
    # apply Baranov catch equation
    for (i in 1:length(ages)){
      y.age[i] <- b.age[i]*f.at.age[i] / (f.at.age[i] + M) * (1-exp(-f.at.age[i] - M))
    }
    maxage <- length(ages)
    ypr <- sum(y.age)
    
  })
}

rzero <- 1E10 # arbitrary
h <- 0.9

sbpr.nofishing <- sbprFN(lh.pars, F=0)
bhb <- (h-0.2)/((1-h)*(0.2*sbpr.nofishing*rzero))
bha <- (1+bhb*sbpr.nofishing*rzero)/(sbpr.nofishing)

fs <- seq(0,5,by=0.001)



ys <- rep(NA,length(fs))
spb <- rep(NA,length(fs))
for (i in 1:length(fs)){
  f <- fs[i]
  sbpr <- sbprFN(lh.pars, F=f)
  ypr <- yprFN(lh.pars, F=f)
  r.tmp <-max(c(0 ,(bha*sbpr-1)/(bhb*sbpr)))
  ys[i]<-r.tmp * ypr
  spb[i] <- sbpr * r.tmp
  }
  
plot(fs,ys,type="l")
(Fmsy <- fs[which(ys == max(ys))] )
