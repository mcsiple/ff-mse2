# Equilibrium Renewal Model
# Sardine Parameters
ages.test <- 0:15
nages.test <- length(ages.test)
selectivity.test <- cbind(age=ages.test,selectivity = c(0.263,1.00,1.000,0.669,0.471,0.390,0.358,0.345,0.339,0.335,0.333,0.332,0.332,0.331,0.331,0.331))  #fishery selectivity from sardine assessment; see Table App.A.1. in Hurtado-Ferro & Punt 2014 
#Appendix A.1 caption: "Fleet-averaged selectivity (computed using the output of model X6e of Hill et al. [2012]). Results are shown for 2011, 2007-2011, and 2002-2011."
lh.pars <- list(M = 0.4,   #from sardine assessment (Hurtado-Ferro & Punt 2014)
                selectivity = selectivity.test,
                ages = ages.test,
                l.at.age = seq(4,28,length.out=nages.test), # Lengths are in mm...? Doesn't actually matter because lengths aren't used in the most basic model
                w.at.age = seq(0.014,0.180,length.out=nages.test)*0.001, #weights at age are in kg; multiplying by 0.001 gives mt
                maturity = c(0,0.4,0.85,0.99,rep(1,times=12)), # Roughly based on maturity ogives from Silva et al. 2006 ICES JMS.
                R0=5e9) 

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

fs <- seq(0,1.5,by=0.001)

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
Fmsy <- fs[which(ys == max(ys))]

sbpr <- sbprFN(lh.pars, F=Fmsy)
r.tmp <-max(c(0 ,(bha*sbpr-1)/(bhb*sbpr)))

# clunky code to simulate
n.years <- 100
F <- Fmsy

selectivity <- lh.pars$selectivity
ages <- lh.pars$ages
l.at.age <- lh.pars$l.at.age
w.at.age <- lh.pars$w.at.age
maturity <- lh.pars$maturity

M <- lh.pars$M


get.n.at.age <- function(lh.pars,R,n.tminus.1,F){
  with(lh.pars, {
  output <- c(R,n.tminus.1[-length(n.tminus.1)]*exp(-M - F*selectivity[-length(n.tminus.1),2]))
  })
}
getyield <- function(bvec,fvec, M) sum(bvec * fvec / (fvec + rep(M ,length(bvec))) * (rep(1,length(bvec)) - exp(-fvec - rep(M,length(bvec)))))

# Set Recruitment Devs, 3 different levels that change regime-like
Rlevels <- c(0.25,1,4)
Rdevs <- rep(NA,n.years)
Rdevs[1] <- 1
p.shift <- 0.1
shift.mag <- log(0.25)
probs <- runif(n.years,0,1)
for (i in 2:n.years) {
  if (probs[i] > p.shift) {
    Rdevs[i] <- Rdevs[i - 1]
  } else {
  Rdevs[i] <- sample(Rlevels[which(!Rlevels==Rdevs[i-1])],size=1)
}
}

# Set up matrix to store output
outmat <- matrix(NA,nrow= n.years, ncol = 16)
# get initial conditions (at F = Fmsy)

outmat[1,1]<-r.tmp
#initial conditions
for (i in 2:length(ages)){
  outmat[1,i]<-outmat[1,i-1]*exp(-F*selectivity[i-1,2]-M)
}
# get yield in first year
y<- rep(NA,n.years)
y[1]<-getyield(bvec=outmat[1,]*w.at.age, F*selectivity[,2],M)

for (i in 2:n.years) {
  sb <- sum(outmat[i-1,]*maturity * w.at.age)
  r.year <- (bha*sb / (1+bhb * sb) )*Rdevs[i-1]
  outmat[i,]<-get.n.at.age(lh.pars,R=r.year,n.tminus.1 = outmat[i-1,],F=Fmsy)
  y[i]<-getyield(bvec = outmat[i,]*w.at.age,F*selectivity[,2],M)
}

plot(1:n.years,y)
y_mean_msy <- mean(y)
print(mean(y))

# repeat at 0.5 * Fmsy
F=0.5 * Fmsy
y<- rep(NA,n.years)
y[1]<-getyield(bvec=outmat[1,]*w.at.age, F*selectivity[,2],M)

for (i in 2:n.years) {
  sb <- sum(outmat[i-1,]*maturity * w.at.age)
  r.year <- (bha*sb / (1+bhb * sb))*Rdevs[i-1]
  outmat[i,]<-get.n.at.age(lh.pars,R=r.year,n.tminus.1 = outmat[i-1,],F=F)
  y[i]<-getyield(bvec = outmat[i,]*w.at.age,F*selectivity[,2],M)
}

plot(1:n.years,y)
print(mean(y))
y_mean_0.5msy <- mean(y)

print(c(y_mean_msy,y_mean_0.5msy))
