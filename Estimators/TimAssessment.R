# The "Delayed Detection" estimator. This function acts like there is a prior on the estimated B, based on the B last year. So it does not "expect" big interannual changes in biomass. 
# Tim's brilliant plan for mimicking obervation errors using something that seems like a Bayesian stock assessment
# 
# source("/Users/mcsiple/Dropbox/Chapter 4 - Harvest Control Rules/Code/Estimators/ObsErrorFns_v2.R")


tim.assessment <- function(Bprev,Bcurr,sigma = 1.5, tau0 = 0.65){
  #' @description This function draws one value of log(B[t+1]/B[t]) from the posterior created by a normal prior for log(B[t+1]/B[t]) and a normal likelihood. The draw is one's "best estimate" of the true change in biomass from the previous year. Then it solves for what the "best estimate" of this year's biomass is, based on the posterior draw and Bprev, the biomass in the previous time step. 
  #' @details The goal of this type of observation error is to mimic what observations might be if we had some prior knowledge about the expected change in B. The level of confidence in our ability to detect change is set by tau0, and sigma^2, the variance in the likelihood (CHECK)
  #' @param Bprev - biomass in the previous year. This can be either the true biomass in the previous year, or the observed biomass in the previous year, if using a survey instead of true values.
  #' @param Bcurr - biomass for the current year. Like Bprev, this can be the observed biomass or the true biomass in the current year.
  #' @param sigma_sq - the variance in the PDF that is in the likelihood
  #' @param tau0 - the amount of confidence in the survey (or assessment)'s ability to detect changes in biomass. A higher tau0 reflects lower confidence in the ability of the survey or assessment to detect changes in B.
  #' @return A single value of "observed" biomass, derived from the prior, Bcurr, and Bprev, to be used with the HCR to determine F in the current year.
  
  Y_obs <- log(Bprev/Bcurr) #Obesrved changed in biomass
  sigma_sq <- sigma^2
  mu_post = Y_obs * ( sigma_sq/((tau0^2) + 1)  )^-1
  tau1_squared = ( (1/tau0^2) + (1/sigma_sq) )^-1
  tau1 <- sqrt(tau1_squared)
  draw <- rnorm(1,mu_post, tau1)
  Obs <- exp(draw) * Bprev
  return(Obs)

}


#EVERYTHING BELOW IS FOR TESTING
# Observation is a sample from the posterior distribution
#nyears=500
# Loop thru values
                  #sigma_sq = 1.5
                  #tau0 = .65
                  Y_obs = 1000 # This is whatever the observation is (can be random or based on the autocorrelated error from before)
                  # Posterior for mu, the true mean change in biomass
                  nyears=24
                  arima.ts <- exp(arima.sim(n=nyears, model = list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488))))  
                  single.shift <- rep(c(25,100),each=nyears/2)
				pdf(file="TimError_example.pdf",height = 6,width=10)
                  par(mfrow=c(2,3))
                  sigtau.vec <- c(20,7,1)
                  for(type in 1:2){ # type 1: single shift 
                    # type 2: arima time series
                  for (st in 1:length(sigtau.vec)){
                  
                  # Params for Bayesian error
                  sig_tau_ratio = sigtau.vec[st]
                  if(type==1) {try.ts = single.shift}
                  if (type==2) {try.ts = arima.ts}
                  sigma = 1.2
                  tau0 <- sigma / sig_tau_ratio
                  
                  # Params for observation error
                  sig.s = 0.02
                  rho = 0.2
    
                  # 
                  est.ts <- obs.ts <- eps.ts <-  vector()
                  est.ts[1] <- obs.ts[1] <- try.ts[1] #Perfect obs in first year
                  eps.ts[1] <- 1 # No change in first year
                  # 
                  for(i in 2:nyears){
                    observation <- add.wied.error(biomass.true = try.ts[i],epsilon.prev = eps.ts[i-1],sig.s = sig.s, rho = rho)
                    obs.ts[i] <- observation$biomass.est
                    eps.ts[i] <- observation$epsilon.curr
                    est.ts[i] <- tim.assessment(Bprev=obs.ts[i-1],Bcurr=obs.ts[i],tau0=tau0,sigma = sigma)
                  }
                  # 
                  obs.color <- adjustcolor( "blue", alpha.f = 0.5)
                  est.color <- adjustcolor( "red", alpha.f = 0.5)
                  yrange=range(c(obs.ts,est.ts,eps.ts))
                  #       #pdf("TimError.pdf",width = 6,height = 5)
                          plot(1:nyears,try.ts,type='o',lwd=2,ylim=yrange,
                               #xaxt="n",yaxt="n",bty="n",bty="n", # these are for making the presentation
                               xlab="Year",ylab="Biomass")
                          lines(1:nyears,obs.ts,col=obs.color,lwd=2)
                          lines(1:nyears,est.ts,col=est.color,lwd=2,type='o')
                          text(5,150,labels = paste("sig/tau =",sig_tau_ratio,"\n sigma =",sigma,"\n tau = ",tau0,"\n obs rho =", rho, "\n obs sigma.s =", sig.s))

                          legend("topleft",lty=c(1,1,1),col=c("black",obs.color,est.color),lwd=c(2,2,2),legend = c("True B","Obs B","Est'd B"),pch=c(1,NA,1))
                  }}
dev.off()


# Another section to look at how close the survey is to the true b --------
# also using sigtau.vec
try.ts <- testie$oneplus.biomass
sig.tau.vec <- seq(1,20,by=1)

# 
    sig_tau_ratio = sigtau.vec[st]
    sigma = 1.2
    tau0 <- sigma / sig_tau_ratio
    est.ts <- obs.ts <- eps.ts <-  vector()
    est.ts[1] <- obs.ts[1] <- try.ts[1] #Perfect obs in first year
    eps.ts[1] <- 1 # No change in first year
    # 
    for(i in 2:nyears){
      observation <- add.wied.error(biomass.true = try.ts[i],epsilon.prev = eps.ts[i-1],sig.s = sig.s, rho = rho)
      obs.ts[i] <- observation$biomass.est
      eps.ts[i] <- observation$epsilon.curr
      est.ts[i] <- tim.assessment(Bprev=obs.ts[i-1],Bcurr=obs.ts[i],tau0=tau0,sigma = sigma)
    }
    # 
    count.huge.peaks <- length(which(est.ts > 5*try.ts))

    