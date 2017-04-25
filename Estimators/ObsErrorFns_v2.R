#############################################################################
# These functions generate observation error within the model
# Most of these assume that whatever decision is made (via the HCR) is based on the biomass detected by the survey
# If you want to plot examples of all the errors, set plot.examples=TRUE

plot.examples = FALSE

add.LN.error <- function(biomass.true, obs.cv, years=1){
  #' @description adds random, lognormal error to biomass.true
  #' @param biomass.true - true biomass
  #' @param obs.cv - observation error
  #' @param years - number of years in time series. Automatically set to 1 for the operating model
  #' @return list: observed biomass and observation error in year y
  obs.error <- rlnorm(years, -obs.cv^2/2, obs.cv)
  biomass.obs <- biomass.true * obs.error
  return(list(biomass.obs = biomass.obs, error = obs.error))
}


# Autocorrelated error with a delay for big changes -----------------------
# This function has been temporarily decommissioned (8/12/16)

add.delay.error <- function(prev.yrs.b.vec, # Years of true biomass before biomass.true
                            biomass.true, rho, epsilon.prev, sig.s, # AC error parameters
                            drama.threshold, detect.delay){         # new parameters
  #' @description Calculates observed biomass based on an error that is lognormal and autocorrelated, but adds a delay anytime the true biomass exceeds threshold "drama.threshold"
  #' @param prev.yrs.b.vec - vector of true biomass up to the current year
  #' @param biomass.true - true biomass
  #' @param rho - amount of autocorrelation in error
  #' @param epsilon.prev - obs error in the previous year
  #' @param sig.s - sd of the deviations around epsilon
  #' @param drama.treshold - the -fold change from the previous year that qualifies a change as a "dramatic change." 
  #' @param detect.delay - number of years that it takes for management to detect a big change
  #' @return list: observed biomass and error in year y 
  nprevyrs <- length(prev.yrs.b.vec)
  dd <- detect.delay
  
  if(nprevyrs < detect.delay){ 
    print("Early in the biomass time series... filling in years up to  detection delay with regular AC error")
        # First, set up biomass in year dd+1
        curly.phi <- rnorm(1,0,sig.s) # the degree of autocorrelation in the estimates
        epsilon.curr <- rho * epsilon.prev + sqrt(1-(rho^2)) * curly.phi # epsilon values are the error
        firstBobs <- biomass.true*exp(epsilon.curr-(0.5*sig.s^2))
        
        return(list(biomass.est = firstBobs, epsilon.curr=epsilon.curr)) # If you are in the early parts of the time series, just use regular autocorrelated error. Then you're done!
        }
  
  else{
  prev.B <- tail(x = prev.yrs.b.vec) # true biomass one year ago
  if ( biomass.true > (1-drama.threshold) * prev.B  &&  biomass.true < (1+drama.threshold) * prev.B){
                         # If biomass.true is within thresholds
                         curly.phi <- rnorm(1,0,sig.s) # the degree of autocorrelation in the estimates
                         epsilon.curr <- rho * epsilon.prev + sqrt(1-(rho^2)) * curly.phi # epsilon values are the error
                         biomass.obs <- biomass.true*exp(epsilon.curr-(0.5*sig.s^2))
  }
    else {
      # If biomass is outside the threshold, base current value on the value dd years ago
      curly.phi <- rnorm(1,0,sig.s) # the degree of autocorrelation in the estimates
      epsilon.curr <- rho * epsilon.prev + sqrt(1-(rho^2)) * curly.phi # epsilon values are the error
      biomass.obs <- prev.yrs.b.vec[nprevyrs-(dd-1)] * exp(epsilon.curr-(0.5*sig.s^2))
    }
  
  return(list(biomass.est=biomass.obs,epsilon.curr = epsilon.curr))
  }
}

# Autocorrelated errors with parameters from Wiedenmann 2015 ------------------------------------------
add.wied.error <- function(biomass.true, epsilon.prev, sig.s, rho){
  #' @description adds autocorrelated obervation error with lag 1, as in Wiedenmann 2015
  #' @param biomass.true - the true biomass in year y
  #' @param epsilon.prev - observation error in previous year
  #' @param sig.s - sd of observation errors 
  #' @param rho - autocorrelation in obs error
  #' @return list: "observed" biomass, and epsilon (observation error) in the current year
  curly.phi <- rnorm(1,0,sig.s) # the degree of autocorrelation in the estimates
  epsilon.curr <- rho * epsilon.prev + sqrt(1-(rho^2)) * curly.phi # epsilon values are the error
  biomass.est <- biomass.true*exp(epsilon.curr-(0.5*sig.s^2))
  return(list(biomass.est=biomass.est,epsilon.curr=epsilon.curr))
}

# Autocorrelation error with worse estimation when there is a big change in biomass --------
error.nis <- function(biomass.true, epsilon.prev, sig.s, rho, biomass.prev = 1, phi.err = 1){
  #' @description This function is similar to Wiedenmann autocorrelated error, but adds error proportionally to the annual change in biomass-- if the biomass changes a lot, the error is large.
  #' @author M. Siple and N. Sand Jacobsen
  #' @param biomass.true - the true biomass in year y
  #' @param epsilon.prev - observation error in previous year
  #' @param sig.s - additional variation around observation errors 
  #' @param rho - autocorrelation in obs error
  #' @param biomass.prev - biomass in the previous year
  #' @param phi.err - how strongly the change in biomass affects the error
  #' @return list: "estimated" biomass, current amount of autocorrelated error, and the total observation error in year y (eps.curr)
  curly.phi <- rnorm(1,0,sig.s) # the degree of autocorrelation in the estimates
  epsilon.curr <- (rho * epsilon.prev + sqrt(1-(rho^2)) * curly.phi) # epsilon values are the error
  # Multiply by a function that takes the variation into account
  B.rel <- biomass.true/biomass.prev
  # obs error becomes equally large on log scales (reduction to 10% is the same as 10 times increase)
  # phi.err determines its strength
  eps.err <- 1+(phi.err*(log(B.rel)^2))   # If biomass[t-1] == biomass[t] the multiplier is 1 (log(1) = 0)
  epsilon.curr <- epsilon.curr*eps.err
  biomass.est <- biomass.true*exp(epsilon.curr-(0.5*sig.s^2)) 
  return(list(biomass.est=biomass.est,epsilon.curr=epsilon.curr,eps.err = eps.err))
}


# “Bayesian” version by Tim ------------------------------------------------
#' Commented out here because it is also in the file "TimAssessment.R"
#' tim.assessment <- function(Bprev,Bcurr,sigma_sq = 1.5, tau0 = 0.65){
#'   #' @description This function draws one value of log(B[t+1]/B[t]) from the posterior created by a normal prior for log(B[t+1]/B[t]) and a normal likelihood. The draw is one's "best estimate" of the true change in biomass from the previous year. Then it solves for what the "best estimate" of this year's biomass is, based on the posterior draw and Bprev, the biomass in the previous time step. 
#'   #' @details The goal of this type of observation error is to mimic what observations might be if we had some prior knowledge about the expected change in B. The level of confidence in our ability to detect change is set by tau0, and sigma^2, the variance in the likelihood (CHECK)
#'   #' @param Bprev - biomass in the previous year. This can be either the true biomass in the previous year, or the observed biomass in the previous year, if using a survey instead of true values.
#'   #' @param Bcurr - biomass for the current year. Like Bprev, this can be the observed biomass or the true biomass in the current year.
#'   #' @param sigma_sq - the variance in the PDF that is in the likelihood
#'   #' @param tau0 - the amount of confidence in the survey (or assessment)'s ability to detect changes in biomass. A higher tau0 reflects lower confidence in the ability of the survey or assessment to detect changes in B.
#'   #' @return A single value of "observed" biomass, derived from the prior, Bcurr, and Bprev, to be used with the HCR to determine F in the current year.
#'   
#'   Y_obs <- log(Bprev/Bcurr) #Obesrved changed in biomass
#'   mu_post = Y_obs * ( sigma_sq/((tau0^2) + 1)  )^-1
#'   tau1_squared = ( (1/tau0^2) + (1/sigma_sq) )^-1
#'   tau1 <- sqrt(tau1_squared)
#'   draw <- rnorm(1,mu_post, tau1)
#'   Obs <- exp(draw) * Bprev
#'   return(Obs)
#' }



# Test and plot all these obs error situations ----------------------------
###########################################################################
###########################################################################
# This code only runs if you set plot.examples==TRUE
if(plot.examples==TRUE){
                setwd("/Users/mcsiple/Dropbox/Chapter 4 - Harvest Control Rules/Figures")
              
              # simulate a time series to test ----------------------------
              nyears <- 40
              try.ts <- exp(arima.sim(n=nyears, model = list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488))))  
              
              # Print a figure showing all three error types
              setwd("/Users/mcsiple/Dropbox/Chapter 4 - Harvest Control Rules/Figures")
              pdf(file = "Error Types.pdf",width =6,height = 12,useDingbats = FALSE)
              par(mfrow=c(4,1))
              
              #  autocorrelated error  (Wiedenmann et al.) --------------------------------------------------
              
              obs.ts <- vector(length=length(try.ts))
              sig.s = 0.95 # 0.35 == This is what Wiedenmann et al. use in the paper; I think realistically it should  be higher.
              rho = 0.7
              eps.prev <- 0.5 # Initialize epsilon value
              for(i in 1:length(obs.ts)){
                eps.prev = ifelse(i==1,0.5,eps.prev) # If it's the first year of the time series, need to set epsilon.prev
                outs <- add.wied.error(biomass.true = try.ts[i],epsilon.prev = eps.prev, sig.s =  sig.s, rho = rho)
                obs.ts[i] <- outs$biomass.est
                eps.prev <- outs$epsilon.curr
              }
              
              par(mai=c(.5,1,1,.5))
              plot(try.ts,type='o',lwd=2, xlab = 'time', ylab = 'Biomass (simulated)',main=paste("Autocorrelated error \n","eps0=0.5","sig.s=",sig.s,"rho=",rho),ylim=c(0.1,300), log = 'y') #c(min(c(try.ts,obs.ts)),max(c(try.ts,obs.ts)))
              lines(1:length(try.ts),obs.ts,lty=3,lwd=2,col="darkgreen")
              legend('topright', lty=c(1,3),lwd=c(2,2),col=c("black","darkgreen"),legend = c('true B','obs B'))
              
              
              # Random lognormal  error -------------------------------------------------
              
              obs.ts <- vector(length=length(try.ts))
              obs.cv = 1.2
              for(i in 1:length(obs.ts)){
                obs.ts[i] <- add.LN.error(biomass.true = try.ts[i],obs.cv = obs.cv,years = 1)$biomass.obs
              }
              plot(try.ts,type='o',lwd=2, log='y', ylim= c(min(c(try.ts,obs.ts)),max(c(try.ts,obs.ts))),
                   xlab = 'time', ylab = 'Biomass (simulated)',main=paste("Lognormal observation error \n","obs.cv=",obs.cv))
              lines(1:length(try.ts),obs.ts,lty=3,lwd=2,col="blue")
              legend('topright', lty=c(1,3),lwd=c(2,2),col=c("black","blue"),legend = c('true B','obs B'))
              
              
              # Nis error ---------------------------------------------------------------
              
              # Plot the true and observed biomass together
              
              obs.ts <- vector(length=length(try.ts))
              obs.err <- vector(length=length(try.ts)) 
              obs.epsR<- vector(length=length(try.ts)) 
              
              rho = 0.7
              eps.prev <- 0.5 #Initialize epsilon value
              phi.err <- 0.9 # How bad is the observation error (0:1)
              
              obs.ts[1] <- error.nis(biomass.true = try.ts[1],epsilon.prev = eps.prev, sig.s =  0.35,rho = rho,
                                     biomass.prev = try.ts[1], phi.err = phi.err)$biomass.est # Perfect knowledge in the first year
              
              obs.err[1] <- error.nis(biomass.true = try.ts[1],epsilon.prev = eps.prev, sig.s =  0.35, rho = rho,
                                      biomass.prev = try.ts[1], phi.err = phi.err)$epsilon.curr
              obs.epsR[1] <- 1 # Assuming starting in equilibrium
              
              for(i in 2:length(obs.ts)){
                ls.tmp <- error.nis(biomass.true = try.ts[i],epsilon.prev = eps.prev, sig.s =  0.35, rho = rho,
                                    biomass.prev = try.ts[i-1], phi.err = phi.err)
                obs.ts[i] <- ls.tmp$biomass.est
                obs.err[i] <- ls.tmp$epsilon.curr
                obs.epsR[i]<- ls.tmp$eps.err
              }
              
              # Biomass 
              plot(try.ts,type='o',lwd=2, xlab = 'time', ylab = 'Biomass (simulated)', 
                   ylim = c(min(c(try.ts,obs.ts)),max(c(try.ts,obs.ts))), log = 'y',
                   main=paste("AC error proportional to change in B \n", "eps0=0.5","sig.s=",sig.s,"rho=",rho))
              lines(1:length(try.ts),obs.ts,lty=3,lwd=2,col="red")
              legend('topright', lty=c(1,3),lwd=c(2,2),col=c("black","red"),legend = c('true B','obs B'))
              
              
              
              
              # autocorrelated error with a delay for large changes ---------------------
              
              obs.ts <- vector(length=length(try.ts))
              sig.s = 0.85 # 0.35 == This is what Wiedenmann et al. use in the paper; I think realistically it should  be higher.
              rho = 0.7
              eps.prev <- 0.5 # Initialize epsilon value
              drama.threshold <- 0.75
              detect.delay <- 4
              
              for(i in 1:(length(obs.ts)-1)){
                eps.prev = ifelse(i==1,0.5,eps.prev) # If it's the first year of the time series, need to set epsilon.prev
                outs <- add.delay.error(prev.yrs.b.vec = try.ts[1:i],biomass.true = try.ts[i+1],epsilon.prev = eps.prev, sig.s =  sig.s, rho = rho, drama.threshold = drama.threshold, detect.delay = detect.delay)
                obs.ts[i] <- outs$biomass.est
                eps.prev <- outs$epsilon.curr
              }
              plot(try.ts,type='o',lwd=2, xlab = 'time', ylab = 'Biomass (simulated)',main=paste("AC error with a delay for big peaks \n","sig.s=",sig.s,"rho=",rho,"detect.delay=",detect.delay,"yrs", "\n drama.threshold=",drama.threshold),log='y')
              lines(1:length(try.ts),obs.ts,lty=3,lwd=2,col="purple")
              legend('topright', lty=c(1,3),lwd=c(2,2),col=c("black","purple"),legend = c('true B','obs B'))
              
              
              # Tim assessment error
              dev.off()
}




