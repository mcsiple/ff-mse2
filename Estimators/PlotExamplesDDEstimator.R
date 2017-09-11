source("~/Dropbox/Chapter4-HarvestControlRules/Code/ff-mse2/Estimators/Estimators.R")
par(mfrow=c(2,1))

# Sine curve example ------------------------------------------------------
tau0 <-  0.1
sigma0 <- 0.2
#tau1 <- (1/tau0^2+1/sigma0^2)^(-0.5)
amp <- 100
freq <- 15
# get a new random seed
rm(.Random.seed)
seed.index <- runif(1,1,100)
set.seed(seed.index)
timelist <- seq(1,60)
B <- amp + amp/1.3 * sin(pi*timelist/freq)
y <- rep(NA, length(timelist))
E <- rep (NA, length(timelist))
E[1]<-B[1]*exp(rnorm(1,0, sigma0))

for(i in 2:length(B)){
  E[i] <- tim.assessment(B = B[i],Eprev = E[i-1],sigma0 = sigma0,tau0 = tau0)
}
plot(B,type='l')
lines(E,col='blue')


#  Actual time series example ---------------------------------------------
B.true <- c(3392.50470568473, 3870.74716918686, 5262.84780497428, 6277.17522686664,  
            7842.31833734085, 4470.77382770654, 8900.63282345029, 5932.59637649889, 
            10353.7303951785, 4255.02095758172, 3245.76935065602, 5026.88222164352, 
            2796.22225220037, 1921.88146228901, 1653.07947724387, 904.507136682943, 
            1344.97037412906, 1643.19549534365, 1371.44393140291, 1736.12985165079, 
            3403.66071706085, 7366.05585022072, 9042.25994491636, 9006.45500582731, 
            3401.21455022866, 2276.45501036882, 2289.70996727027, 2719.47023746305, 
            2316.42653647645, 1900.51158054185, 3850.24566803188, 4319.37004493524, 
            4621.79148841807, 6437.07779077995, 2886.15536407626, 3800.35269399225, 
            2372.0504565025, 3156.81765700592, 2554.18726858172, 4356.81372251249, 
            5022.60022859825, 4548.56317976613, 5612.05880504775, 2855.28529302988, 
            2261.65686250796, 4772.99684749926, 3661.46656527244, 1643.56268533365, 
            1641.65435056853, 1672.37277047132, 2807.75463551049, 3347.00344357892, 
            4098.29010637602, 4922.94228974542, 8452.46928649048, 5793.51812263745, 
            4416.60642293191, 1892.68921288849, 2710.52490095019, 1754.64167878724
            )
#This is some example biomass from simulations
plot(timelist,B.true,type='l')
B.obs <- B.true

B.obs[1] <- B.true[1]*exp(rnorm(1,0, sigma0))

for(i in 2:length(B.obs)){
  B.obs[i] <- tim.assessment(B = B.true[i],Eprev = B.obs[i-1],sigma0 = sigma0,tau0 = tau0)
}

lines(B.obs,col='blue')

# Now compare to AC error:
B.obs2 <- B.obs
        sig.s = 0.41 # This value comes from running the delay function a bunch of times, getting a target sd(log) - this is the target for AC
                     # I changed this from the Wiedenmann et al. value because the error didn't look big enough
        rho = 0.5

# eps.prev = ifelse(yr==1,0.5,eps.prev) # Initialize epsilon. Need to do this every time this type of error is used. #2 is NEW
# outs <-  add.wied.error(biomass.true = oneplus.biomass[yr],
#                         epsilon.prev = eps.prev, 
#                         sig.s =  sig.s, rho = rho)
# biomass[yr] <- outs$biomass.est
# eps.prev <- outs$epsilon.curr

eps.vec <- vector(length=length(B.obs))
eps.vec[1] <- 0.5

for(i in 2:length(B.obs2)){
    outs <- add.wied.error(biomass.true = B.true[i],epsilon.prev = eps.vec[i-1],sig.s = sig.s,rho=rho)
    B.obs2[i] <- outs$biomass.est
    eps.vec[i] <-outs$epsilon.curr
}

lines(B.obs2,col='red')

