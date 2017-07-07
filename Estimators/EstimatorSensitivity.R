# File to test the error types, sensitivity, etc. 
basedir <- "/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Code/ff-mse2"
resultsdir <- "/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Results"


# Load harvest rules
source(file.path(basedir,"Control Rules/smith_oceana.R"))
source(file.path(basedir,"Control Rules/cfp.R"))
source(file.path(basedir,"Control Rules/hockey-stick.R"))
source(file.path(basedir,"Control Rules/trend-based-rule.R"))

# Load reshape2 and ggplot2 for plotting examples
library(reshape2)
library(ggplot2)

# Load rev devs generating fxn, MSE main model, estimator fxns
toplot=FALSE
source(file.path(basedir,"Recruitment/GenerateDevs.R")) # Don't plot examples of rec trajectories
source(file.path("/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Code/HCR_Trajectory.R"))
source(file.path(basedir,"Estimators/Estimators.R"))

# 
Types <- c("Sardine","Anchovy","Menhaden")
nyrs <- 20
startyr <- 250-nyrs
#yrs <- startyr:250
yrs=1:20
par(mfcol=c(2,3))

for(t in 1:3){
  subDir <- Types[t] # Name of ff type
    
    # Load control files
    # Sardines
    if(subDir == "Sardine"){
      source(file.path(basedir,"Ctl/Sardine_LHControl.R"))
      source(file.path(basedir,"Ctl/Sardine_FisheryControl.R"))
      # Sardine params
      recruit.sd <- 0.6
      recruit.rho <- 0.9
    }
    # Anchovy/Herring
    if(subDir == "Anchovy"){
      source(file.path(basedir,"Ctl/Anchovy_LHControl.R"))
      source(file.path(basedir,"Ctl/Anchovy_FisheryControl.R"))
      #Anchovy recruitment dev params
      recruit.sd <- 0.6
      recruit.rho <- 0.5
    }
    
    # Menhaden
    if(subDir == "Menhaden"){
      source(file.path(basedir,"Ctl/Menhaden_LHControl.R"))
      source(file.path(basedir,"Ctl/Menhaden_FisheryControl.R"))
      #Anchovy recruitment dev params
      recruit.sd <- 0.8
      recruit.rho <- 0.2
    }
    
    # Other params
    years.test <- 250
    obs.type <- "AC"
    steepness  <- 0.6
    tim.params <- list(sigma0 = 0.2, tau0 = 0.1)
    set.seed(123)
    rec.dev.test <- generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd)
    (equilib <- getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=steepness) )
    testie <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "constF",const.f.rate=0,equilib = equilib,steepness=steepness,obs.type = obs.type,R0.traj = NA, tim.params = tim.params)
    yrange <- range(c(testie$biomass.oneplus.obs[yrs],testie$biomass.oneplus.true[yrs]))
    plot(testie$biomass.oneplus.true[yrs],type='l',ylim=yrange,
         col="darkblue",lwd=1.8,
         xlab="Year",
         ylab="Biomass or B_est",
         main=paste(obs.type, subDir, sep = " - "))
    lines(testie$biomass.oneplus.obs[yrs],col='lightblue',lwd=1.8) #obs 1+ biomass
    # lines(testie$total.catch,col='red')
    obs.type = "Tim"
    set.seed(123)
    #tim.params <- list(sigma=1.2,tau0=1.2/5)
    rec.dev.test <- generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd)
    (equilib <- getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=steepness) )
    testie <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "constF",const.f.rate=0,equilib = equilib,steepness=steepness,obs.type = obs.type,R0.traj = NA, tim.params = tim.params)
    yrange <- range(c(testie$biomass.oneplus.obs[yrs],testie$biomass.oneplus.true[yrs]))
    plot(testie$biomass.oneplus.true[yrs],type='l',ylim=yrange,
         col="darkblue",lwd=1.8,
         xlab="Year",
         ylab="Biomass or B_est",
         main=paste("DD", subDir, sep = " - "))
    lines(testie$biomass.oneplus.obs[yrs],col='lightblue',lwd=1.8)
}


# 
# testie <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = 20,hcr.type = "constF",const.f.rate=0,equilib = equilib,steepness=steepness,obs.type = obs.type,R0.traj = NA, tim.params = tim.params)
# testie$biomass.oneplus.obs
# testie$biomass.oneplus.true
# plot(testie$biomass.oneplus.true,type='l')
# lines(testie$biomass.oneplus.obs,col="lightblue")


# Get sigma from the variability in the delay detection model - sigma is different for each forage fish type
# The code below runs the DD ("Tim") version of the operating model, gets the SD from that, and compares it to the actual SD when there is just autocorrelated error.

tim.params
obs.type = "Tim"
target.sd.vec <- AC.sd.vec <- vector(length=3)
#     Type     TargetSD
# 1  Sardine 0.5117551
# 2  Anchovy 0.3221512
# 3 Menhaden 0.3009440
for(t in 1:3){
      type <- Types [t]
      if(type == "Sardine"){
        source(file.path(basedir,"Ctl/Sardine_LHControl.R"))
        source(file.path(basedir,"Ctl/Sardine_FisheryControl.R"))
        recruit.sd <- 0.6
        recruit.rho <- 0.9
        #sig.s <- 0.511
      }
      if(type == "Anchovy"){
        source(file.path(basedir,"Ctl/Anchovy_LHControl.R"))
        source(file.path(basedir,"Ctl/Anchovy_FisheryControl.R"))
        recruit.sd <- 0.6
        recruit.rho <- 0.5
        #sig.s <- 0.322
      }
      if(type == "Menhaden"){
        source(file.path(basedir,"Ctl/Menhaden_LHControl.R"))
        source(file.path(basedir,"Ctl/Menhaden_FisheryControl.R"))
        recruit.sd <- 0.8
        recruit.rho <- 0.2
        #sig.s <- 0.300
      }
      nreps <- 100
      sd.vec <- sd.vec2 <- 1:nreps
      set.seed(123)
      #sig.s <- (tim.params$sigma0^2 + 1/(tim.params$tau0^2))^0.5
      sig.s <- (1/(tim.params$tau0^2) + 1/(tim.params$sigma0^2))^(-0.5)
      rec.dev.test <- generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd)
      for(i in 1:nreps){
        (equilib <- getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=steepness) )
        testie <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, 
                                  rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "constF",
                                  const.f.rate=0,equilib = equilib,steepness=steepness,obs.type = "Tim",R0.traj = NA, tim.params = tim.params,sig.s=sig.s)
        testie2 <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, 
                                  rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "constF",
                                  const.f.rate=0,equilib = equilib,steepness=steepness,obs.type = "AC",R0.traj = NA, tim.params = tim.params,sig.s=sig.s)
        sd.vec[i] <- sd(log(testie$biomass.oneplus.obs))
        sd.vec2[i] <- sd(log(testie2$biomass.oneplus.obs))
      }
      sigma.vec <- sd.vec*sqrt(1-0.5^2)         # Corrections for autocorrelation? Because rho = 0.5
      sigma.vec2 <- sd.vec2*sqrt(1-0.5^2)
      test.sigma <- median(sigma.vec)
      test.sigma2 <- median(sigma.vec2)
      target.sd.vec[t] <- test.sigma
      AC.sd.vec[t] <- test.sigma2
      #sd.check.vec
}

#SDs of the simulated time series are very similar
sd(sd.vec)
sd(sd.vec2)

hist(sd.vec,xlim=c(0,0.75),col=rgb(1,0,0,0.5)) # DD = red; AC = blue
hist(sd.vec2, col=rgb(0,0,1,0.5), add=T)
print(data.frame("Type"=Types, "TargetSDlog" = target.sd.vec, "ActualSDlog"=AC.sd.vec))

#        
#         
