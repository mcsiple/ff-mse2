# File to test the error types, sensitivity, etc. 
basedir <- "/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Code/ff-mse2"
resultsdir <- "/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Results"
par(mfcol=c(2,3))
Types <- c("Sardine","Anchovy","Menhaden")
for(t in 1:3){
  subDir <- Types[t] # Name of ff type

    # Load reshape2 and ggplot2 for plotting examples
    library(reshape2)
    library(ggplot2)
    # Load rev devs generating fxn, MSE main model, estimator fxns
    toplot=FALSE
    source(file.path(basedir,"Recruitment/GenerateDevs.R")) # Don't plot examples of rec trajectories
    source(file.path("/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Code/HCR_Trajectory.R"))
    source(file.path(basedir,"Estimators/Estimators.R"))
    
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
    
    
    # Load harvest rules
    source(file.path(basedir,"Control Rules/smith_oceana.R"))
    source(file.path(basedir,"Control Rules/cfp.R"))
    source(file.path(basedir,"Control Rules/hockey-stick.R"))
    source(file.path(basedir,"Control Rules/trend-based-rule.R"))
    
    # Other params
    years.test <- 250
    obs.type <- "AC"
    steepness  <- 0.6
    set.seed(123)
    rec.dev.test <- generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd)
    (equilib <- getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=steepness) )
    testie <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "constF",const.f.rate=0,equilib = equilib,steepness=steepness,obs.type = obs.type,R0.traj = NA, tim.params = tim.params)
    
    plot(testie$oneplus.biomass[100:250],type='l',#ylim=c(0,2e6),
         col="darkblue",lwd=1.8,
         xlab="Year",
         ylab="Biomass or B_est",
         main=paste(obs.type))
    lines(testie$biomass[100:250],col='lightblue',lwd=1.8) #obs 1+ biomass
    # lines(testie$total.catch,col='red')
    obs.type = "Tim"
    set.seed(123)
    tim.params <- list(sigma=1.2,tau0=1.2/5)
    rec.dev.test <- generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd)
    (equilib <- getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=steepness) )
    testie <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "constF",const.f.rate=0,equilib = equilib,steepness=steepness,obs.type = obs.type,R0.traj = NA, tim.params = tim.params)
    plot(testie$oneplus.biomass[100:250],type='l',#ylim=c(0,2e6),
         col="darkblue",lwd=1.8,
         xlab="Year",
         ylab="Biomass or B_est",
         main=paste(obs.type))
    lines(testie$biomass[100:250],col='lightblue',lwd=1.8)
}
# Get sigma from the variability in the delay detection model
nreps <- 100
sd.vec <- 1:nreps
for(i in 1:nreps){
  rec.dev.test <- generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd)
  (equilib <- getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=steepness) )
  testie <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, 
                            rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "constF",
                            const.f.rate=0,equilib = equilib,steepness=steepness,obs.type = obs.type,R0.traj = NA, tim.params = tim.params)
  sd.vec[i] <- sd(log(testie$biomass))
}
sigma.vec <- sd.vec*sqrt(1-0.8^2)
test.sigma <- median(sigma.vec)
#        
#         
