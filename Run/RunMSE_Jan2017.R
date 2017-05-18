# Set directories
basedir <- "/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Code/ff-mse2"
resultsdir <- "/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Results"
subDir <- "Sardine" # Name of ff type

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

# If a sardine folder doesn't exist already, create one!
dir.create(file.path(resultsdir, subDir))
setwd(file.path(resultsdir, subDir))

h = c(0.9, 0.6)
obs.error.type=c("AC","Tim")
HCR <- c("cfp","constF","oceana","lenfest","trend")

scenarios <- expand.grid(h,obs.error.type,HCR,recruit.sd,recruit.rho)
colnames(scenarios) <- c("h","obs.error.type","HCR","recruit.sd","recruit.rho")
nscenarios <- nrow(scenarios)
scenarios$scenario <- 1:nscenarios #Label them so it's easier to find/index em later

#Set up other simulation params
#nyrs.to.use = 150
years.test = 250
nsims = 1000
tim.params = list(sigma = 1.2,tau0=1.2/5) # The value for tau comes from the "sensitivity tests" for the error function-- this sig/tau ratio ensures that there won't be too many giant peaks (>5 x the max biomass)
                                      # sigma = 1.2 and tau0 = 2 means this is the "high-tau" scenario. If tau 
const.f.rate = 0.6
R0.sens = NA #NO DYNAMIC R0 in these cases!!

write.table(scenarios,file = "Scenario_Table.txt")

CFP <- oceana <- constF <- lenfest <- trend <-
            list(biomass=matrix(nrow = nsims,ncol = years.test),
            total.catch=matrix(nrow = nsims,ncol = years.test),
            fishing= matrix(nrow = nsims,ncol = years.test),
            rec= matrix(nrow = nsims,ncol = years.test),
            depl= matrix(nrow = nsims,ncol = years.test),
            obs.biomass = matrix(nrow = nsims,ncol = years.test),
            total.true.biomass = matrix(nrow = nsims,ncol = years.test))


# Test params and runs to make sure they look good ------------------------
        steepness = scenarios$h[1]
        obs.type <- scenarios$obs.error.type[2]
        HCR <- scenarios$HCR[1]
        recruit.sd = .6 #scenarios$recruit.sd[1]
        recruit.rho = .9 #scenarios$recruit.rho[1]
        equilib = getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=steepness)
        rec.dev.test <- generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd)
        test.constF <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "constF", const.f.rate = 0, steepness = steepness,obs.type = obs.type,equilib=equilib,R0.traj = R0.sens, tim.params = tim.params)
        nofish <- melt(test.constF[-c(1,4,9,10)])
        nofish$year <- rep(1:years.test,times=length(unique(nofish$L1)))
        #colnames(nofish) <- c("age","year","value","variable")
        ggplot(nofish,aes(x=year,y=value)) + geom_line() + facet_wrap(~L1,scales = "free_y") + xlim(c(150,250))

        #plot(test.constF$biomass,test.constF$oneplus.biomass)
# --- ---------------------------------------------------------------------


for(s in 1:nscenarios){
  steepness = scenarios$h[s]
  obs.type <- scenarios$obs.error.type[s]
  HCR <- scenarios$HCR[s]
  recruit.sd = scenarios$recruit.sd[s]
  recruit.rho = scenarios$recruit.rho[s]
  # results.cfp <- results.constF <- results.oceana <- list()
  # Still need to fix recruitment deviations function -- ???? is this still true??
  equilib = getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=steepness) # NO recruitment devs used in the equilibrium calculations, so don't need to embed in the loop
  
  if(HCR=="cfp"){    
  set.seed(123) # IMPORTANT! Start each round of sims at the same random seed
  for (sim in 1:nsims){
        rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
        expt.cfp <- calc.trajectory(lh = lh.test,obs.cv = NA, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "cfp",equilib = equilib,steepness=steepness,obs.type = obs.type,R0.traj = R0.sens, tim.params = tim.params)
        
    CFP[["biomass"]][sim,] <- expt.cfp$oneplus.biomass
    CFP[["total.catch"]][sim,] <- expt.cfp$total.catch
    CFP[["fishing"]][sim,] <- expt.cfp$fishing
    CFP[["rec"]][sim,] <- expt.cfp$rec
    CFP[["depl"]][sim,] <- expt.cfp$depl
    CFP[["obs.biomass"]][sim,] <- expt.cfp$biomass
    CFP[["total.true.biomass"]][sim,] <- expt.cfp$biomass.true
  }
    save(CFP,file=paste("All",s,"CFP",".RData",sep="_"))
    }
    
  if(HCR=="constF"){
    set.seed(123) # same seed
    for (sim in 1:nsims){
          rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
          expt.constF <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "constF", const.f.rate = const.f.rate, steepness = steepness,obs.type = obs.type,equilib=equilib,R0.traj = R0.sens, tim.params = tim.params)
        
    constF[["biomass"]][sim,] <- expt.constF$oneplus.biomass
    constF[["total.catch"]][sim,] <- expt.constF$total.catch
    constF[["fishing"]][sim,] <- expt.constF$fishing
    constF[["rec"]][sim,] <- expt.constF$rec
    constF[["depl"]][sim,] <- expt.constF$depl
    constF[["obs.biomass"]][sim,] <- expt.constF$biomass
    constF[["total.true.biomass"]][sim,] <- expt.constF$biomass.true
    }
    save(constF,file=paste("All",s,"constF",".RData",sep="_"))
    }
    
  if(HCR=="oceana"){
    set.seed(123) # same seed
    for (sim in 1:nsims){
        rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
        expt.oceana <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "oceana",equilib = equilib,steepness=steepness,obs.type = obs.type,R0.traj = R0.sens, tim.params = tim.params)
       
    oceana[["biomass"]][sim,] <- expt.oceana$oneplus.biomass
    oceana[["total.catch"]][sim,] <- expt.oceana$total.catch
    oceana[["fishing"]][sim,] <- expt.oceana$fishing
    oceana[["rec"]][sim,] <- expt.oceana$rec
    oceana[["depl"]][sim,] <- expt.oceana$depl
    oceana[["obs.biomass"]][sim,] <- expt.oceana$biomass
    oceana[["total.true.biomass"]][sim,] <- expt.oceana$biomass.true
    }
    save(oceana,file=paste("All",s,"oceana",".RData",sep="_")) 
  }
  
  if(HCR=="lenfest"){
    set.seed(123) # same seed
    for (sim in 1:nsims){
      rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
      expt.lenfest <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "lenfest",equilib = equilib,steepness=steepness,obs.type = obs.type,R0.traj = R0.sens, tim.params = tim.params)
      
      lenfest[["biomass"]][sim,] <- expt.lenfest$oneplus.biomass
      lenfest[["total.catch"]][sim,] <- expt.lenfest$total.catch
      lenfest[["fishing"]][sim,] <- expt.lenfest$fishing
      lenfest[["rec"]][sim,] <- expt.lenfest$rec
      lenfest[["depl"]][sim,] <- expt.lenfest$depl
      lenfest[["obs.biomass"]][sim,] <- expt.lenfest$biomass
      lenfest[["total.true.biomass"]][sim,] <- expt.lenfest$biomass.true
    }
    save(lenfest,file=paste("All",s,"lenfest",".RData",sep="_")) 
  }
  
  if(HCR=="trend"){
    set.seed(123) # same seed
    for (sim in 1:nsims){
      rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
      expt.trend <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "trend",const.f.rate = 0.6,equilib = equilib,steepness=steepness,obs.type = obs.type,R0.traj = R0.sens, tim.params = tim.params)
      
      trend[["biomass"]][sim,] <- expt.trend$oneplus.biomass
      trend[["total.catch"]][sim,] <- expt.trend$total.catch
      trend[["fishing"]][sim,] <- expt.trend$fishing
      trend[["rec"]][sim,] <- expt.trend$rec
      trend[["depl"]][sim,] <- expt.trend$depl
      trend[["obs.biomass"]][sim,] <- expt.trend$biomass # Observed one-plus biomass
      trend[["total.true.biomass"]][sim,] <- expt.trend$biomass.true
    }
    save(trend,file=paste("All",s,"trend",".RData",sep="_")) 
  }
}


  
  