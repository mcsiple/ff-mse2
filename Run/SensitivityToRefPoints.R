# Sensitivity to differences in reference points
# The estimator within the model is not a full assessment, so it does not return reference points every iteration
# We can show the sensitivity of the results to over- or under-estimation of reference points:
# This code:
# - chooses a normal distribution around the "true" Fmsy, Bmsy, B0, which are used in the harvest control rules
# - the error in that distribution reflects the error in the estimators 
# - selects random values from that distribution
# - determines the variation in performance measures and/or tradeoffs, based on those runs


# Libraries
library(plyr)

types <- c("Anchovy"="Anchovy",
           "Menhaden" = "Menhaden",
           "Sardine" = "Sardine")

fftype = types[1] # Test with anchovy

# Set directories
basedir <- "~/Dropbox/Chapter4-HarvestControlRules/Code/ff-mse2"
resultsdir <- "~/Dropbox/Chapter4-HarvestControlRules/Results/Sensitivity"
subDir <- fftype

#Set up other simulation params
years.test = 250
nsims = 1000
tim.params = list(sigma0 = 0.2,tau0 = 0.1)
sig.s = 0.3
R0.sens = NA # NO DYNAMIC R0 anymore-- ignore

# Load packages
library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)

# Load rev devs generating fxn, MSE main model, estimator fxns
#toplot=FALSE      # Don't plot examples of recruitement trajectories
source(file.path(basedir,"Recruitment/GenerateDevs.R")) 
source(file.path(basedir,"Estimators/CalcFTrue.R"))
source(file.path(basedir,"Run/HCR_Trajectory_NEW.R"))
source(file.path(basedir,"Estimators/Estimators.R"))
source(file.path(basedir,"Run/generate_M.R"))

# Load control files & set parameter values
# Anchovy/Herring
if(subDir == "Anchovy"){
  source(file.path(basedir,"Ctl/Anchovy_LHControl.R"))
  source(file.path(basedir,"Ctl/Anchovy_FisheryControl.R"))
  #Anchovy recruitment dev params
  recruit.sd <- 0.6
  recruit.rho <- 0.5
}

# Load harvest rules
source(file.path(basedir,"Control Rules/smith_oceana.R"))
source(file.path(basedir,"Control Rules/cfp.R"))
source(file.path(basedir,"Control Rules/hockey-stick.R"))
source(file.path(basedir,"Control Rules/trend-based-rule.R"))


# Scenarios - for sensitivity, limiting the number of scenarios.

setwd(resultsdir)
h = c(0.6) #0.9, 
obs.error.type = c("AC") #,"noerror","Tim"
HCR = c("cfp","constF","C1","C2","C3","constF_HI") 
M.type = c("constant") # took out "regimeshift" and "time-varying" to save time but can be added back in for sensitivity analyses

scenarios <- expand.grid(h,obs.error.type,HCR,recruit.sd,recruit.rho,M.type)
colnames(scenarios) <- c("h","obs.error.type","HCR","recruit.sd","recruit.rho","M.type")
nscenarios <- nrow(scenarios)
scenarios$scenario <- 1:nscenarios #Label them so it's easier to find/index em later


write.table(scenarios,file = "Scenario_Table.txt")

CFP <- C1 <- C2 <- C3 <- constF <- trend <- constF_HI <- 
  list(biomass.oneplus.true=matrix(nrow = nsims,ncol = years.test), 
       total.catch=matrix(nrow = nsims,ncol = years.test),
       fishing= matrix(nrow = nsims,ncol = years.test),
       intended.f=matrix(nrow = nsims,ncol = years.test),
       rec= matrix(nrow = nsims,ncol = years.test),
       biomass.oneplus.obs = matrix(nrow = nsims,ncol = years.test),
       biomass.total.true = matrix(nrow = nsims,ncol = years.test),
       no.fishing.tb = matrix(nrow = nsims,ncol = years.test))

# Test params and runs to make sure they look good ------------------------
steepness = scenarios$h[1]
obs.type <- scenarios$obs.error.type[1]
HCR <- scenarios$HCR[1]
equilib = getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=steepness)
rec.dev.test <- generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd)
test.constF <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test,rec.ram = NA, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "constF", const.f.rate = 0, steepness = steepness,obs.type = obs.type,equilib=equilib,R0.traj = R0.sens, tim.params = tim.params,time.var.m = NA,sig.s = sig.s)
nofish <- melt(test.constF[c("biomass.oneplus.obs","biomass.total.obs","total.catch","fishing","biomass.total.true","rec","biomass.oneplus.true")])
nofish$year <- rep(1:years.test,times=length(unique(nofish$L1)))
ggplot(nofish,aes(x=year,y=value)) + geom_line() + facet_wrap(~L1,scales = "free_y") #+ xlim(c(150,250))



# Get appropriate SD for variation in F -----------------------------------
nsims = 100
sd.sim = vector(length = nsims)
set.seed(123)
for(sim in 1:nsims){
  rec.dev.test <- generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd)
  test.F.sd <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test,rec.ram = NA, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "constF", const.f.rate = 0.7, steepness = steepness,obs.type = "AC",equilib=equilib,R0.traj = R0.sens, tim.params = tim.params,time.var.m = NA,sig.s = sig.s)
  sd.sim[sim] <- sd(test.F.sd$fishing[150:250]) #last 100 years of time series
}
median(sd.sim)

# With sig.s = 0.3 and F = 0.7, this is the sd(fishing mortality) that you would expect with AC error
# F = 1.1    0.5602482
# F = 0.7    0.2956493
# F = 0.3   0.1047076
# F = 0.1   0.03187569
# So sd of F depends on the constant F rate that you set... I think there's an adjustment but I'm not sure what it is... might as well set sd high and see what happens.
#Going to use the high value of sd(F) just to be precautionary. sd = 0.56

Fmsy.sd <- 0.56
B0.sd <- 0.3

# Plot a normal curve 

# Test sensitivity to over-or under-estimating B0 ---------------------------------------------------

for(s in 1:nscenarios){  #
  steepness = scenarios$h[s]
  obs.type <- scenarios$obs.error.type[s]
  HCR <- scenarios$HCR[s]
  recruit.sd = scenarios$recruit.sd[s]
  recruit.rho = scenarios$recruit.rho[s]
  M.type = scenarios$M.type[s]
  equilib.true = getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=steepness) # NO recruitment devs used in the equilibrium calculations, so don't need to embed in the loop
  equilib <- equilib.true
  no.fishing <- matrix(NA, nrow = nsims, ncol = years.test)
  set.seed(123) # Start each round of sims at same random seed
  time.var.m <- NA # Base case: M constant 
  # if(M.type == "timevar"){time.var.m <- rw.M(Mbar = lh.test$M, rho.m = 0.6, sigma.m = 0.2,n = years.test)}
  # if(M.type == "regimeshift"){time.var.m <- regime.M(Mbar = lh.test$M,cutoff.yr = 201,n = years.test)}
  for (sim in 1:nsims){

    rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
    F0.Type <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "constF", const.f.rate = 0, steepness = steepness,obs.type = obs.type,equilib=equilib,R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s)$biomass.total.true # Only need to do this 1x for each simulation (not repeat for each CR) because the seed is the same and there is no fishing.
    no.fishing[sim,] <- F0.Type # This is the time series 
  }
  
  B0.error <- rnorm(1000,0,sd = B0.sd)
  const.f.rate = 0.5*equilib$Fmsy
  
  if(HCR=="cfp"){    
    set.seed(123) # Start each round of sims at same random seed
    for (sim in 1:nsims){
      equilib$B0 <- equilib.true$B0 * exp(B0.error[sim]) # error in B0 will only affect control rules. Error is added here (instead of before with the unfished part) but want unfished dynamics determined with the true equilib #s 
      rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
      expt.cfp <- calc.trajectory(lh = lh.test,obs.cv = NA, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "cfp",equilib = equilib,steepness=steepness,obs.type = obs.type,R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s,rec.ram = NA)
      
      CFP[["biomass.oneplus.true"]][sim,] <- expt.cfp$biomass.oneplus.true # This is the true one-plus biomass
      CFP[["total.catch"]][sim,] <- expt.cfp$total.catch  
      CFP[["fishing"]][sim,] <- expt.cfp$fishing
      CFP[["intended.f"]][sim,] <- expt.cfp$intended.f    
      CFP[["rec"]][sim,] <- expt.cfp$rec
      CFP[["biomass.oneplus.obs"]][sim,] <- expt.cfp$biomass.oneplus.obs        # This is the observed one-plus biomass
      CFP[["biomass.total.true"]][sim,] <- expt.cfp$biomass.total.true          # This is the true total biomass
      CFP[["no.fishing.tb"]] <- no.fishing    # True total biomass with no fishing
      
    }
    save(CFP,file=paste("All",s,"CFP",".RData",sep="_"))
  }
  
  if(HCR=="constF"){
    set.seed(123) # same seed
    for (sim in 1:nsims){
      equilib$B0 <- equilib.true$B0 * exp(B0.error[sim])
      rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
      expt.constF <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "constF", const.f.rate = const.f.rate, steepness = steepness,obs.type = obs.type,equilib=equilib,R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s)
      
      constF[["biomass.oneplus.true"]][sim,] <- expt.constF$biomass.oneplus.true
      constF[["total.catch"]][sim,] <- expt.constF$total.catch
      constF[["fishing"]][sim,] <- expt.constF$fishing
      constF[["intended.f"]][sim,] <- expt.constF$intended.f   
      constF[["rec"]][sim,] <- expt.constF$rec
      constF[["biomass.oneplus.obs"]][sim,] <- expt.constF$biomass.oneplus.obs
      constF[["biomass.total.true"]][sim,] <- expt.constF$biomass.total.true
      constF[["no.fishing.tb"]] <- no.fishing
    }
    save(constF,file=paste("All",s,"constF",".RData",sep="_"))
  }
  
  if(HCR=="C1"){
    set.seed(123) # same seed
    for (sim in 1:nsims){
      equilib$B0 <- equilib.true$B0 * exp(B0.error[sim])
      rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
      expt.c1 <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "C1",equilib = equilib,steepness=steepness,obs.type = obs.type,R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s)
      
      C1[["biomass.oneplus.true"]][sim,] <- expt.c1$biomass.oneplus.true
      C1[["total.catch"]][sim,] <- expt.c1$total.catch
      C1[["fishing"]][sim,] <- expt.c1$fishing
      C1[["intended.f"]][sim,] <- expt.c1$intended.f    
      C1[["rec"]][sim,] <- expt.c1$rec
      C1[["biomass.oneplus.obs"]][sim,] <- expt.c1$biomass.oneplus.obs
      C1[["biomass.total.true"]][sim,] <- expt.c1$biomass.total.true
      C1[["no.fishing.tb"]] <- no.fishing
    }
    save(C1,file=paste("All",s,"C1",".RData",sep="_")) 
  }
  
  if(HCR=="C2"){
    set.seed(123) # same seed
    for (sim in 1:nsims){
      equilib$B0 <- equilib.true$B0 * exp(B0.error[sim])
      rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
      expt.c2 <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "C2",equilib = equilib,steepness=steepness,obs.type = obs.type,R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s)
      
      C2[["biomass.oneplus.true"]][sim,] <- expt.c2$biomass.oneplus.true
      C2[["total.catch"]][sim,] <- expt.c2$total.catch
      C2[["fishing"]][sim,] <- expt.c2$fishing
      C2[["intended.f"]][sim,] <- expt.c2$intended.f    
      C2[["rec"]][sim,] <- expt.c2$rec
      #C2[["depl"]][sim,] <- expt.c2$depl
      C2[["biomass.oneplus.obs"]][sim,] <- expt.c2$biomass.oneplus.obs
      C2[["biomass.total.true"]][sim,] <- expt.c2$biomass.total.true
      C2[["no.fishing.tb"]] <- no.fishing
    }
    save(C2,file=paste("All",s,"C2",".RData",sep="_")) 
  }
  
  if(HCR=="C3"){
    set.seed(123) # same seed
    for (sim in 1:nsims){
      equilib$B0 <- equilib.true$B0 * exp(B0.error[sim])
      rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
      expt.c3 <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "C3",equilib = equilib,steepness=steepness,obs.type = obs.type,R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s)
      
      C3[["biomass.oneplus.true"]][sim,] <- expt.c3$biomass.oneplus.true
      C3[["total.catch"]][sim,] <- expt.c3$total.catch
      C3[["fishing"]][sim,] <- expt.c3$fishing
      C3[["intended.f"]][sim,] <- expt.c3$intended.f    
      C3[["rec"]][sim,] <- expt.c3$rec
      C3[["biomass.oneplus.obs"]][sim,] <- expt.c3$biomass.oneplus.obs
      C3[["biomass.total.true"]][sim,] <- expt.c3$biomass.total.true
      C3[["no.fishing.tb"]] <- no.fishing
    }
    save(C3,file=paste("All",s,"C3",".RData",sep="_")) 
  }
  
  
  if(HCR=="trend"){
    set.seed(123) # same seed
    for (sim in 1:nsims){
      equilib$B0 <- equilib.true$B0 * exp(B0.error[sim])
      rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
      expt.trend <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "trend",const.f.rate = 0.6,equilib = equilib,steepness=steepness,obs.type = obs.type,R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s)
      
      trend[["biomass.oneplus.true"]][sim,] <- expt.trend$biomass.oneplus.true
      trend[["total.catch"]][sim,] <- expt.trend$total.catch
      trend[["fishing"]][sim,] <- expt.trend$fishing
      trend[["intended.f"]][sim,] <- expt.trend$intended.f    # True total biomass with no fishing
      trend[["rec"]][sim,] <- expt.trend$rec
      trend[["biomass.oneplus.obs"]][sim,] <- expt.trend$biomass.oneplus.obs # Observed one-plus biomass
      trend[["biomass.total.true"]][sim,] <- expt.trend$biomass.total.true
      trend[["no.fishing.tb"]] <- no.fishing
    }
    save(trend,file=paste("All",s,"trend",".RData",sep="_")) 
  }
  
  if(HCR=="constF_HI"){
    const.f.rate = equilib$Fmsy
    set.seed(123) # same seed
    for (sim in 1:nsims){
      equilib$B0 <- equilib.true$B0 * exp(B0.error[sim])
      rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
      expt.constF_HI <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "constF", const.f.rate = const.f.rate, steepness = steepness,obs.type = obs.type,equilib=equilib,R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s)
      
      constF_HI[["biomass.oneplus.true"]][sim,] <- expt.constF_HI$biomass.oneplus.true
      constF_HI[["total.catch"]][sim,] <- expt.constF_HI$total.catch
      constF_HI[["fishing"]][sim,] <- expt.constF_HI$fishing
      constF_HI[["intended.f"]][sim,] <- expt.constF_HI$intended.f   
      constF_HI[["rec"]][sim,] <- expt.constF_HI$rec
      constF_HI[["biomass.oneplus.obs"]][sim,] <- expt.constF_HI$biomass.oneplus.obs
      constF_HI[["biomass.total.true"]][sim,] <- expt.constF_HI$biomass.total.true
      constF_HI[["no.fishing.tb"]] <- no.fishing
    }
    save(constF_HI,file=paste("All",s,"constF_HI",".RData",sep="_"))
  }
} # end of scenarios loop



#  Now load data (if needed) and summarize --------------------------------

hist()
