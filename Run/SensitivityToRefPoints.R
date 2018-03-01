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
metric = "Fmsy" # Choose which rp to test sensitivity
resultsdir <- paste("~/Dropbox/Chapter4-HarvestControlRules/Results/Sensitivity/",metric,sep="")
subDir <- fftype

#Set up other simulation params
years.test = 250
nsims = 500
R0.sens = NA # NO DYNAMIC R0 anymore-- ignore
tim.params = list(sigma0 = 0.2,tau0 = 0.1)
tau1 = (1/tim.params$tau0^2 + 1/tim.params$sigma0^2)^(-0.5)
sig.s = 0.3 

# Load packages
library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library(truncnorm) #For truncated normal (adding bias)

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

# Create list of matrices with error for the delayed detection scenario. This is important because the random seed will offset otherwise - this is new
        set.seed(123)
        tim.rands.list <- list() #for all ensuring random values
        tim.inits.vec <- rnorm(n.ages,0,tim.params$sigma0)  # just for initial values
        n.ages = length(lh.test$ages)
        for(sim in 1:nsims){
          tim.mat <- matrix(NA,nrow=n.ages,ncol=years.test)
          for(i in 1:years.test){
            tim.mat[,i] <- rnorm(n.ages,0,tau1)
          }
          tim.rands.list[[sim]] <- tim.mat
        }
        
        # Create errors for AC scenario
        set.seed(123)
        curly.phi.mat <- matrix(NA, nrow = nsims,ncol = years.test)
        for(sim in 1:nsims){
          curly.phi.mat[sim,] <- rnorm(years.test,0,sig.s)
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

# FIRST: make sure you chose which metric you want sensitivity for!
accuracy = c("over","under","exact")

scenarios <- expand.grid(h,obs.error.type,HCR,recruit.sd,recruit.rho,M.type, accuracy)
colnames(scenarios) <- c("h","obs.error.type","HCR","recruit.sd","recruit.rho","M.type","accuracy")
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
      # steepness = scenarios$h[1]
      # obs.type <- scenarios$obs.error.type[1]
      # HCR <- scenarios$HCR[1]
      # equilib = getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=steepness)
      # rec.dev.test <- generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd)
      # test.constF <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test,rec.ram = NA, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "constF", const.f.rate = 0, steepness = steepness,obs.type = obs.type,equilib=equilib,R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s, tim.rand.inits = tim.inits.vec, tim.rands = tim.rands.list[[sim]],curly.phi.vec = curly.phi.mat[sim,]) # Need to fix the rest of these
      # nofish <- melt(test.constF[c("biomass.oneplus.obs","biomass.total.obs","total.catch","fishing","biomass.total.true","rec","biomass.oneplus.true")])
      # nofish$year <- rep(1:years.test,times=length(unique(nofish$L1)))
      # ggplot(nofish,aes(x=year,y=value)) + geom_line() + facet_wrap(~L1,scales = "free_y") #+ xlim(c(150,250))
      # 

                  # Get appropriate SD for variation in F -----------------------------------
                  #nsims = 100
                  # sd.sim = vector(length = nsims)
                  # set.seed(123)
                  # for(sim in 1:nsims){
                  #   rec.dev.test <- generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd)
                  #   test.F.sd <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test,rec.ram = NA, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "constF", const.f.rate = 0.7, steepness = steepness,obs.type = "AC",equilib=equilib,R0.traj = R0.sens, tim.params = tim.params,time.var.m = NA,sig.s = sig.s)
                  #   sd.sim[sim] <- sd(test.F.sd$fishing[150:250]) #last 100 years of time series
                  # }
                  # median(sd.sim)
                  
                  # With sig.s = 0.3 and F = 0.7, this is the sd(fishing mortality) that you would expect with AC error
                  # F = 1.1    0.5602482
                  # F = 0.7    0.2956493
                  # F = 0.3   0.1047076
                  # F = 0.1   0.03187569
                  # So sd of F depends on the constant F rate that you set... I think there's an adjustment but I'm not sure what it is... might as well set sd high and see what happens.
                  #Going to use the high value of sd(F) just to be precautionary. sd = 0.56

Fmsy.sd <- 0.56
B0.sd <- 0.3

# Make vectors of modifier to B0, to test over- and under-estimates of B0 against accurate estimates
B0.over <- rtruncnorm(nsims,mean = 0,sd = B0.sd,a = 0.0001, b = Inf) # overestimates
B0.under <- rtruncnorm(nsims,mean = 0,sd = B0.sd,a = -Inf, b = -0.001) # underestimates

Fmsy.over <- rtruncnorm(nsims,mean = 0,sd = Fmsy.sd,a = 0.0001, b = Inf) 
Fmsy.under <- rtruncnorm(nsims,mean = 0,sd = Fmsy.sd,a = -Inf, b = -0.001) 
# Does it makes more sense to just draw from a uniform distribution?



# Test sensitivity to over-or under-estimating B0 ---------------------------------------------------

for(s in 1:nscenarios){  #
  steepness = scenarios$h[s]
  obs.type <- scenarios$obs.error.type[s]
  HCR <- scenarios$HCR[s]
  recruit.sd = scenarios$recruit.sd[s]
  recruit.rho = scenarios$recruit.rho[s]
  M.type = scenarios$M.type[s]
  acc = scenarios$accuracy[s] # ****this is new****
  equilib.true = getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=steepness) # NO recruitment devs used in the equilibrium calculations, so don't need to embed in the loop
  equilib <- equilib.true
  const.f.rate = 0.5*equilib$Fmsy
  no.fishing <- matrix(NA, nrow = nsims, ncol = years.test)
  set.seed(123) # Start each round of sims at same random seed
  time.var.m <- NA # Base case: M constant 
  
  for (sim in 1:nsims){
    rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
    F0.Type <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "constF", const.f.rate = 0, steepness = steepness,obs.type = obs.type,equilib=equilib,R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s, tim.rand.inits = tim.inits.vec, tim.rands = tim.rands.list[[sim]],curly.phi.vec = curly.phi.mat[sim,])$biomass.total.true # Only need to do this 1x for each simulation (not repeat for each CR) because the seed is the same and there is no fishing.
    no.fishing[sim,] <- F0.Type # This is the time series 
  }
  
  
  if(HCR=="cfp"){    
    set.seed(123) # Start each round of sims at same random seed
    for (sim in 1:nsims){    # error in B0 will only affect control rules. Error is added here (instead of before with the unfished part) but want unfished dynamics determined with the true equilib #s 
      if(metric=="B0"){
      if(acc=="over"){equilib$B0 <- equilib.true$B0 * exp(B0.over[sim])}else{
      if(acc=="under"){equilib$B0 <- equilib.true$B0 * exp(B0.under[sim])}else{
        equilib$B0 <- equilib.true$B0
      }}
      }
      if(metric=="Fmsy"){
        if(acc=="over"){equilib$Fmsy <- equilib.true$Fmsy * exp(Fmsy.over[sim])}else{
          if(acc=="under"){equilib$Fmsy <- equilib.true$Fmsy * exp(Fmsy.under[sim])}else{
            equilib$Fmsy <- equilib.true$Fmsy
      }}
      }
      rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
      expt.cfp <- calc.trajectory(lh = lh.test,obs.cv = NA, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "cfp",equilib = equilib,steepness=steepness,obs.type = obs.type,R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s, tim.rand.inits = tim.inits.vec, tim.rands = tim.rands.list[[sim]],curly.phi.vec = curly.phi.mat[sim,])
      
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
      if(metric=="B0"){
        if(acc=="over"){equilib$B0 <- equilib.true$B0 * exp(B0.over[sim])}else{
          if(acc=="under"){equilib$B0 <- equilib.true$B0 * exp(B0.under[sim])}else{
            equilib$B0 <- equilib.true$B0
          }}
      }
      if(metric=="Fmsy"){
        if(acc=="over"){equilib$Fmsy <- equilib.true$Fmsy * exp(Fmsy.over[sim])}else{
          if(acc=="under"){equilib$Fmsy <- equilib.true$Fmsy * exp(Fmsy.under[sim])}else{
            equilib$Fmsy <- equilib.true$Fmsy
          }}
      }
      rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
      expt.constF <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "constF", const.f.rate = const.f.rate, steepness = steepness,obs.type = obs.type,equilib=equilib,R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s, tim.rand.inits = tim.inits.vec, tim.rands = tim.rands.list[[sim]],curly.phi.vec = curly.phi.mat[sim,])
      
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
      if(metric=="B0"){
        if(acc=="over"){equilib$B0 <- equilib.true$B0 * exp(B0.over[sim])}else{
          if(acc=="under"){equilib$B0 <- equilib.true$B0 * exp(B0.under[sim])}else{
            equilib$B0 <- equilib.true$B0
          }}
      }
      # if(metric=="Fmsy"){
      #   if(acc=="over"){equilib$Fmsy <- equilib.true$Fmsy * exp(Fmsy.over[sim])}else{
      #     if(acc=="under"){equilib$Fmsy <- equilib.true$Fmsy * exp(Fmsy.under[sim])}else{
      #       equilib$Fmsy <- equilib.true$Fmsy
      #     }}
      # }
      rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
      expt.c1 <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "C1",equilib = equilib,steepness=steepness,obs.type = obs.type,R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s, tim.rand.inits = tim.inits.vec, tim.rands = tim.rands.list[[sim]],curly.phi.vec = curly.phi.mat[sim,])
      
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
    #*******
    # c1.over <- expt.c1
    # c1.under <- expt.c1
    # c1.under.2 <- expt.c1
    # plot(c1.over$biomass.total.true,type='l')
    # lines(c1.under$biomass.total.true,col='red')
    # lines(c1.under.2$biomass.total.true,col='blue')
  }
  
  if(HCR=="C2"){
    set.seed(123) # same seed
    for (sim in 1:nsims){
      if(metric=="B0"){
        if(acc=="over"){equilib$B0 <- equilib.true$B0 * exp(B0.over[sim])}else{
          if(acc=="under"){equilib$B0 <- equilib.true$B0 * exp(B0.under[sim])}else{
            equilib$B0 <- equilib.true$B0
          }}
      }
      # if(metric=="Fmsy"){
      #   if(acc=="over"){equilib$Fmsy <- equilib.true$Fmsy * exp(Fmsy.over[sim])}else{
      #     if(acc=="under"){equilib$Fmsy <- equilib.true$Fmsy * exp(Fmsy.under[sim])}else{
      #       equilib$Fmsy <- equilib.true$Fmsy
      #     }}
      # }
      rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
      expt.c2 <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "C2",equilib = equilib,steepness=steepness,obs.type = obs.type,R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s, tim.rand.inits = tim.inits.vec, tim.rands = tim.rands.list[[sim]],curly.phi.vec = curly.phi.mat[sim,])
      
      C2[["biomass.oneplus.true"]][sim,] <- expt.c2$biomass.oneplus.true
      C2[["total.catch"]][sim,] <- expt.c2$total.catch
      C2[["fishing"]][sim,] <- expt.c2$fishing
      C2[["intended.f"]][sim,] <- expt.c2$intended.f    
      C2[["rec"]][sim,] <- expt.c2$rec
      C2[["biomass.oneplus.obs"]][sim,] <- expt.c2$biomass.oneplus.obs
      C2[["biomass.total.true"]][sim,] <- expt.c2$biomass.total.true
      C2[["no.fishing.tb"]] <- no.fishing
    }
    save(C2,file=paste("All",s,"C2",".RData",sep="_")) 
  }
  
  if(HCR=="C3"){
    set.seed(123) # same seed
    for (sim in 1:nsims){
      if(metric=="B0"){
        if(acc=="over"){equilib$B0 <- equilib.true$B0 * exp(B0.over[sim])}else{
          if(acc=="under"){equilib$B0 <- equilib.true$B0 * exp(B0.under[sim])}else{
            equilib$B0 <- equilib.true$B0
          }}
      }
      if(metric=="Fmsy"){
        if(acc=="over"){equilib$Fmsy <- equilib.true$Fmsy * exp(Fmsy.over[sim])}else{
          if(acc=="under"){equilib$Fmsy <- equilib.true$Fmsy * exp(Fmsy.under[sim])}else{
            equilib$Fmsy <- equilib.true$Fmsy
          }}
      }
      rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
      expt.c3 <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "C3",equilib = equilib,steepness=steepness,obs.type = obs.type,R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s, tim.rand.inits = tim.inits.vec, tim.rands = tim.rands.list[[sim]],curly.phi.vec = curly.phi.mat[sim,])
      
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
      if(metric=="B0"){
        if(acc=="over"){equilib$B0 <- equilib.true$B0 * exp(B0.over[sim])}else{
          if(acc=="under"){equilib$B0 <- equilib.true$B0 * exp(B0.under[sim])}else{
            equilib$B0 <- equilib.true$B0
          }}
      }
      # if(metric=="Fmsy"){
      #   if(acc=="over"){equilib$Fmsy <- equilib.true$Fmsy * exp(Fmsy.over[sim])}else{
      #     if(acc=="under"){equilib$Fmsy <- equilib.true$B0 * exp(Fmsy.under[sim])}else{
      #       equilib$Fmsy <- equilib.true$Fmsy
      #     }}
      # }
      rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
      expt.trend <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "trend",const.f.rate = 0.6,equilib = equilib,steepness=steepness,obs.type = obs.type,R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s, tim.rand.inits = tim.inits.vec, tim.rands = tim.rands.list[[sim]],curly.phi.vec = curly.phi.mat[sim,])
      
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
      if(metric=="B0"){
        if(acc=="over"){equilib$B0 <- equilib.true$B0 * exp(B0.over[sim])}else{
          if(acc=="under"){equilib$B0 <- equilib.true$B0 * exp(B0.under[sim])}else{
            equilib$B0 <- equilib.true$B0
          }}
      }
      if(metric=="Fmsy"){
        if(acc=="over"){equilib$Fmsy <- equilib.true$Fmsy * exp(Fmsy.over[sim])}else{
          if(acc=="under"){equilib$Fmsy <- equilib.true$Fmsy * exp(Fmsy.under[sim])}else{
            equilib$Fmsy <- equilib.true$Fmsy
          }}
      }
      rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
      expt.constF_HI <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "constF", const.f.rate = const.f.rate, steepness = steepness,obs.type = obs.type,equilib=equilib,R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s, tim.rand.inits = tim.inits.vec, tim.rands = tim.rands.list[[sim]],curly.phi.vec = curly.phi.mat[sim,])
      
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

library(gtools)
library(RcppRoll)
library(plyr); library(dplyr)
library(scales)
library(gridExtra)
library(extrafont)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
source("/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Code/ff-mse2/Plots/Megsieggradar.R")
source("/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Code/ff-mse2/Plots/SummaryFxns.R")
#source("/Users/mcsiple/Dropbox/ChapterX-synthesis/Theme_Black.R")

path <- paste("~/Dropbox/Chapter4-HarvestControlRules/Results/Sensitivity/",metric,sep="")
files <- list.files(path=path)
rm <- grep(files,pattern = ".txt") # Don't load the text table summary
files <- files[-rm]
files <- files[grep(files, pattern = ".RData")] #only load rdata files
files <- mixedsort(files) # IMPORTANT: need this line to order in same order as scenario table!
results <- sapply(files, function(x) mget(load(x)),simplify=TRUE) # This is a giant list of all results
nscenarios <- length(results)
scenarios <- raw.table <- read.table("Scenario_Table.txt") 
nsims <- nrow(results[[1]]$biomass.oneplus.true) # count nrows to know how many sims there are
years.test <- ncol(results[[1]]$biomass.oneplus.true) # count cols to know how many years there are
nyrs.to.use <- 100 # How many years you want to use to calculate all your metrics - There are no big differences btwn 50 and 100 yrs
calc.ind <- tail(1:years.test, nyrs.to.use) # Which years to calculate median depletion over (length = nyrs.to.use)

# Add performance measure columns to table
performance.measures <- c("LTmeancatch","LTnonzeromeancatch","SDcatch","n.5yrclose","n.10yrclose","nyrs0catch","meanbiomass","good4preds","SDbiomass","very.bad4preds","meanDepl","overallMaxCollapseLength","overallMaxBonanzaLength","BonanzaLength","CollapseLength","Prob.Collapse","Collapse.Severity","CV.Catch","Bonafide.Collapse")
pm.type <- c(rep("Fishery",times=6),rep("Ecosystem",times=11),"Fishery","Ecosystem") # for distinguishing types of PMs (mostly for plotting...)
raw.table[,performance.measures] <- NA


# ------------------------------------------------------------------------
# Summarize everything in one giant "outputs" table
# ------------------------------------------------------------------------

# Fxns for summarizing and plotting ---------------------------------------

all.summaries <- NA

#Indices for which results were over- and under-estimates of B0
o.i <- subset(scenarios,accuracy=="over")$scenario
u.i <- subset(scenarios,accuracy=="under")$scenario
e.i <- subset(scenarios,accuracy=="exact")$scenario

# Summarize data (this takes a little while)
overs <- lapply(results[o.i], FUN = summ.tab) %>%  # Summarize results (medians and quantiles)
  rbind.fill() %>% # turn list to dataframe
  mutate(scenario = rep(o.i, each=length(performance.measures))) %>%     # add column with scenario identifier
  left_join(scenarios, by = "scenario") #join with scenarios to get info about control rules etc. (tight!!!!)

unders <- lapply(results[u.i], FUN = summ.tab) %>% 
  rbind.fill() %>% 
  mutate(scenario = rep(u.i, each=length(performance.measures))) %>%
  left_join(scenarios, by = "scenario")

exacts <- lapply(results[e.i], FUN = summ.tab) %>% 
  rbind.fill()%>% 
  mutate(scenario = rep(e.i, each=length(performance.measures))) %>% 
  left_join(scenarios, by = "scenario")
  
head(exacts)
head(unders)
head(overs)


# Compare the performance metrics between exact and over- or under --------
all.results <- rbind.fill(exacts,unders,overs)
# Subset to the control rules that use metric m:
if(metric == "B0"){plot.results <- subset(all.results, HCR %in% c("C1","C2","C3"))
                    col.to.use <- hcr.colors[1:3]}
if(metric == "Fmsy"){plot.results <- subset(all.results, HCR %in% c("constF","constF_HI","cfp","C3"))
                    col.to.use <- hcr.colors[c(3,6,4,5)]}

# Recode performance measures for plotting and take out PMs that aren't used in the paper
remove.these <- c("n.10yrclose","SDbiomass","meanDepl","LTnonzeromeancatch","good4preds","very.bad4preds","CV.Catch","overallMaxCollapseLength","overallMaxBonanzaLength","Bonafide.Collapse")
plot.results <- plot.results %>% filter(!PM %in% remove.these) %>%
  mutate(HCR = recode(HCR, 'cfp' = 'Stability-favoring',
             'constF' = 'Constant F',
             'C1' = 'Basichockey',
             'C2' = 'Low Blim',
             'C3' = 'High Fmax',
             'trend' = "Trend-based",
             'constF_HI' = "Constant F - High"),
             PM = recode(PM, 'LTmeancatch' = "Mean catch", # These are slightly different than in the Summarize_MSE file, because these aren't re-scaled for "negative metrics" or scaled to max
                             "meanbiomass" = "Mean biomass",
                             "BonanzaLength" = "Bonanza length",
                             'SDcatch' = "Catch variation (SD)",
                             'nyrs0catch' = "Years with 0 catch",
                             'n.5yrclose' = "P(5 yr closure|closure)",
                             "CollapseLength" = "Collapse length",
                             "Prob.Collapse" = "P(collapse)",
                             "Collapse.Severity" = "Collapse severity"))

plots <- ggplot(plot.results,aes(x=HCR,y=med,fill=HCR,shape=accuracy,colour=HCR)) + 
  geom_point(size=3,position = position_dodge(1.1)) + 
  geom_pointrange(aes(ymin=loCI,ymax=hiCI),position = position_dodge(1.1)) +
  facet_wrap(~PM,scales="free_y") + 
  ylab("Value") +
  scale_colour_manual(values = col.to.use) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(paste(metric,"Sensitivity.pdf",sep="_"),width = 10,height = 6,useDingbats = FALSE)
plots
dev.off()




# Individual plots/notes/messy stuff --------------------------------------------------------

par(mfrow=c(3,1))
sim = 150
# These plots are for C1
plot(results[[15]]$biomass.oneplus.true[sim,],type='l',ylab="True 1+ Biomass") #exact
lines(results[[3]]$biomass.oneplus.true[sim,],col='red') #over
lines(results[[9]]$biomass.oneplus.true[sim,],col='blue') #under

plot(results[[15]]$biomass.oneplus.obs[sim,],type='l',ylab="Observed 1+ Biomass") #exact
lines(results[[3]]$biomass.oneplus.obs[sim,],col='red') #over
lines(results[[9]]$biomass.oneplus.obs[sim,],col='blue') #under

plot(results[[15]]$total.catch[sim,],type='l',ylab="Total Catch")
lines(results[[3]]$total.catch[sim,],col='red')
lines(results[[9]]$total.catch[sim,],col='blue')

plot(results[[15]]$intended.f[sim,],type='l',ylab="F_target")
lines(results[[3]]$intended.f[sim,],col='red')
lines(results[[9]]$intended.f[sim,],col='blue')

plot(apply(results[[17]]$biomass.oneplus.true,MARGIN = 1,FUN = mean,na.rm=TRUE),
     apply(results[[17]]$total.catch,MARGIN = 1,FUN = mean,na.rm=TRUE),pch=19,col=alpha(colour = 'black',alpha = 0.5),
     xlim=c(2e3,9e3),ylim=c(0,2e3),xlab="Biomass",ylab="Catch")
points(apply(results[[5]]$biomass.oneplus.true,MARGIN = 1,FUN = mean,na.rm=TRUE),
       apply(results[[5]]$total.catch,MARGIN = 1,FUN = mean,na.rm=TRUE),pch=19,col=alpha(colour = 'red',alpha = 0.5))     
points(apply(results[[11]]$biomass.oneplus.true,MARGIN = 1,FUN = mean,na.rm=TRUE),
       apply(results[[11]]$total.catch,MARGIN = 1,FUN = mean,na.rm=TRUE),pch=19,col=alpha(colour = 'blue',alpha = 0.5))   

toplot <- rbind.fill(subset(B0.exacts,HCR %in% c("C1","C2","C3")),
                     subset(B0.overs,HCR %in% c("C1","C2","C3")),
                     subset(B0.unders,HCR %in% c("C1","C2","C3")))

# Take inverse of the metrics that are "bad" --- i.e., SD(catch)
bad.pms <- c("SDcatch","n.5yrclose","n.10yrclose","nyrs0catch","SDbiomass","very.bad4preds","overallMaxCollapseLength","CollapseLength","Prob.Collapse","Collapse.Severity","CV.Catch","Bonafide.Collapse")
bi <- which(toplot$PM %in% bad.pms)
toplot[bi,c("loCI","med","hiCI")] <- apply(toplot[bi,c("loCI","med","hiCI")], MARGIN = 2,FUN = function(x) ifelse(x==0,1,1/x))
# Take out metrics not being used in pairs plots:
remove.these <- c("n.10yrclose","SDbiomass","meanDepl","LTnonzeromeancatch","good4preds","very.bad4preds","CV.Catch","overallMaxCollapseLength","overallMaxBonanzaLength","Bonafide.Collapse")

# Remove quantiles because now they're a mess and I don't need em anyway
toplot <- toplot[,-c(2,4)]

# Scale all performance measures to max (1 = best performance)
toplot2 <- toplot %>% group_by(PM,B0.accuracy) %>% 
  mutate(scaled_med = med/max(med,na.rm=T))%>% 
  select(c(PM,HCR,B0.accuracy,scaled_med)) %>% 
  filter(!PM %in% remove.these) %>%
  as.data.frame()

if(length(which(toplot2$scaled_med==1))<5){print("Stop! Problems with summarizing in plyr!")}

# These are to check that there are maxes for everything...
subset(toplot2,PM=="meanbiomass")
subset(toplot2,PM=="SDcatch")

combo <- dcast(toplot2,HCR+B0.accuracy~PM, value.var = "scaled_med") # Can use ggplot if you want

# An easier way to do a pairs plot:
palette <- brewer.pal(6,"Spectral")
hcr.colors <- palette[c(6,5,4,3,1,2)]
combo$color <- rep(hcr.colors[1:3],each=length(unique(combo$B0.accuracy))) # colours for C1, C2, C3: hcr.colors[1:3]
pairs(combo[4:ncol(combo)-1],pch=rep(c(21,24,25),times=3),bg=combo$color)






# DEPRACATED: Test fitting a 1/x function to the exact estimates ---------------------

test <- subset(combo,B0.accuracy=="exact") %>% arrange(nyrs0catch)
x = test$nyrs0catch
y = test$LTmeancatch
fit <- lm(data = test,y~I(1/x))
plot(x,y,xlim=c(0,1),ylim=c(0,1))
newdata <- data.frame(x = seq(0.05,1,by=0.1))
pt <- predict(fit,newdata=newdata)
lines(newdata$x,pt,col="blue")


# Try a new thing: quantify tradeoffs with 1/x fitting --------------------

toplot.raw <- toplot %>% group_by(PM,B0.accuracy) %>% 
  #mutate(scaled_med = med/max(med,na.rm=T))%>% 
  select(c(PM,HCR,B0.accuracy,med)) %>% 
  filter(!PM %in% remove.these) %>%
  as.data.frame()

test2 <- subset(toplot.raw,B0.accuracy=="under")
x = subset(test2,PM=="BonanzaLength")$med
y = subset(test2,PM=="LTmeancatch")$med
fit <- lm(data = test2,y~I(1/x))
plot(x,y,xlim=c(0,max(x)*2),ylim=c(0,max(y)*2))
newdata <- data.frame(x = seq(0.05,max(x)*4,by=1))
pt <- predict(fit,newdata=newdata)
lines(newdata$x,pt,col="blue")
fit$coefficients


#  using raw outputs instead of tplot - think I will keep this one --------
par(mfrow=c(1,1))
raw <- B0.overs %>% group_by(PM,B0.accuracy) %>%
  select(c(PM,HCR,B0.accuracy,med)) %>% 
  filter(!PM %in% remove.these) %>%
  as.data.frame()
x = subset(raw,PM=="BonanzaLength")$med
y = subset(raw,PM=="meanbiomass")$med
fit <- lm(y~I(1/x))
plot(x,y,xlim=c(0,max(x)*2),ylim=c(0,max(y)*2))
newdata <- data.frame(x = seq(0.05,max(x)*4,by=1))
pt <- predict(fit,newdata=newdata)
lines(newdata$x,pt,col="blue")
fit$coefficients


# Maybe later transfer for summary code: quantify tradeoffs, redo tileplot --------
# using B0.exacts as a test case for now
totile <- B0.exacts %>% group_by(PM,B0.accuracy) %>%
  select(c(PM,HCR,B0.accuracy,med)) %>% 
  filter(!PM %in% remove.these) %>%
  as.data.frame()
totile <- dcast(totile,HCR~PM, value.var = "med")

# Need to get inverse of "bad.pms"
bi <- which(colnames(totile) %in% bad.pms)
totile[,bi] <- 1/totile[,bi]


to.mat <- matrix(NA,nrow=ncol(totile),ncol=ncol(totile))
for(i in 2:ncol(totile)){
  x <- totile[,i]
  for(j in 2:ncol(totile)){
  y <- totile[,j]
  if(any(is.na(c(x,y)))){next}
  fit <- lm(y~I(1/x))
  to.mat[i,j] <- fit$coefficients[2]
}}

