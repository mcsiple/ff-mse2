# Sensitivity to differences in life history params
# We want to know which parameters exactly are having the strongest influence on tradeoffs and/or differences in control rule performance. 
# This code:
# - Takes one life history type - anchovy
# - changes one lh param at a time: 
  # - natural mortality 
  # - max age
  # - selectivity-maturity combo


# Libraries
library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)

fftype = "Anchovy"  

# Set directories
basedir <- "~/Dropbox/Chapter4-HarvestControlRules/Code/ff-mse2"
lhparam <- "mortality"
resultsdir <- paste("~/Dropbox/Chapter4-HarvestControlRules/Results/Sensitivity/",lhparam,sep="")
subDir <- fftype

#Set up other simulation params
years.test = 250
nsims = 1000
R0.sens = NA # no dynamic R0
tim.params = list(sigma0 = 0.2,tau0 = 0.1)
tau1 = (1/tim.params$tau0^2 + 1/tim.params$sigma0^2)^(-0.5)
sig.s = 0.3 


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
  # lh params go below so they can reset in the for loop
  # Rec dev params
  recruit.sd <- 0.6
  recruit.rho <- 0.5
}


# Create errors for AC scenario
set.seed(123)
curly.phi.mat <- matrix(NA, nrow = nsims,ncol = years.test)
for(sim in 1:nsims){
  curly.phi.mat[sim,] <- rnorm(years.test,0,sig.s)
}

# Create list of matrices with error for the delayed detection scenario. Important bc of random seed issues.
set.seed(123)
tim.rands.list.short <- list() #for all ensuring random values
n.ages = 7
tim.inits.vec.short <- rnorm(n.ages,0,tim.params$sigma0)  # just for initial values
for(sim in 1:nsims){
  tim.mat <- matrix(NA,nrow=n.ages,ncol=years.test)
  for(i in 1:years.test){
    tim.mat[,i] <- rnorm(n.ages,0,tau1)
  }
  tim.rands.list.short[[sim]] <- tim.mat
}


# another set of rands, for when there are older ages ---------------------
set.seed(123)
tim.rands.list.long <- list() #for all ensuring random values
n.ages = 16
tim.inits.vec.long <- rnorm(n.ages,0,tim.params$sigma0)  # just for initial values
for(sim in 1:nsims){
  tim.mat <- matrix(NA,nrow=n.ages,ncol=years.test)
  for(i in 1:years.test){
    tim.mat[,i] <- rnorm(n.ages,0,tau1)
  }
  tim.rands.list.long[[sim]] <- tim.mat
}
############# end random number generation


# Load harvest rules
source(file.path(basedir,"Control Rules/smith_oceana.R"))
source(file.path(basedir,"Control Rules/cfp.R"))
source(file.path(basedir,"Control Rules/hockey-stick.R"))
source(file.path(basedir,"Control Rules/trend-based-rule.R"))


# Scenarios - for sensitivity, limiting the number of scenarios.
setwd(resultsdir)
h = c(0.6) 
obs.error.type = c("AC") #,"noerror","Tim"
HCR = c("constF","C1","C2","C3","constF_HI") 
M.type = c("constant") # took out "regimeshift" and "time-varying" to save time but can be added back in for sensitivity analyses
LH <- c("base","lowM","select.old","olderfish","sard.devs")

scenarios <- expand.grid(h,obs.error.type,HCR,recruit.sd,recruit.rho,M.type,LH)
colnames(scenarios) <- c("h","obs.error.type","HCR","recruit.sd","recruit.rho","M.type","LH")
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



# Test sensitivity to over-or under-estimating B0 ---------------------------------------------------

for(s in 1:nscenarios){  #
  source(file.path(basedir,"Ctl/Anchovy_LHControl.R"))
  source(file.path(basedir,"Ctl/Anchovy_FisheryControl.R"))
  
  steepness = scenarios$h[s]
  obs.type <- scenarios$obs.error.type[s]
  HCR <- scenarios$HCR[s]
  recruit.sd = scenarios$recruit.sd[s]
  recruit.rho = scenarios$recruit.rho[s]
  M.type = scenarios$M.type[s]
  LH = scenarios$LH[s]
  tim.rands.list = tim.rands.list.short
  tim.inits.vec = tim.inits.vec.short
  # LH changes --------------------------------------------------------------
  # First: base case.
  # Second: low natural mortality
  if(LH=="lowM"){
    lh.test$M = 0.5 
  }
  # Third: Selectivity shifted to the right (fish mature before they are taken by the fishery)
  if(LH=="select.old"){
  lh.test$selectivity[3:7,2] <- lh.test$selectivity[1:5,2]
  lh.test$selectivity[1:2,2] <- 0
  }
              #plot(lh.test$ages,lh.test$maturity)
              #lines(lh.test$ages,lh.test$selectivity[,2])
  # Fourth: More ages... 
  if(LH=="olderfish"){
      lh.test$selectivity <- rbind(lh.test$selectivity, matrix(c(7:15,rep(1,times=9)),nrow=9,ncol=2))
      lh.test$ages <- 0:15
      lh.test$l.at.age <- c(lh.test$l.at.age,rep(20, times=9)) # made this up because length isn't used in this model
      lh.test$w.at.age <- c(lh.test$w.at.age, c(8.2,8.3,8.4,rep(8.5, times=6)))
      lh.test$maturity <- c(lh.test$maturity,rep(1,times=9))
      tim.rands.list = tim.rands.list.long
      tim.inits.vec = tim.inits.vec.long
      
      tot <- sum(c(6.221,4.170,2.795,1.873,1.256,0.842,0.564,rep(0.6,times=9)))
      init.prop <- c(6.221,4.170,2.795,1.873,1.256,0.842,0.564,rep(0.6,times=9))/tot 
      init.B <- init.prop*20000 
      init.test <- init.B / lh.test$w.at.age 
  }
  # Fifth: Everything the same except rec devs are different
  if(LH=="sard.devs"){
   recruit.sd <- 0.6
   recruit.rho <- 0.9
  }
  
  
  
  equilib = getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=steepness) # NO recruitment devs used in the equilibrium calculations, so don't need to embed in the loop
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



