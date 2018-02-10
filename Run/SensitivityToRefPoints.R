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
nsims = 100
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
B0.accuracy = c("over","under","exact")

scenarios <- expand.grid(h,obs.error.type,HCR,recruit.sd,recruit.rho,M.type, B0.accuracy)
colnames(scenarios) <- c("h","obs.error.type","HCR","recruit.sd","recruit.rho","M.type","B0.accuracy")
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

# B0.over <- rtruncnorm(nsims,mean = 0,sd = B0.sd,a = 0.0001, b = Inf) #overestimates
# B0.under <- rtruncnorm(nsims,mean = 0,sd = B0.sd,a = -Inf, b = -0.001) # underestimates

# Think it makes more sense to just draw from a uniform distribution...
B0.over <- runif(nsims,0,1) # overestimates
B0.under <- runif(nsims,-1,0) # underestimates


# Test sensitivity to over-or under-estimating B0 ---------------------------------------------------

for(s in 1:nscenarios){  #
  steepness = scenarios$h[s]
  obs.type <- scenarios$obs.error.type[s]
  HCR <- scenarios$HCR[s]
  recruit.sd = scenarios$recruit.sd[s]
  recruit.rho = scenarios$recruit.rho[s]
  M.type = scenarios$M.type[s]
  B0.acc = scenarios$B0.accuracy[s] #****this is new
  equilib.true = getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=steepness) # NO recruitment devs used in the equilibrium calculations, so don't need to embed in the loop
  equilib <- equilib.true
  const.f.rate = 0.5*equilib$Fmsy
  no.fishing <- matrix(NA, nrow = nsims, ncol = years.test)
  set.seed(123) # Start each round of sims at same random seed
  time.var.m <- NA # Base case: M constant 
  
  for (sim in 1:nsims){
    rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
    F0.Type <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "constF", const.f.rate = 0, steepness = steepness,obs.type = obs.type,equilib=equilib,R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s)$biomass.total.true # Only need to do this 1x for each simulation (not repeat for each CR) because the seed is the same and there is no fishing.
    no.fishing[sim,] <- F0.Type # This is the time series 
  }
  
  
  if(HCR=="cfp"){    
    set.seed(123) # Start each round of sims at same random seed
    for (sim in 1:nsims){    # error in B0 will only affect control rules. Error is added here (instead of before with the unfished part) but want unfished dynamics determined with the true equilib #s 
      if(B0.acc=="over"){equilib$B0 <- equilib.true$B0 * exp(B0.over[sim])}else{
      if(B0.acc=="under"){equilib$B0 <- equilib.true$B0 * exp(B0.under[sim])}else{
        equilib$B0 <- equilib.true$B0
      }}
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
      if(B0.acc=="over"){equilib$B0 <- equilib.true$B0 * exp(B0.over[sim])}else{
        if(B0.acc=="under"){equilib$B0 <- equilib.true$B0 * exp(B0.under[sim])}else{
          equilib$B0 <- equilib.true$B0
        }}
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
      if(B0.acc=="over"){equilib$B0 <- equilib.true$B0 * exp(B0.over[sim])}else{
        if(B0.acc=="under"){equilib$B0 <- equilib.true$B0 * exp(B0.under[sim])}else{
          equilib$B0 <- equilib.true$B0
        }}
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
      if(B0.acc=="over"){equilib$B0 <- equilib.true$B0 * exp(B0.over[sim])}else{
        if(B0.acc=="under"){equilib$B0 <- equilib.true$B0 * exp(B0.under[sim])}else{
          equilib$B0 <- equilib.true$B0
        }}
      rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
      expt.c2 <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "C2",equilib = equilib,steepness=steepness,obs.type = obs.type,R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s)
      
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
      if(B0.acc=="over"){equilib$B0 <- equilib.true$B0 * exp(B0.over[sim])}else{
        if(B0.acc=="under"){equilib$B0 <- equilib.true$B0 * exp(B0.under[sim])}else{
          equilib$B0 <- equilib.true$B0
        }}
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
      if(B0.acc=="over"){equilib$B0 <- equilib.true$B0 * exp(B0.over[sim])}else{
        if(B0.acc=="under"){equilib$B0 <- equilib.true$B0 * exp(B0.under[sim])}else{
          equilib$B0 <- equilib.true$B0
        }}
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
      if(B0.acc=="over"){equilib$B0 <- equilib.true$B0 * exp(B0.over[sim])}else{
        if(B0.acc=="under"){equilib$B0 <- equilib.true$B0 * exp(B0.under[sim])}else{
          equilib$B0 <- equilib.true$B0
        }}
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
source("/Users/mcsiple/Dropbox/ChapterX-synthesis/Theme_Black.R")

path <- "~/Dropbox/Chapter4-HarvestControlRules/Results/Sensitivity"
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
o.i <- subset(scenarios,B0.accuracy=="over")$scenario
u.i <- subset(scenarios,B0.accuracy=="under")$scenario
e.i <- subset(scenarios,B0.accuracy=="exact")$scenario

# Summarize data (this takes a little while)
B0.overs <- lapply(results[o.i], FUN = summ.tab) %>%  # Summarize results (medians and quantiles)
  rbind.fill() %>% #turn list to dataframe
  mutate(scenario = rep(o.i, each=length(performance.measures))) %>%     # add column with scenario identifier
  left_join(scenarios, by = "scenario") #join with scenarios to get info about control rules etc. (tight!!!!)

B0.unders <- lapply(results[u.i], FUN = summ.tab) %>% 
  rbind.fill() %>% 
  mutate(scenario = rep(u.i, each=length(performance.measures))) %>%
  left_join(scenarios, by = "scenario")

B0.exacts <- lapply(results[e.i], FUN = summ.tab) %>% 
  rbind.fill()%>% 
  mutate(scenario = rep(e.i, each=length(performance.measures))) %>% 
  left_join(scenarios, by = "scenario")
  
head(B0.exacts)
head(B0.unders)
head(B0.overs)

plot(results[[15]]$biomass.oneplus.true[100,],type='l')
lines(results[[3]]$biomass.oneplus.true[100,],col='red')
lines(results[[9]]$biomass.oneplus.true[100,],col='blue')

plot(results[[15]]$total.catch[100,],type='l')
lines(results[[3]]$total.catch[100,],col='red')
lines(results[[9]]$total.catch[100,],col='blue')

plot(results[[15]]$intended.f[100,],type='l')
lines(results[[3]]$intended.f[100,],col='red')
lines(results[[9]]$intended.f[100,],col='blue')

plot(apply(results[[17]]$biomass.oneplus.true,MARGIN = 1,FUN = mean,na.rm=TRUE),
     apply(results[[17]]$total.catch,MARGIN = 1,FUN = mean,na.rm=TRUE),pch=19,col=alpha(colour = 'black',alpha = 0.5),
     xlim=c(2e3,9e3),ylim=c(0,2e3),xlab="Biomass",ylab="Catch")
points(apply(results[[5]]$biomass.oneplus.true,MARGIN = 1,FUN = mean,na.rm=TRUE),
       apply(results[[5]]$total.catch,MARGIN = 1,FUN = mean,na.rm=TRUE),pch=19,col=alpha(colour = 'red',alpha = 0.5))     
points(apply(results[[11]]$biomass.oneplus.true,MARGIN = 1,FUN = mean,na.rm=TRUE),
       apply(results[[11]]$total.catch,MARGIN = 1,FUN = mean,na.rm=TRUE),pch=19,col=alpha(colour = 'blue',alpha = 0.5))   

plot(B0.over,apply(results[[17]]$total.catch,MARGIN = 1,FUN = mean,na.rm=TRUE),pch=19,col=alpha(colour = 'black',alpha = 0.5))


toplot <- rbind.fill(subset(B0.exacts,HCR %in% c("C1","C2","C3")),
                     subset(B0.overs,HCR %in% c("C1","C2","C3")),
                     subset(B0.unders,HCR %in% c("C1","C2","C3")))

# Take inverse of the metrics that are "bad" --- i.e., SD(catch)
bad.pms <- c("SDcatch","n.5yrclose","n.10yrclose","nyrs0catch","SDbiomass","very.bad4preds","overallMaxCollapseLength","CollapseLength","Prob.Collapse","Collapse.Severity","CV.Catch","Bonafide.Collapse")
bi <- which(toplot$PM %in% bad.pms)



# Use old code to get raw performance metrics -----------------------------

for (s in 1:nscenarios){
  #**N** indicate metrics for which higher values mean worse performance (like SD(catch)) - these metrics are in scen.table as 1/x
  result.to.use <- results[[s]]
  raw.table[s,performance.measures[1]] <- median(rowMeans(result.to.use$total.catch[,calc.ind],na.rm = TRUE)) 
  #calculate mean B over years to use in the index - the final number is the median (across all simulations) mean B
  nonzero.catch <- result.to.use$total.catch[,calc.ind]
  nonzero.catch <- ifelse(nonzero.catch<0.1,NA,nonzero.catch)
  mnz.catches <- rowMeans(nonzero.catch,na.rm=TRUE)
  raw.table[s,performance.measures[2]] <- median(mnz.catches,na.rm = TRUE) #"LTnonzeromeancatch"
  raw.table[s,performance.measures[3]] <- median(apply(X = result.to.use$total.catch[,calc.ind],FUN = sd,MARGIN = 1)) 
  ######
  five.yr.closure.given1 <- n.multiyr.closures(result.to.use$total.catch[,calc.ind],threshold=0)$count5 / 
    n.multiyr.closures(result.to.use$total.catch[,calc.ind],threshold=0)$count1
  median.P5 <- median(five.yr.closure.given1,na.rm=TRUE) # This is now the # of 5-year closures given a 1-year closure
  if(all(n.multiyr.closures(result.to.use$total.catch[,calc.ind],threshold=0)$count5==0) &
     all(n.multiyr.closures(result.to.use$total.catch[,calc.ind],threshold=0)$count1==0)){
    median.P5 = 0
  }
  
  raw.table[s,performance.measures[4]] <- median.P5 # Mean number of 5-yr closures given a 1-yr closure
  
  ######
  ten.yr.closure.given5 <- n.multiyr.closures(result.to.use$total.catch[,calc.ind],threshold=0)$count10 / 
    n.multiyr.closures(result.to.use$total.catch[,calc.ind],threshold=0)$count5
  median.P10 <- median(ten.yr.closure.given5,na.rm=TRUE) 
  if(all(n.multiyr.closures(result.to.use$total.catch[,calc.ind],threshold=0)$count10==0) &
     all(n.multiyr.closures(result.to.use$total.catch[,calc.ind],threshold=0)$count5==0)){
    median.P10 = 0
  }
  
  raw.table[s,performance.measures[5]] <- median.P10 # Mean number of 10-yr closures given a 5-yr closure
  ######
  raw.table[s,performance.measures[6]] <- median(apply(X = result.to.use$total.catch[,calc.ind],FUN = nzeroes,MARGIN = 1))
  
  raw.table[s,performance.measures[7]] <- median(rowMeans(result.to.use$biomass.total.true[,calc.ind])) # "meanbiomass"
  
  raw.table[s,performance.measures[8]] <- median(apply(X = result.to.use$biomass.total.true[,calc.ind],
                                                       FUN = good4pred,MARGIN = 1, 
                                                       F0.x = result.to.use$no.fishing.tb[,calc.ind])) # "good4preds" - this is calculated from TOTAL biomass (including age 0)
  
  raw.table[s,performance.measures[9]] <- median(apply(X = result.to.use$biomass.total.true[,calc.ind],FUN = sd,MARGIN = 1))  #Actual raw SD of Biomass
  
  yrs.bad <- apply(X = result.to.use$biomass.total.true[,calc.ind],FUN = bad4pred,MARGIN = 1, F0.x = result.to.use$no.fishing.tb[,calc.ind] ) # Number of years that are very bad for preds - length of vector is nsims 
  
  raw.table[s,performance.measures[10]] <- median(yrs.bad)
  
  sw <- subset(all.summaries,scenario==s)
  for(pm in 11:19){
    select.ind <- which(sw$PM == performance.measures[pm]) 
    raw.table[s,performance.measures[pm]] <- sw[select.ind,'med'] 
  }
  # 1. Depletion
  # 2. max collapse length
  # 3. **N**max bonanza length
  # 4. Mean bonanza length
  # 5. Mean collapse length
  # 6. # probability of collapse 
  # 7. Collapse severity
  # 8. CV.Catch
  # 9. Bonafide.collapse
} 

