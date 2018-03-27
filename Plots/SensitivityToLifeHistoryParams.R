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


###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
# Read them all in; plot and stuff ----------------------------------------


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
Type = "Anchovy" #FF type to summarize

# Set path to wherever the simulation results are:
path <- paste("/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Results/Sensitivity/mortality")
setwd(path)  

# Read all files into giant list
files <- list.files(path=path)
rm <- grep(files,pattern = ".txt") # Don't load the text table summary
files <- files[-rm]
files <- files[grep(files, pattern = ".RData")] #only load rdata files
files <- mixedsort(files) # IMPORTANT: need this line to order in same order as scenario table!
results <- sapply(files, function(x) mget(load(x)),simplify=TRUE) # This is a giant list of all results
nscenarios <- length(results)
raw.table <- read.table("Scenario_Table.txt")  #This empty table is generated when the simulations are run - then you fill it in after everything runs
nsims <- nrow(results[[1]]$biomass.oneplus.true) # count nrows to know how many sims there are
years.test <- ncol(results[[1]]$biomass.oneplus.true) # count cols to know how many years there are
nyrs.to.use <- 100 # How many years you want to use to calculate all your metrics - There are no big differences btwn 50 and 100 yrs
calc.ind <- tail(1:years.test, nyrs.to.use) # Which years to calculate median depletion over (length = nyrs.to.use)

# Add performance measure columns to table
performance.measures <- c("LTmeancatch","LTnonzeromeancatch","SDcatch","n.5yrclose","n.10yrclose","nyrs0catch","meanbiomass","good4preds","SDbiomass","very.bad4preds","meanDepl","overallMaxCollapseLength","overallMaxBonanzaLength","BonanzaLength","CollapseLength","Prob.Collapse","Collapse.Severity","CV.Catch","Bonafide.Collapse")
pm.type <- c(rep("Fishery",times=6),rep("Ecosystem",times=11),"Fishery","Ecosystem") # for distinguishing types of PMs (mostly for plotting...)
raw.table[,performance.measures] <- NA


# ------------------------------------------------------------------------
# Summarize everything in one giant "outputs" table - this is ugly, sorry
# ------------------------------------------------------------------------

# Fxns for summarizing and plotting ---------------------------------------

all.summaries <- NA
all.summaries <- lapply(results,FUN = summ.tab, calc.ind = calc.ind, ou.ind = NA)   # This will take a while
all.summaries <- do.call(rbind.data.frame, all.summaries)
all.summaries$scenario <- rep(1:nscenarios,each=length(performance.measures))

# Match the scenarios to type of error, etc.
all.summaries <- merge(all.summaries,raw.table[,1:8],by="scenario")    # all.summaries is a giant table with 1080 rows = 72 scenarios * 15 PMs
all.summaries2 <- all.summaries %>% mutate(obs.error.type = recode_factor(obs.error.type, 
                                                                          "Tim"="Delayed change detection",
                                                                          "AC" = "Autocorrelated"),
                                           HCR = recode_factor(HCR, 'cfp' = 'Stability-favoring',
                                                               'constF' = 'Low F',
                                                               'C1' = 'Basic hockey stick',
                                                               'C2' = 'Low Blim',
                                                               'C3' = 'High Fmax',
                                                               'trend' = "Trend-based",
                                                               'constF_HI' = "High F"))

subset(all.summaries2, h==0.6 & obs.error.type=="Autocorrelated" & PM == "Collapse.Severity")
write.csv(all.summaries2, file = paste(Type,"_AllSummaries.csv",sep=""))

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

write.csv(raw.table, file=paste(Type,Sys.Date(),"_outputs.csv",sep=""))


#  Some plots -------------------------------------------------------------
palette <- brewer.pal(6,"Spectral")
hcr.colors <- palette[c(6,5,4,1,3,2)]
raw.table <- mutate(raw.table, HCR = recode_factor(HCR, 'cfp' = 'Stability-favoring',
                                                   'constF' = 'Low F',
                                                   'C1' = 'Basic hockey stick',
                                                   'C2' = 'Low Blim',
                                                   'C3' = 'High Fmax',
                                                   'constF_HI' = "High F"))

# Kite plots showing tradeoffs --------------------------------------------

nice.pms <- data.frame(original = colnames(raw.table)[-(1:8)],
                       polished = c("Mean catch","Mean nonzero catch",
                                    "Catch stability","Minimize \n P(5-yr closure)",
                                    "Minimize \n P(10-yr closure)","Minimize years \n w/ zero catch",
                                    "Mean biomass","Number of yrs \n above pred threshold",
                                    "SD(Biomass)","Number of yrs \n below low pred threshold",
                                    "Mean depletion","Max collapse length","Max bonanza length",
                                    "Bonanza length","Minimize \n collapse length", "Minimize \n P(collapse)",
                                    "Minimize \n collapse severity","CV(Catch)","Minimize \n long collapses"))

plotnames <- list()
tileplots <- list()
all.scaled <- vector()
plots <- data.frame("steepness"=0.6,"obs.error.type" = "AC","LH"=unique(raw.table$LH))

for(p in 1:5){
  tab <- subset(raw.table,obs.error.type == as.character(plots$obs.error.type[p]) & h == plots$steepness[p] & LH == plots$LH[p])
  tab.metrics <- tab[,-c(1:2,4:6,8)]
  crs <- tab[,"HCR"]
  remove.these <- c("n.10yrclose","SDbiomass","meanDepl","LTnonzeromeancatch","good4preds","very.bad4preds","CV.Catch","overallMaxCollapseLength","overallMaxBonanzaLength","Bonafide.Collapse")
  # Removed "Bonafide collapse" metric bc all CRs were performing similarly on it (in the paper this is called an "extended collapse")
  remove.ind <- which(colnames(tab.metrics) %in% remove.these)
  tab.metrics <- tab.metrics[-remove.ind]
  
  # Sometimes there won't be collapses, so for those turn NA's into zeroes (collapse length and severity)
  tab.metrics[,c("CollapseLength","Collapse.Severity")][is.na(tab.metrics[,c("CollapseLength","Collapse.Severity")])] <- 0
  #Here we deal with the PMs that are NEGATIVES (i.e., a high value for these is bad news) - I chose Option 3 below
  bad.pms <- c("SDcatch","n.5yrclose","n.10yrclose","nyrs0catch","SDbiomass","very.bad4preds","overallMaxCollapseLength","CollapseLength","Prob.Collapse","Collapse.Severity","CV(Catch)","Bonafide.Collapse")
  which.bad <- which(colnames(tab.metrics) %in% bad.pms)
  
  # For PMs with non-decimal values:
  tab.metrics[,c("SDcatch","CollapseLength")] <- 1 / tab.metrics[,c("SDcatch","CollapseLength")]
  tab.metrics$CollapseLength[is.infinite(tab.metrics$CollapseLength)] <- 1 # for when there are no collapses (get rid of infinite values)
  tab.metrics$SDcatch[is.infinite(tab.metrics$SDcatch)] <- 1 # there is also a case where SDcatch is zero because catches crash before the final 100 yrs
  
  # PMs with decimal values:
  tab.metrics[,c("Prob.Collapse","Collapse.Severity","n.5yrclose")] <- 1 - tab.metrics[,c("Prob.Collapse","Collapse.Severity","n.5yrclose")]
  tab.metrics[,c("nyrs0catch")] <- 1-(tab$nyrs0catch/nyrs.to.use) # This is essentially now the proportion of years when there *wasn't* 0 catch
  
  #tab.metrics[,-which.bad]
  props <- tab.metrics[,-(1:2)]
  props <- apply(props, MARGIN = 2,FUN = function(x) x/max(x,na.rm=T))
  all(props<=1) # check to make sure everything worked
  
  # OPTION #1: this scales everything between 0 and 1... which is kind of weird because it can make little differences look HUGE
  #tab.metrics[,which.bad] <- apply(tab.metrics[,which.bad],MARGIN = 2,FUN = function(x) 1-(x/ifelse(all(x==0),1,max(x,na.rm=T)))) # Turn all the "bad" PMs to their inverse
  #tab.metrics[,-which.bad] <- apply(tab.metrics[,-which.bad],MARGIN = 2,FUN = function(x) (x-min(x))/(max(x)-min(x)))
  
  # OPTION #2: if you want to scale to the best performance but not necessarily between 0 (worst) and 1 (best)
  # PRscaled = 1- min(1, PR / Z)
  
  # OPTION #3: scale differently for different metrics. 
  # 1/x for all the "bad" performance measures that come in numbers >= 1: PRscaled = 1/PR
  # SD catch
  # Collapse length
  # 1-x for all the "bad" performance measures that come in numbers <1: PRscaled = 1-PR
  # Probability of collapse
  # Probability of a 5-yr closure given that there was a closure
  # Proportion of years with zero catch
  # Collapse severity
  # if the whole column is zero for any of these, make all PRscaled=1
  
  final.tab <- data.frame(group = crs,props)
  test.nas <- apply(X = final.tab,FUN = anyNA,MARGIN = 2)
  na.metrics <- names(which(test.nas))
  
  if(length(na.metrics>0)){
    print(paste("The following performance metrics had NAs and were removed from the figure: ",na.metrics))
    final.tab <- final.tab[,-which(test.nas)]
    rm.metrics <- which(nice.pms$original %in% na.metrics)
    axis.labels <- nice.pms[-rm.metrics,'polished']
  }
  
  legend.presence <- ifelse(p == 1,TRUE,FALSE)
  axis.ind <- match(x = colnames(final.tab)[-1],table = nice.pms$original)
  axis.labels <- nice.pms$polished[axis.ind]
  final.tab$group <- factor(final.tab$group,levels=c("Basic hockey stick","Low Blim","High Fmax","High F","Low F","Stability-favoring"))
  plotnames[[p]] <- ggradar(final.tab,font.radar = "Helvetica",   # Add "_b" to fxn name if making w black background
                            grid.label.size=3,axis.label.size=8, 
                            legend.text.size = 4,
                            axis.labels = axis.labels,
                            plot.legend=FALSE,
                            palette.vec = hcr.colors[1:5],
                            manual.levels = levels(final.tab$group)[1:5])
  
  # ftm <- melt(final.tab,id.vars="group")
  # ftm$name <- ftm$variable
  # ftm$name <- factor(ftm$name, levels=c("LTmeancatch", "meanbiomass", "BonanzaLength","SDcatch","nyrs0catch","n.5yrclose","CollapseLength","Prob.Collapse","Collapse.Severity"))
  # ftm$group <- factor(ftm$group,c("Basic hockey stick","Low Blim","High Fmax","Low F","High F","Stability-favoring"))
  # ftm <- mutate(ftm,name = recode_factor(name, 'LTmeancatch' = "Mean catch",
  #                                        "meanbiomass" = "Mean biomass",
  #                                        "BonanzaLength" = "Bonanza Length",
  #                                        'SDcatch' = "Minimize \n catch variation",
  #                                        'nyrs0catch' = "Minimize \n years with 0 catch",
  #                                        'n.5yrclose' = "Minimize \n P(5 yr closure|closure)",
  #                                        "CollapseLength" = "Minimize \n collapse length",
  #                                        "Prob.Collapse" = "Minimize P(collapse)",
  #                                        "Collapse.Severity" = "Minimize collapse severity"))
  # tileplots[[p]] <- ggplot(ftm,aes(x=name,y=group)) + geom_tile(aes(fill=value)) + scale_fill_distiller(palette="RdYlBu",trans="reverse") + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5)) + geom_text(aes(label=round(value,digits = 1)))
  all.scaled <- rbind(all.scaled,final.tab)
}

pdf(file = paste(Type,Sys.Date(),"SENSITIVITY_KitePlots_v2.pdf",sep=""),width = 15,height=27,useDingbats = FALSE)
grid.arrange(plotnames[[1]],plotnames[[1]],plotnames[[2]],plotnames[[3]],plotnames[[4]],plotnames[[5]],ncol=2)
dev.off()

