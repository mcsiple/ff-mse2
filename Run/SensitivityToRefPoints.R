# Sensitivity to differences in reference points

# Libraries
library(plyr)

types <- c("Anchovy"="Anchovy",
           "Menhaden" = "Menhaden",
           "Sardine" = "Sardine")

# Parallelize so each simulation is on a different core:
# registerDoParallel(4)
# llply(.data=types,.fun = run.mod,.parallel = TRUE)
# stopCluster()

fftype = types[f]

#run.mod <- function(fftype){
# Set directories
basedir <- "~/Dropbox/Chapter4-HarvestControlRules/Code/ff-mse2"
resultsdir <- "~/Dropbox/Chapter4-HarvestControlRules/Results"
#subDir <- "Anchovy" # Name of ff type
subDir <- fftype

#Set up other simulation params
years.test = 250
nsims = 1000
tim.params = list(sigma0 = 0.2,tau0 = 0.1)
sig.s = 0.3 #0.30 (changed to 0.001 to check whether catches were limited at higher F by observation error)
R0.sens = NA #NO DYNAMIC R0 anymore-- ignore

# Load packages
library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)

# Load rev devs generating fxn, MSE main model, estimator fxns
#toplot=FALSE      # Don't plot examples of rec trajectories
source(file.path(basedir,"Recruitment/GenerateDevs.R")) 
source(file.path(basedir,"Estimators/CalcFTrue.R"))
source(file.path(basedir,"Run/HCR_Trajectory_NEW.R"))
source(file.path(basedir,"Estimators/Estimators.R"))
source(file.path(basedir,"Run/generate_M.R"))

# Load control files & set parameter values
#     Type     TargetSD
# 1  Sardine 0.5117551
# 2  Anchovy 0.3221512
# 3 Menhaden 0.3009440
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
  #Menhaden recruitment dev params
  recruit.sd <- 0.8
  recruit.rho <- 0.2
}


# Load harvest rules
source(file.path(basedir,"Control Rules/smith_oceana.R"))
source(file.path(basedir,"Control Rules/cfp.R"))
source(file.path(basedir,"Control Rules/hockey-stick.R"))
source(file.path(basedir,"Control Rules/trend-based-rule.R"))


# Scenarios
h = c(0.9, 0.6)
obs.error.type = c("AC","Tim","noerror")

HCR = c("cfp","constF","C1","C2","C3","constF_HI") # Took out trend because it was unrealistic-- but using trend in CPUE as adjustment (data-poor method) might be a good idea!
M.type = c("constant") # took out "regimeshift" and "time-varying" to save time but can be added back in for sensitivity

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

