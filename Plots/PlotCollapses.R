# Set directories
                  basedir <- "/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Code/ff-mse2"
                  resultsdir <- "/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Results"
                  subDir <- "Sardine" # Name of ff type
                  
                  #Set up other simulation params
                  years.test = 250
                  nsims = 1000
                  tim.params = list(sigma0 = 0.2,tau0 = 0.1)
                  R0.sens = NA #NO DYNAMIC R0 anymore-- ignore
                  
                  # Load packages
                  library(reshape2)
                  library(ggplot2)
                  # Load rev devs generating fxn, MSE main model, estimator fxns
                  toplot=FALSE      # Don't plot examples of rec trajectories
                  source(file.path(basedir,"Recruitment/GenerateDevs.R")) 
                  source(file.path("/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Code/HCR_Trajectory.R"))
                  source(file.path(basedir,"Estimators/Estimators.R"))
                  source(file.path(basedir,"Run/generate_M.R"))
                  
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
                  

ffdat <- read.csv("/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Datasets/allforagedata.csv")

library(ggplot2)
library(plyr)
library(dplyr)

ggplot(ffdat,aes(x=Years,y=R)) + geom_line()+geom_point() +facet_wrap(~ASSESSID,scales="free_y")
tslength <- ffdat %>% group_by(ASSESSID) %>% subset(!is.na(R)) %>% summarize(tslength = max(Years)-min(Years),minyr=min(Years),maxyr=max(Years)) %>% as.data.frame()

max(tslength$tslength)
subset(ffdat,ASSESSID == "HERRNOSS")


#rec.ram.test <- subset(ffdat,ASSESSID == "HERRNOSS")$R
rec.ram.test <- subset(ffdat,ASSESSID == "CMACKPCOAST")$R
rec.ram.test <- rec.ram.test[-80]

# Test to see if collapse thresholds are similar for RAM species. Use rec estimates from RAM models directly, to see if they crash the population as frequently.
x <- ffdat %>% subset(!is.na(R)) 
stocks <- unique(x$ASSESSID)
nstocks <- length(stocks)
bigframe <- vector()
for(i in 1:nstocks){
  stock <- ffdat %>% subset(ASSESSID == stocks[i] & !is.na(R)) %>% as.data.frame()
  rec.ram.test <- stock$R
  years.to.plot <- length(rec.ram.test)
  steepness = 0.6
  obs.type <- "AC"
  HCR <- "constF"
  recruit.sd = .6 #scenarios$recruit.sd[1]
  recruit.rho = .9 #scenarios$recruit.rho[1]
  equilib = getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=steepness)
  #rec.dev.test <- generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd)
  test.constF <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = NA,rec.ram = rec.ram.test, F0 = F0.test, cr = cr.test, years = years.to.plot,hcr.type = "constF", const.f.rate = 0, steepness = steepness,obs.type = obs.type,equilib=equilib,R0.traj = R0.sens, tim.params = tim.params,time.var.m = NA)
  nofish <- melt(test.constF[-c(1,4,9,10)])
  nofish$year <- rep(1:years.to.plot,times=length(unique(nofish$L1)))
  nofish$stock <- stocks[i]
  mf <- subset(nofish,L1=="biomass.oneplus.true")
  mf$thresh <- (mean(mf$value))*0.2
  mf$bm <- mean(mf$value)
  bigframe <- rbind(bigframe,mf)
}

#ltp<- unique(bigframe[,c("stock","thresh","bm")]) # thresholds
pdf("B_oneplus_true_RAM.pdf",width=14,height =9)
ggplot(bigframe,aes(x=year,y=value)) + geom_line() + geom_hline(aes(yintercept = thresh),col="red")  + #geom_hline(aes(yintercept = bm),col="blue") +
  facet_wrap(~stock,scales = "free_y")  #+ xlim(c(150,250))
dev.off()