# This code is to figure out whether the rec devs we made up generate similar biomass dynamics to real populations, without fishing. I.e., do they seem to generate the same number and frequency of collapses without fishing, as real populations would experience w/o fishing? I fed RAM recruitment estimates straight into the age-structured model to see.

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
id.table <- unique(x[,c('ASSESSID','commonname')])
id.table$type <- NA
# Anchovy and herring are together - I know this is a little weird
 id.table$type[grep(pattern = "anch",ignore.case = T,x = id.table$commonname)] <- "Anchovy"
 id.table$type[grep(pattern = "herring",ignore.case = T,x = id.table$commonname)] <- "Anchovy"
 
 id.table$type[grep(pattern = "menhaden",ignore.case = T,x = id.table$commonname)] <- "Menhaden"
 
 id.table$type[grep(pattern = "sardine",ignore.case = T,x = id.table$commonname)] <- "Sardine"
 id.table$type[grep(pattern = "pilchard",ignore.case = T,x = id.table$commonname)] <- "Sardine"

na.types <- id.table$ASSESSID[which(is.na(id.table$type))] # these stocks are removed because they're types we're not including in the analysis (.e.g, mackerel)
 
bigframe <- vector()
for(i in 1:nstocks){
  stock <- ffdat %>% subset(ASSESSID == stocks[i] & !is.na(R)) %>% as.data.frame()
  type <- id.table$type[which(id.table$ASSESSID == stock$ASSESSID[1])]
  if(is.na(type)){next}
  print(type)
  rec.ram.test <- stock$R
  years.to.plot <- length(rec.ram.test)
  actual.years <- stock$Years
  steepness = 0.6
  obs.type <- "AC"
  HCR <- "constF"
  # Get LH characteristics for that ff type
            # Sardines
            if(type == "Sardine"){
              source(file.path(basedir,"Ctl/Sardine_LHControl.R"))
              source(file.path(basedir,"Ctl/Sardine_FisheryControl.R"))
            }else{
            # Anchovy/Herring
            if(type == "Anchovy"){
              source(file.path(basedir,"Ctl/Anchovy_LHControl.R"))
              source(file.path(basedir,"Ctl/Anchovy_FisheryControl.R"))
            }else{
            # Menhaden
            if(subDir == "Menhaden"){
              source(file.path(basedir,"Ctl/Menhaden_LHControl.R"))
              source(file.path(basedir,"Ctl/Menhaden_FisheryControl.R"))
            }
            }}
  equilib = getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=steepness)
  #rec.dev.test <- generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd)
  test.constF <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = NA,rec.ram = rec.ram.test, F0 = F0.test, cr = cr.test, years = years.to.plot,hcr.type = "constF", const.f.rate = 0, steepness = steepness,obs.type = obs.type,equilib=equilib,R0.traj = R0.sens, tim.params = tim.params,time.var.m = NA)
  nofish <- melt(test.constF[-c(1,4,9,10)])
  nofish$year <- rep(actual.years,times=length(unique(nofish$L1)))
  nofish$stock <- stocks[i]
  mf <- subset(nofish,L1=="biomass.oneplus.true")
  mf$thresh <- (mean(mf$value))*0.2           #"Collapse" threshold - always based on unfished biomass
  mf$bm <- mean(mf$value)                     # mean biomass, in case you want to plot it
  mf$type <- type
  test.coll <- mf$value < mf$thresh
  coll.years <- length(which(test.coll))      # Total # years collapsed
  y <- rle(test.coll)
  n.colls <- length(which(y$values==TRUE))    # number of "collapse" periods (consecutive years with low biomass)
  mf$ncolls <- n.colls
  mf$collyears <- coll.years
  bigframe <- rbind(bigframe,mf)
}

# For adding summary stats to plots
labels <- unique(bigframe[,c("stock","ncolls","collyears")])
labels <- bigframe %>% group_by(stock, ncolls, collyears) %>% summarize(ycoord = max(value) * 0.9) %>% as.data.frame()
print(labels)

baseplot <- ggplot(bigframe,aes(x=year,y=value,colour=type)) + geom_line() + geom_hline(aes(yintercept = thresh),col="red")  + 
  scale_color_brewer(type="qual",palette=2) +
  geom_hline(aes(yintercept = bm,colour=type),alpha=0.2) +
  facet_wrap(~stock,scales = "free_y") + ylab("1+ Biomass") + 
  xlab("Year") + theme_classic(base_size = 12) 

# Plot thresholds and stuff - doesn't make a difference what the obs model is.
pdf("B_oneplus_true_RAM_AC.pdf",width=14,height =9)
baseplot + 
  geom_label(x = 2000, aes(y=ycoord, label = ncolls),size=3,colour="black",data=labels) +
  theme(strip.background = element_rect(colour="white"))
dev.off()

