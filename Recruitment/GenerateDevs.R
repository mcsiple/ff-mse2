
cv <- function(ts){
  mean = mean(ts)
  sd = sd(ts)
  CV <- sd/mean
  return(CV)
}

# A la Thorson et al.  ----------------------------------------------------

generate.devs <- function(N, rho, sd.devs, burnin=100, plot=FALSE){
  #' @description Function to generate recruitment deviations using the "time series" method 
  #' @param N - number of years for which to generate recruitment deviations 
  #' @param rho - autocorrelation
  #' @param sd.devs - standard deviation of the recruitment devs
  
  #Generate noise
  innov <- rnorm(N+burnin,mean = (-sd.devs^2)/2,sd = sd.devs) # Square root of 1-rho corrects it so it doesn't just wander unidirectionally
  
  #Make a vector to hold simulated deviations
  dev.ts <- numeric(length=N+burnin)
  for(yr in 2:(N+burnin)){
    dev.ts[yr] <- rho * dev.ts[yr-1] + sqrt(1-rho^2) * innov[yr] #The 0.7 just keeps it low... I think this needs to be fixed
  }
  
  #Trim off burn-in years (this is so we get the stationary part, if it's stationary)
  dev.ts <- dev.ts[-(1:burnin)]
  dev.ts <- exp(dev.ts) #Since these are log deviations
  if(plot==TRUE){ plot(1:N, dev.ts, type='l',yaxt="n",ylab='',xlab='Year')
                  title(ylab="Recruit deviations", line=0.4, cex.lab=1.2)
                  CV <- round(cv(dev.ts),digits = 2)
                  text(x = 46,y = 0.99*max(dev.ts),labels = paste("CV =",CV))
                  }
  return(dev.ts)
}




# Test and plot the devs! -------------------------------------------------
              # if(toplot == TRUE){
              # 
              #         par(mfcol=c(3,3))
              #         # "Sardine"
              #         for(i in 1:3){
              #           generate.devs(N=50,rho = 0.9,sd.devs = 0.1,burnin = 100,plot = TRUE)
              #           if(i==1){mtext(side = 3,text = "Sardine-like",outer = FALSE)}
              #           }
              #         
              #         # "Herring/anchovy" 
              #         for(i in 1:3){
              #           generate.devs(N=50,rho = 0.5,sd.devs = 0.3,burnin = 100,plot = TRUE)
              #           if(i==1){mtext(side = 3,text = "Herring/anchovy-like",outer = FALSE)}
              #         }
              #         
              #         # "Menhaden"
              #         for(i in 1:3){
              #           generate.devs(N=50,rho = 0,sd.devs = 0.7,burnin = 100,plot = TRUE)
              #           if(i==1){mtext(side = 3,text = "Menhaden-like",outer = FALSE)}
              #         }
              #   
              #   # Now plot unfished biomass with the rec devs -- ****fix this part to show biomass over time
              #           # Load all the contents of the control rules folder:
              #           setwd("/Users/mcsiple/Dropbox/Chapter 4 - Harvest Control Rules/Code/Control Rules")
              #           lapply(list.files(pattern = "[.]R$", recursive = TRUE), source)
              #           par(mfcol=c(3,3))
              #         
              #         
              #         # "Sardine"
              #         for(i in 1:3){
              #           devs <- generate.devs(N=200,rho = 0.9,sd.devs = 0.1,burnin = 100,plot = FALSE)
              #           biomass <- calc.trajectory(lh = lh.test, obs.cv = NULL, init = init.test, rec.dev = devs, F0, cr, years = years.test, hcr.type = "hockeystick",obs.type="LN",const.f.rate = 0.2,equilib = NULL, buffer = NULL,steepness = 0.9, R0.traj = 0,tim.params=NULL)$biomass
              #           
              #           if(i==1){mtext(side = 3,text = "Sardine-like",outer = FALSE)}
              #         }
              #         
              #         # "Herring/anchovy" 
              #         for(i in 1:3){
              #           generate.devs(N=200,rho = 0.5,sd.devs = 0.3,burnin = 100,plot = FALSE)
              #           if(i==1){mtext(side = 3,text = "Herring/anchovy-like",outer = FALSE)}
              #         }
              #         
              #         # "Menhaden"
              #         for(i in 1:3){
              #           generate.devs(N=200,rho = 0,sd.devs = 0.7,burnin = 100,plot = FALSE)
              #           if(i==1){mtext(side = 3,text = "Menhaden-like",outer = FALSE)}
              #         }
              # }



# Combine with population model to get time series of recruitment ---------
# source("~/Dropbox/Chapter 4 - Harvest Control Rules/Code/MSY_Fxns.R")
# source("~/Dropbox/Chapter 4 - Harvest Control Rules/Code/HCR_Trajectory.R")
# 
# 
# testA <- calc.trajectory(lh = lh.test,obs.cv = 0.75, init = init.test, rec.dev = test1, F0 = 0, cr = NULL, years = 200,hcr.type = "constF", const.f.rate = 0, equilib = NULL, steepness = 0.9,obs.type = "AC")
# 



##################################################################
##################################################################

# Function used to generate recruitment deviations
generate.devs.OLD <- function(N, rho, sd.devs, burnin=100, plot=FALSE){
  # N == length of the series of recruitment deviations
  # dev0 == the deviation in the first year
  # rho == the autocorrelation
  # sd.devs == the standard deviation of the recruitment devs
  #Generate noise
  log.innov <- rnorm(N+burnin,mean = 0,sd = sd.devs*sqrt(1-rho)) # Square root of 1-rho corrects it so it doesn't just wander unidirectionally
  innov <- exp(log.innov)
  
  #Make a vector for simulated rec devs
  dev.ts <- numeric(length=N+burnin)
  for(yr in 2:(N+burnin)){
    dev.ts[yr] <- rho * dev.ts[yr-1] + 0.7*innov[yr] #The 0.7 just keeps it low... I think this needs to be fixed
  }
  
  #Trim off burn-in years (this is so we get the stationary part, if it's stationary)
  dev.ts <- dev.ts[-(1:burnin)]
  if(plot==TRUE){ plot(1:N, dev.ts, type='l')}
  return(dev.ts)
}

# 
# pdf("~/Dropbox/Chapter4-HarvestControlRules/Figures/RecruitmentDevsExample.pdf",height = 3, width = 10)
# par(mfcol=c(1,3))
# # "Sardine"
#   generate.devs(N=50,rho = 0.9,sd.devs = 0.1,burnin = 100,plot = TRUE)
#   mtext(side = 3,text = "Sardine-like",outer = FALSE)
# 
# # "Herring/anchovy" 
#   generate.devs(N=50,rho = 0.5,sd.devs = 0.3,burnin = 100,plot = TRUE)
#   mtext(side = 3,text = "Herring/anchovy-like",outer = FALSE)
# 
# # "Menhaden"
# 
#   generate.devs(N=50,rho = 0,sd.devs = 0.7,burnin = 100,plot = TRUE)
#   mtext(side = 3,text = "Menhaden-like",outer = FALSE)
# dev.off()
# 
