#############################################################################
# Takes life history characteristics, observation error, a trajectory of recruitment, a fishing rate, and a number of years, and returns time series of biomass, catches, and abundance
# The calc.trajectory function is modified from Christine's getTrajectory() function. It calls different control rules, rec devs, and types of observation error to run the MSE.
# 
#############################################################################

#rm(list=ls()[-which(ls()=="workdir")]) # Clear workspace if you need it


# Calculate SBPR0 to use in stock recruit fn
getSBPR<-function(nat.mort, maturity, fecun, n.ages){
  #' @description Function to calculate spawners per recruit
  #' @param nat.mort - natural mortality
  #' @param maturity - maturity ogive
  #' @param fecun - fecundity matrix
  #' @param n.ages - number of age classes
  #' @return spawning biomass per recruit
  N <- rep(NA, n.ages)
  N[1] <- 1
  for(age in 2:n.ages)
    N[age] <- N[age-1]*exp(-nat.mort)
  N[n.ages] <- N[n.ages]/(1-exp(-nat.mort))
  SBPR <- sum(N*maturity*fecun)
  return(SBPR)
}

# Beverton-Holt stock recruit fn
bevHolt <- function(h, R0 = 1000, SBPR0, S){
  #' @description Function to calculate recruitment
  #' @param h - steepness
  #' @param R0 - recruitment in unfished population
  #' @param SBPR0 - spawning biomass produced by one recrut in its lifetime
  #' @param S - spawning biomass
  # alpha and beta are regular BevHolt parameters
  alpha <- SBPR0 * ((1-h)/(4*h))
  beta <- (5*h-1)/(4*h*R0)
  Rec <- S / (alpha + beta * S)
  return (Rec)
}

calc.trajectory <- function(lh, obs.cv = NULL, init, rec.dev, rec.ram=NA, F0, cr, years, hcr.type = "hockeystick",obs.type="LN",const.f.rate = NULL,equilib = NULL, buffer = NULL,steepness = 0.9, R0.traj = 0,tim.params=NULL,time.var.m=NA, sig.s = NA){
  #' @description Function to calculate a trajectory of biomass, catch, and surplus production
  #' @param lh is a list containing the following:
        # M - natural mortality (scalar)
        # selectivity - sel. of fishery (matrix; column 1 is ages, col 2 is selectivities)
        # l.at.age - lengths at age (vector)
        # w.at.age - weights at age (vector)
        # maturity - vector of maturities at each age (vector of length nages)
  #' @param init - initial age distribution  (vector of length nages)
  #' @param obs.cv - observation error (scalar)
  #' @param rec.dev - vector of recruitment deviations. These will be simulated and multiplied by the outputs of the recruitment function.
  #' @param rec.ram - for testing; put in recruitment estimates from a RAM stock assessment to see if number of collapses etc is realistic; default is 
  #' @param F0 - fishing in initial year
  #' @param cr - list of control rule parameters, Blim, Btarget, Fmax
  #' @param years - number of years to simulate (should be dictated by length of recruitment time series)
  #' @param hcr.type - type of HCR. Currently there are "hockeystick" and "constF"
  #' @param const.f.rate - constant fishing rate if hcr.type == "constF"
  
  #' @param equilib - list of equilbrium values for the HCR
  #' @param buffer - the buffer used in the 40-10 rule
  #' @param steepness - h parameter
  #' @param R0.traj - a time series of R0, if you are testing if R0 is time-varying
  #' @param tim.params - a list of parameters (sigma0, tau0) describing the level of belief that the assessment or surveyor has in the possibility of large peaks or collapses
  #' @param time.var.m - a time series of natural mortality (M) values for cases where M is time-varying
  #' @return list: obs biomass, true biomass, catches, recruitment, and a bundle of other stuff
  #' IMPORTANT: 
  #' biomass.true - vector of true biomass-at-age
  #' biomass - vector of observed biomass-at-age
  #' oneplus.biomass - vector of true one-plus biomass
  #' ** Function outputs:
  #'      biomass.oneplus.obs - total observed B1+
  #'      biomass.total.true - total true B1+
  #'      biomass.oneplus.true - true B1+
  #'      biomass.total.obs
  
  if(is.null(lh$M)){
    print("M is missing")
    break
  }
  
  if(is.null(lh$selectivity)){
    print("Selectivity is missing")
    break
  }
  
  if(is.null(lh$ages)){
    print("Age vector is missing")
    break
  }
  
  if(any(is.na(lh$l.at.age))){
    print("NAs in length-at-age")
    break
  }
  
  if(any(is.na(lh$w.at.age))){
    print("NAs in weight-at-age")
    break
  }
  
  # if(any(is.na(rec.dev))){
  #   print("NAs in recruitment deviation time series")
  #   break
  # }
  
  n.ages<-length(lh$ages)
  sizes <- list(length.at.age = matrix(nrow=n.ages,ncol=years),
                weight.at.age = matrix(nrow=n.ages,ncol=years))  
  
  # Fill all years with the constant length- and weight-at-age
  for(yr in 1:years){  
    sizes$length.at.age[1:n.ages,yr] <- lh$l.at.age  
    sizes$weight.at.age[1:n.ages,yr] <- lh$w.at.age 
  }
  
  # Initiate population size and catch at age
  popn <- catch.at.age <- biomass <- biomass.true <- oneplus.biomass  <- matrix(nrow=n.ages,ncol=years+1)
  sp <- catch <- rep(NA, years)
   spawners <- sp.true <- R <- fishing <- intended.f <-  rep(NA, years+1)
  obs.err <- obs.epsR <- ac.obs <- rep(NA,years+1)
  
  # Initialize starting values
  pop.curr <- pop.next <- init
  biomass.true[,1] <- pop.curr*sizes$weight.at.age[,1]
  fishing[1] <- F0 # This is the true fishing rate
  intended.f[1] <- F0  # This is what the managers *think* F is (i.e., this is the F from the control rule) - 
  popn[,1] <- pop.curr 
  # Calculate sbpr
  sbpr <- getSBPR(lh$M, lh$maturity,fecun = lh$w.at.age, n.ages)   #lh$maturity is the maturity at age - it's a vector of length n.ages. Need to remember why fecundity at age = weight at age...
  if(all(!is.na(time.var.m))){sbpr <- getSBPR(time.var.m[1], lh$maturity,fecun = lh$w.at.age, n.ages)}
  
  # Selectivity
  sel.at.age <- matrix(ncol=years,nrow=nrow(lh$selectivity))
  for(i in 1:years){
    sel.at.age[,i] <- lh$selectivity[,2]   #lh$selectivity is a nages x 2 matrix, the first column is  ages, 2nd column is selectivity at age. 
                                            # This case is constant selectivity. 
  }
  
  # Start year loop
  for(yr in 1:(years)){ 
    popn[,yr] <- pop.curr   # pop.curr is a vector of nums at age
    
    # ESTIMATOR ##########################################
    if(obs.type == "LN"){
      biomass[,yr] <- add.LN.error(biomass.true = biomass.true[,yr], obs.cv = obs.cv, years = 1)$biomass.obs  #determine obs biomass based on true b and cv (this can be any error function)
    }   
    
    
    # Autocorrelated error  ###########################
    # Set sig.s and rho: These values are best estimates from Wiedenmann et al. 2015 Table 5: Median estimates of sd and autocorrelated in biomass observation error. For high steepness, slightly lower rho and sig.s
    #sig.s = ifelse(steepness>0.5,0.30,0.38) 
          rho = 0.5 # 0.5 #ifelse(steepness>0.5,0.82,0.87)
        if(is.na(sig.s)){ # if sig.s isn't provided at the beginning of the fxn, provide it here. sig.s=0.51 is for sardine.
          sig.s = 0.3 # This value comes from running the delay function a bunch of times, getting a target sd(log). This matches the sd from the autocorrelated error to the delay detection error.
    }
    if(obs.type == "AC"){
      eps.prev = ifelse(yr==1,1,eps.prev) # Initialize epsilon. Need to do this every time this type of error is used. #2 is NEW
      outs <-  add.wied.error(biomass.true = biomass.true[,yr],
                              epsilon.prev = eps.prev, 
                              sig.s =  sig.s, rho = rho)
      biomass[,yr] <- outs$biomass.est
      eps.prev <- outs$epsilon.curr
    }
          
    if(obs.type == "Tim"){
      sigma0 = tim.params$sigma0
      tau0 = tim.params$tau0
      # Add random error to first year, bc no prior information 
      # (i.e., in the first year, the "estimate" is just the biomass + some LN observation error):
      if(yr==1){biomass[,yr] <- biomass.true[,yr] * exp(rnorm(n.ages,0,sigma0)) 
                } else{
        biomass[,yr] <- tim.assessment(Eprev = biomass[,yr-1],
                                      B = biomass.true[,yr],
                                      sigma0 = sigma0,
                                      tau0 = tau0)
                }} 
          
      if(obs.type == "noerror"){
        biomass[,yr] <- biomass.true[,yr]
      }
    
    # HARVEST CONTROL RULE ################################
    if(hcr.type=="constF"){ 
      if(is.na(const.f.rate)){print("Constant fishing rate not provided in inputs")}else{
        imp.rate <- const.f.rate}
       } # true fishing rate based on observed biomass
    if(hcr.type=="cfp"){
      imp.rate <- calc.F.cfp(prevCatch = ifelse(yr==1,F0*sum(biomass[,yr]),catch[yr-1]),   #If it's year 1, use F0. Otherwise, determine f from hockey stick and previous year's catch, à la CFP's rule
                                Bobs = biomass[,yr], 
                                Btru = biomass.true[,yr], #Total biomass
                                Blim = 0.5*equilib$Bmsy, 
                                Btarget = equilib$Bmsy, 
                                Fmax = equilib$Fmsy,
                                lh = lh,
                                sel.at.age = sel.at.age,
                                sizes = sizes)
      
    }
    if(hcr.type=="C1"){  
      imp.rate <- calc.F.oceana(Bt = sum(biomass[,yr]),
                                   Blim = 0.4*equilib$B0,
                                   Btarget = 0.8*equilib$B0,
                                   M = lh$M)            # Fmax will be 0.5*lh$M
    }
    if(hcr.type=="C2"){ # This is C1 but with Low Blim -- Fmax = 0.5M
        imp.rate <- calc.F.oceana(Bt = sum(biomass[,yr]),
                                    Blim = 0.1*equilib$B0, 
                                    Btarget = 0.8*equilib$B0, 
                                    M = lh$M)
    }
    
    if(hcr.type=="C3"){ # This is C1 but with Hi Fmax
        imp.rate <- calc.F.stick(Bt = sum(biomass[,yr]),
                                     Blim = 0.4*equilib$B0, 
                                     Btarget = 0.8*equilib$B0, 
                                     Fmax = equilib$Fmsy
                                     )
      }
    if(hcr.type=="trend"){
      imp.rate <- trend.rule(B.2yr = ifelse(yr %in% 1:2,sum(biomass[,yr]),sum(biomass[,yr-2])),
                                B.prev = ifelse(yr == 1,sum(biomass[,yr]), sum(biomass[,yr-1])),
                                B.curr = sum(biomass[,yr]), 
                                F.const = const.f.rate)
    }
    
    # Get total catches
    tac <- biomass[,yr] *(1-exp(-(imp.rate*sel.at.age[,1]+lh$M)))*imp.rate*sel.at.age[,1] / (imp.rate*sel.at.age[,1] + lh$M) 
    
    # USE BARANOV CATCH EQUATION TO GET *TRUE F*
    if(obs.type=="noerror"){fishing[yr] = imp.rate}else{
        if(imp.rate==0){fishing[yr] = 0} else{
          fishing[yr] <- calc.true.f(tac.fn = tac,M.fn = lh$M,sel.fn = sel.at.age[,yr],Btrue = biomass.true[,yr], w.at.age = sizes$weight.at.age[,1]) # Double check that selectivity shouldn't be zero for age 0 fish
        } }
        intended.f[yr] <- imp.rate # This is what the 'managers' think the F is each year
    ##########
    death.rate <- lh$M + sel.at.age[,yr] * fishing[yr] 
    if(all(!is.na(time.var.m))){death.rate <- time.var.m[yr] + sel.at.age[,yr] * fishing[yr]}
    
    catch.at.age[,yr] <- (sizes$weight.at.age[,yr] * sel.at.age[,yr] * fishing[yr]*(1-exp(-death.rate)) * popn[,yr])/death.rate 
    catch[yr]<-sum(catch.at.age[,yr])  #Total catch
    next.year.S <- sum(pop.curr*lh$maturity*sizes$weight.at.age[,yr])
    
    #Recruitment fed into the model comes from real or simulated recruitment estimates
    R0 <- ifelse(length(R0.traj)>1 , R0.traj[yr], lh$R0) # R0 is usually set by life history traits. If it's NOT, make it whatever it should be in that year, based on R0.traj 
    
    pop.next[1] <- bevHolt(h = steepness, S = next.year.S, SBPR0 = sbpr, R0 = R0)*rec.dev[yr]  # Beverton-Holt recruitment * recruitment deviations
    if(!is.na(rec.ram[1])){pop.next[1] <- rec.ram[yr]}
    
    if(all(!is.na(time.var.m))){
      sbpr <- getSBPR(time.var.m[yr], lh$maturity,fecun = lh$w.at.age, n.ages)  
      pop.next[1] <-bevHolt(h = steepness, S = next.year.S, SBPR0 = sbpr, R0 = R0)*rec.dev[yr] 
      if(!is.na(rec.ram[1])){pop.next[1] <- rec.ram[yr]}
      }
    
    if(is.na(pop.next[1])){   
      print("Recruitment is NA for some reason! Below are parameters")
                                  print(c(yr,next.year.S,sbpr,R0,rec.dev[yr]))
                                  }
    
    R[yr] <- pop.next[1]
    survival <- pop.curr*exp(-death.rate)
    
    pop.next[2:n.ages] <- survival[1:(n.ages-1)]
    pop.next[n.ages] <- pop.next[n.ages] + survival[n.ages]  #Plus group pop size
    pop.curr <- pop.next
    biomass.true[,yr+1] <- pop.curr*sizes$weight.at.age[,yr] # These are from the original model
    oneplus.biomass[-1,yr+1] <- pop.curr[-1]*sizes$weight.at.age[-1,yr] # These are from the original model
    sp.true[yr] <- sum(biomass.true[,yr+1]) - sum(biomass.true[,yr]) + catch[yr]      # Surplus production

    } # end of years loop

  for(yr in 1:(years-1)){
    sp[yr] <- sum(biomass[,yr+1]) - sum(biomass[,yr]) + catch[yr]    # Like observed surplus production (depracated)
  }
  

  # GETMSY CONDITION - If calc.trajectory is being used in the getEquilibriumConditions() function
  if(is.null(equilib) && length(R0.traj) == 1 ){B0.traj = rep(1,times=years) } # Basically, depletion in this case just has to get filled in; values don't mean anything

  # ALL OTHER CONDITIONS - if you are running the full MSE
  if(!is.null(equilib) && length(R0.traj) > 1){B0.traj <- (R0.traj[1:years])*0.110548} # weird linear relationship between R0 and B0
  if(!is.null(equilib) && length(R0.traj) == 1 ){B0.traj <- rep(equilib$B0,times=years)}
  
return(list(popn=popn[,1:years], 
            biomass.oneplus.obs=apply(biomass[-1,1:years],2,sum,na.rm=T), # observed one plus biomass of the population
            sp=sp, 
            catch.at.age=catch.at.age[,1:years], 
            total.catch=catch, 
            fishing=fishing[1:years], 
            biomass.total.true=apply(biomass.true[,1:years],2,sum,na.rm=T),         # biomass.true is the true total biomass
            sp.true=sp.true[1:years],
            wts=sizes$weight.at.age,lengths=sizes$length.at.age,
            "rec"=R[1:years], intended.f=intended.f[1:years],
            biomass.oneplus.true=apply(biomass.true[-1,1:years],2,sum,na.rm=T)))    # oneplus.biomass is the true one-plus biomass
} #end of calc.trajectory function



# The function for calculating MSY (uses getTrajectory, so has to be loaded last)
source("~/Dropbox/Chapter4-HarvestControlRules/Code/ff-mse2/Run/MSY_Fxns.R")


# testie2 <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "cfp",equilib = equilib,steepness=steepness,obs.type = "AC", tim.params = tim.params,const.f.rate=0.6, sig.s = .3,rec.ram=NA)
# testie3 <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "C2",equilib = equilib,steepness=steepness,obs.type = "AC", tim.params = tim.params,const.f.rate=0.6, sig.s = .3,rec.ram=NA)
############################### ############################### ###############################
############################### ############################### ###############################
############################### ############################### ###############################

# For testing:
# basedir <- "/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Code/ff-mse2"
# source(file.path(basedir,"Recruitment/GenerateDevs.R")) 
# source(file.path(basedir,"Estimators/CalcFTrue.R"))
# source(file.path(basedir,"Estimators/Estimators.R"))
# source(file.path(basedir,"Run/generate_M.R"))
# # Load harvest rules
# source(file.path(basedir,"Control Rules/smith_oceana.R"))
# source(file.path(basedir,"Control Rules/cfp.R"))
# source(file.path(basedir,"Control Rules/hockey-stick.R"))
# source(file.path(basedir,"Control Rules/trend-based-rule.R"))
# #source(file.path(basedir,"Ctl/Sardine_LHControl.R"))
# #source(file.path(basedir,"Ctl/Sardine_FisheryControl.R"))
# source(file.path(basedir,"Ctl/Anchovy_LHControl.R"))
# source(file.path(basedir,"Ctl/Anchovy_FisheryControl.R"))
# 
# 
# # # # Menhaden recruitment dev params
# recruit.sd <- 0.8
# recruit.rho <- 0.2
# steepness = 0.6
# ages.test <- 0:6
# nages.test <- length(ages.test)
# time.var.m=NA
# toplot = FALSE
# years.test <- 250
# tim.params <- list(sigma0 = 0.2,tau0 = 0.1)
# obs.cv = 1.2
# #init = init.test
# #F0 = F0.test
# cr = NA
#       years = years.test
#       hcr.type = "cfp"
#       equilib = getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=steepness)
#       obs.type = "AC"
#       tim.params = tim.params
#       #const.f.rate=0.6
#       sig.s = .3
#       rec.ram=NA
# 
#       # tot <- sum(c(6.221,4.170,2.795,1.873,1.256,0.842,0.564))
#       # init.prop <- c(6.221,4.170,2.795,1.873,1.256,0.842,0.564)/tot
#       # init.B <- init.prop*200000 # Biomass at age -- 110.54 is the B0 when R0 = 1000
#       # init.test <- init.B / lh.test$w.at.age # These are all the numbers at age (spawning and non-spawning biomass)
#       # init = init.test
#       # F0.test <- 0.3 # This is arbitrary
# 
#       par(mfrow=c(3,3))
#       #set.seed(123)
#       for(i in 1:9){
#       rec.dev.test <- generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd)
#       rec.dev = rec.dev.test
#       testie4 <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0, cr = cr, years = years.test,hcr.type = hcr.type,equilib = equilib,steepness=steepness,obs.type = obs.type, tim.params = tim.params,const.f.rate=equilib$Fmsy, sig.s = .3,rec.ram=NA,time.var.m = NA)
# #
# #      par(mfrow=c(2,1))
# plot(testie4$biomass.oneplus.true,type='l',ylim=c(0,6e4))
# lines(testie4$total.catch,col='red')
# }
# any(testie4$total.catch == 0)
# plot(testie4$intended.f,type='l')
# lines(testie4$fishing,col='red')
# which(testie4$biomass.oneplus.true < testie4$total.catch)
