# Testing to see whether move to Dropbox worked.
#' @title "stability-favoring" harvest control rule. Uses catch in the previous year, as well as traditional "hockey-stick" HCR parameters, to determine what the fishing rate will be in the following year.
#' \code{calc.F.cfp} returns the fishing rate, given biomass and catch in the previous year.
#' This function is essentially the same as other "hockey-stick" rules but does not allow the catch to change by more than 15%. It is a rule intended to stabilize catches.
#' @param prevCatch The catch in the previous year
#' @param Btru True biomass in the current year (vector of B_age - used to be single value)
#' @param Bobs Observed biomass in the current year
#' @param Btarget Biomass at which Fmax occurs - if B>Btarget, F=Fmax
#' @param Blim Biomass at which F = 0 - if B<Blim, F=0
#' @param Fmax Max fishing rate
#' @param lh List of life history traits, used for Baranov catch eqn
#' @param sel.at.age Matrix of ages, fishery selectivity
#' 
#' 
source("~/Dropbox/Chapter4-HarvestControlRules/Code/ff-mse2/Estimators/CalcFTrue.R") # need this to get F  that gives the proper catch
calc.F.cfp <- function(prevCatch, Bobs, Btru, Btarget, Blim, Fmax, lh = NA, sel.at.age = NA, sizes = NA){
  # prevCatch is total catch in yr-1
  # Bobs is a vector of biomass at age
  # Btru is is a vector of biomass at age
  # Btarget
  # Blim
  # Fmax
  f <- NA
  slope = Btarget/(Btarget-Blim)      # slope of the diagonal part of the hockey stick
  adj.constant <- Btarget/Fmax        # scales y axis to max fishing mortality
  if (sum(Bobs) <= Blim) {f = 0}
  if (sum(Bobs) > Blim & sum(Bobs) <= Btarget) {f <- slope * (sum(Bobs)-Blim) / adj.constant} 
                                      # adj.constant scales so the CR is linear btwn Blim and Btarget
  if (sum(Bobs) > Btarget) {f = Fmax}
  
  # Result of above is F from the basic hockey stick rule. Now find out if catch would change >15%
  possible.catch <- sum(Bobs *(1-exp(-(f*sel.at.age[,1]+lh$M)))*f*sel.at.age[,1] / (f*sel.at.age[,1] + lh$M) ) # Baranov catch eqn
  #print(c("Blim=",Blim,"Bobs=",sum(Bobs),"possible.catch=",possible.catch,"sum of Bobs > Blim?",sum(Bobs)>Blim,"PrevCatch=",prevCatch,"Btru",Btru))
  newcatch <- possible.catch
  # CHECK RESULTS
  if(sum(Bobs) > Blim){ # Per reviewer 2, this catch within 15% of previous catch requirement is sometimes only implemented above Btarget
      if(possible.catch != 0 && possible.catch < 0.85*prevCatch) {newcatch <- 0.85*prevCatch}
      if(possible.catch != 0 && possible.catch > 1.15*prevCatch) {newcatch <- 1.15*prevCatch}
      #print("Possible catch is >15% higher!")
      }else{newcatch <- possible.catch}
      if(newcatch > sum(Btru)){newcatch = sum(Btru) # is this problematique?
      #print("Suggested catch higher than biomass")
      }
  
  #print(c("new catch",newcatch))
  
  # Get the f from the control rule
  f.new <- calc.true.f(tac.fn = newcatch,M.fn = lh$M,sel.fn = sel.at.age[,1],Btrue.fn = Btru, w.at.age = sizes$weight.at.age[,1])
  return(f.new)
}

 # calc.F.cfp(prevCatch = 1000,
 #            Bobs = results[[1]]$biomass.oneplus.obs[2,191],
 #             Btru = results[[1]]$biomass.oneplus.true[2,191],
 #           Blim = 0.5*equilib$Bmsy,
 #           Btarget = equilib$Bmsy,
 #           Fmax = equilib$Fmsy,lh = lh.test,
 #           sel.at.age = lh.test$selectivity,
 #           sizes = sizes)

# Test
# First: get params from one of the ff types
        # basedir <- "~/Dropbox/Chapter4-HarvestControlRules/Code/ff-mse2"
        # source("/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Code/ff-mse2/Ctl/Anchovy_LHControl.R")
        # source(file.path(basedir,"Estimators/Estimators.R"))
        #  test.b <- c(7021.04, 4706.28, 3154.44, 2113.87, 1417.52,  950.28,  636.53)/10
        #  test.bobs <- test.b #(for now)
        #  sum(test.b)
        #  calc.F.cfp(prevCatch = 100, Bobs = test.bobs,Btru=test.b,Btarget = 1000,Blim = 100,Fmax = 0.6,lh = lh.test,sel.at.age = lh.test$selectivity)

        # test.frame <- data.frame(Bcurr = 300,  #seq(10,10000,by=100)
        #                          Btarget = 2000,
        #                          Blim = 100,
        #                          Fmax = 0.6,
        #                          prevCatch = seq(0,1000,by=50),
        #                          calc.f = NA)

# for(i in 1:nrow(test.frame)){
#   test.frame$calc.f[i] <- calc.F.cfp(prevCatch = test.frame$prevCatch[i],
#                                      Bt = test.frame$Bcurr[i],
#                                      Btarget = test.frame$Btarget[i],
#                                      Blim = test.frame$Blim[i],
#                                      Fmax = test.frame$Fmax[i],lh = lh.test,
#                                      sel.at.age = lh.test$selectivity)
# 
# }
# plot(test.frame$prevCatch,test.frame$calc.f,type='l')

