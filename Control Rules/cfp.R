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
  f <- NA
  slope = Btarget/(Btarget-Blim)      # slope of the diagonal part of the hockey stick
  adj.constant <- Btarget/Fmax        # scales y axis to max fishing mortality
  if (sum(Bobs) <= Blim) {f = 0}
  if (sum(Bobs) > Blim & sum(Bobs) <= Btarget) {f <- slope * (sum(Bobs)-Blim) / adj.constant} 
                                      # adj.constant scales so the CR is linear btwn Blim and Btarget
  if (sum(Bobs) > Btarget) {f = Fmax }
  # possible.catch <- f*Bt # old possible.catch
  # Result of above is F from the basic hockey stick rule. Now find out if catch would change >15%
  possible.catch <- sum(Bobs *(1-exp(-(f*sel.at.age[,1]+lh$M)))*f*sel.at.age[,1] / (f*sel.at.age[,1] + lh$M) ) # Baranov catch eqn
  newcatch <- possible.catch
  if(possible.catch != 0 && possible.catch < 0.85*prevCatch) {newcatch <- 0.85*prevCatch}
  if(possible.catch != 0 && possible.catch > 1.15*prevCatch) {newcatch <- 1.15*prevCatch}  # I checked whether 30% change (instead of 15%) works similarly; it didn't make a huge difference in relative performance.
  if(newcatch > sum(Btru[-1])){newcatch = sum(Btru[-1]) * 0.5}
  # This used to allow catches not to increase either-this was causing issues so I changed it.
  # Get the f from the control rule
  f.new <- calc.true.f(tac.fn = newcatch,M.fn = lh$M,sel.fn = sel.at.age[,1],Btrue = Btru, w.at.age = sizes$weight.at.age[,1])
  #f <- newcatch / Bt
  return(f.new)
}

# calc.F.cfp(prevCatch = results[[1]]$total.catch[2,190],
#            Bt = results[[1]]$biomass.oneplus.obs[2,191],
#            Blim = 0.5*equilib$Bmsy, 
#            Btarget = equilib$Bmsy, 
#            Fmax = equilib$Fmsy,lh = lh.test,
#            sel.at.age = lh.test$selectivity,
#            sizes = sizes)

# Test
# First: get params from one of the ff types
# source("/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Code/ff-mse2/Ctl/Menhaden_LHControl.R")
#  test.b <- c(7021.0485, 4706.2807, 3154.4495, 2113.8762, 1417.5272,  950.2850,  636.5326)
#  calc.F.cfp(prevCatch = 100, Bt = test.b,Btarget = 1000,Blim = 100,Fmax = 0.6,lh = lh.test,sel.at.age = lh.test$selectivity)
# 
# test.frame <- data.frame(Bcurr = 300,  #seq(10,10000,by=100)
#                          Btarget = 2000,
#                          Blim = 100,
#                          Fmax = 0.6,
#                          prevCatch = seq(0,1000,by=50),
#                          calc.f = NA)
# 
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

