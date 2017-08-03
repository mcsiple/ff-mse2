# Testing to see whether move to Dropbox worked.
#' @title "stability-favoring" harvest control rule. Uses catch in the previous year, as well as traditional "hockey-stick" HCR parameters, to determine what the fishing rate will be in the following year.
#' \code{calc.F.cfp} returns the fishing rate, given biomass and catch in the previous year.
#' This function is essentially the same as other "hockey-stick" rules but does not allow the catch to change by more than 15%. It is a rule intended to stabilize catches.
#' @param prevCatch The catch in the previous year
#' @param Bt Biomass in the current year
#' @param Btarget Biomass at which Fmax occurs - if B>Btarget, F=Fmax
#' @param Blim Biomass at which F = 0 - if B<Blim, F=0
#' @param Fmax Max fishing rate
#' @param lh List of life history traits, used for Baranov catch eqn
#' @param sel.at.age Matrix of ages, fishery selectivity
#' 
#' 
source("~/Dropbox/Chapter4-HarvestControlRules/Code/ff-mse2/Estimators/CalcFTrue.R") # need this to get F  that gives the proper catch
calc.F.cfp <- function(prevCatch, Bt, Btarget, Blim, Fmax, lh = NA, sel.at.age = NA, sizes = NA){
  f <- NA
  slope = Btarget/(Btarget-Blim)      # slope of the diagonal part of the hockey stick
  adj.constant <- Btarget/Fmax        # scales y axis to max fishing mortality
  if (sum(Bt) <= Blim) {f = 0}
  if (sum(Bt) > Blim & sum(Bt) <= Btarget) {f <- slope * (sum(Bt)-Blim) / adj.constant} 
                                      # adj.constant scales so the CR is linear btwn Blim and Btarget
  if (sum(Bt) > Btarget) {f = Fmax }
  # possible.catch <- f*Bt # old possible.catch
  # Result of above is F from the basic hockey stick rule. Now find out if catch would change >15%
  possible.catch <- sum(Bt *(1-exp(-(f*sel.at.age[,1]+lh$M)))*f*sel.at.age[,1] / (f*sel.at.age[,1] + lh$M) ) # Baranov catch eqn
  newcatch <- possible.catch
  if(possible.catch >= prevCatch * 1.15) {newcatch <- prevCatch*1.15 }else{
  if(possible.catch != 0 && possible.catch < 0.85*prevCatch) {newcatch <- 0.85*prevCatch } }
  # Get the f from the control rule
  f <- calc.true.f(tac.fn = newcatch,M.fn = lh$M,sel.fn = sel.at.age[,1],Btrue = Bt, w.at.age = sizes$weight.at.age[,1])
  #f <- newcatch / Bt
  return(f)
}


# Test
test.b <- c(7021.0485, 4706.2807, 3154.4495, 2113.8762, 1417.5272,  950.2850,  636.5326)
 calc.F.cfp(prevCatch = 100, Bt = test.b,Btarget = 1000,Blim = 100,Fmax = 0.6,lh = lh.test,sel.at.age = lh.test$selectivity)
# calc.F.cfp(prevCatch = 35.3, Bt = 117.7,Btarget = 21,Blim = 100,Fmax = 0.6)

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

