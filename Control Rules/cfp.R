# Testing to see whether move to Dropbox worked.
#' @title "stability-favoring" harvest control rule. Uses catch in the previous year, as well as traditional "hockey-stick" HCR parameters, to determine what the fishing rate will be in the following year.
#' \code{calc.F.cfp} returns the fishing rate, given biomass and catch in the previous year.
#' This function is essentially the same as other "hockey-stick" rules but does not allow the catch to change by more than 15%. It is a rule intended to stabilize catches.
#' @param prevCatch The catch in the previous year
#' @param Bt Biomass in the current year
#' @param Btarget Biomass at which Fmax occurs - if B>Btarget, F=Fmax
#' @param Blim Biomass at which F = 0 - if B<Blim, F=0
#' @param Fmax Max fishing rate

calc.F.cfp <- function(prevCatch, Bt, Btarget, Blim, Fmax){
  f <- NA
  slope = Btarget/(Btarget-Blim)      # slope of the diagonal part of the hockey stick
  adj.constant <- Btarget/Fmax        # scales y axis to max fishing mortality
  if (Bt <= Blim) {f = 0}
  if (Bt > Blim & Bt <= Btarget) {f <- slope * (Bt-Blim) / adj.constant}
  if (Bt > Btarget) {f = Fmax }
  possible.catch <- f*Bt
  if(possible.catch >= prevCatch * 1.15) {newcatch <- prevCatch*1.15 }
  if(possible.catch != 0 && possible.catch < 0.85*prevCatch) {newcatch <- 0.85*prevCatch } # This used to say Bt
  else(newcatch <- possible.catch)
  #if(newcatch > Bt) (newcatch = Bt) # This step shouldn't be necessary, so I have to figure out why it's happening!
  f <- newcatch / Bt
  return(f)
}


# Test
# calc.F.cfp(prevCatch = 39.6, Bt = 95.8,Btarget = 21,Blim = 100,Fmax = 0.6)
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
#                                      Fmax = test.frame$Fmax[i])
# }
# plot(test.frame$prevCatch,test.frame$calc.f,type='l')
# plot(test.frame$Bcurr,test.frame$calc.f*test.frame$Bcurr, type='l')
