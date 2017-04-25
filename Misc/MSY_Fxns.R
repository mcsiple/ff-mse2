#' @description function which returns Bmsy, Fmsy given life history parameters
#' @param lh = list of life history parameters for one stock
#' @param fish = sequence of all possible fishing values
#' @param years = number of years to run at each fishing value
#' @return a list equil with elements fmsy, bmsy, binf

getEquilibriumConditions<-function(lh, fish, years, steepness){
  # Run get trajectory function across fishing vector
  # and return final biomass and surplus production
  if(is.null(steepness)){print("Steepness value is missing")}
  equilib <- sapply(fish, function(x)
    with(calc.trajectory(lh=lh, obs.cv = 0, init = rep(1000,length(lh$ages)), 
                         rec.dev = rep(1,times=years), F0=x, cr=NULL, 
                         years=years, hcr.type = "constF",const.f.rate = x,
                         equilib=NULL,steepness=steepness),
    c(x,tail(biomass, 1), tail(sp, 2)[1],total.catch[years])))  # tail(sp, 2)[1] is because the last year of surplus production is NA (because it's calculated from B[t+1]) so it's hacky but it's fine.
  
  # The value of F which maximizes surplus production is Fmsy
  f.msy <- fish[which.max(equilib[3,])]
  # The value of B which maximizes biomass is Binf
  b.msy <- equilib[2, which.max(equilib[3,])]
  b.inf <- equilib[2, which.max(equilib[2,])]
  
  unfished <- calc.trajectory(lh, obs.cv = 0, init = rep(10000,length(lh$ages)), 
                              rec.dev = rep(1,times=years), F0=0, cr=NULL, 
                              years=years, hcr.type = "constF",const.f.rate = 0, equilib=NULL,steepness=steepness)
  b.0 <- tail(unfished$biomass,1)

  #Plot equilibrium yield curve
  #plot(equilib[2,]/b.inf,equilib[4,],type="l",xlim=c(0,1),main="Final year catch",xlab="Relative depletion",ylab="Catch (t)")
  
  return(list("Fmsy"=f.msy,"Bmsy"=b.msy,"Binf"=b.inf,"B0" = b.0))
}


# years.test <- 150
# getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 10,steepness=0.9)
# getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=0.5)
