
getEquilibriumConditions <- function(lh, fish, years, steepness){
  #' @description function which returns Bmsy, Fmsy given life history parameters
  #' @param lh = list of life history parameters for one stock
  #' @param fish = sequence of all possible fishing values
  #' @param years = number of years to run at each fishing value
  #' @return a list equil with elements fmsy, bmsy, binf
  #
  # Run get trajectory function across fishing vector
  # and return final biomass and surplus production
  if(is.null(steepness)){print("Steepness value is missing")}
  equilib <- sapply(fish, function(x)
    with(calc.trajectory(lh=lh, obs.cv = 0, init = rep(1000,length(lh$ages)), 
                         rec.dev = rep(1,times=years), F0=x, cr=NULL,R0.traj = NA, 
                         years=years, hcr.type = "constF",const.f.rate = x,
                         equilib=NULL,steepness=steepness,obs.type = "noerror",rec.ram = NA,time.var.m = NA),
    c(x,tail(biomass.total.true,1), # B0
      tail(sp,2)[1],  # Equilibrium SP. tail(sp, 2)[1] is because the last year of surplus production is NA (because it's calculated from B[t+1]) so it's hacky but it's fine.
      tail(total.catch,1))))  # Total catch, for getting msy 
  #print(equilib)
  #The value of F which maximizes surplus production is Fmsy
  f.msy <- fish[which.max(equilib[3,])]
  #The value of B which maximizes biomass is Binf
  b.msy <- equilib[2, which.max(equilib[3,])]
  b.inf <- equilib[2, which.max(equilib[2,])]
  
  unfished <- calc.trajectory(lh, obs.cv = 0, init = rep(10000,length(lh$ages)), 
                              rec.dev = rep(1,times=years), F0=0, cr=NULL, 
                              years=years, hcr.type = "constF",obs.type="noerror",const.f.rate = 0, 
                              equilib=NULL,steepness=steepness)
  b.0 <- tail(unfished$biomass.oneplus.true,1)
  
  #Plot equilibrium yield curve
  x.vec <- fish
  y.vec <- equilib[3,]
  
  # x.vec <- equilib[2,]/b.inf
  # y.vec <- equilib[4,]
  
  #lines(equilib[2,]/b.inf,equilib[4,],type="l",xlim=c(0,1),main="Final year catch",xlab="Relative depletion",ylab="Catch (t)",col="blue")
  
  
  return(list("Fmsy"=f.msy,"Bmsy"=b.msy,"Binf"=b.inf,"B0" = b.0,"x.vec"=x.vec,"y.vec"=y.vec))
}


              # 
              # i=1
              # h = c(0.6,0.9)
              # lh.test$M <- plot.mat[i,1] 
              # xx <- getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=h[1])$x.vec
              # yy <- getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=h[1])$y.vec
              # plot(xx,yy,type='l',ylim=c(0,2e5),main="Final year catch",xlab="Relative depletion",ylab="Catch (t)")
              # 
              # for(i in 2:4){
              #   lh.test$M <- plot.mat[i,1]
              #   xx <- getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=h[1])$x.vec
              #   yy <- getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=h[1])$y.vec
              #   lines(xx,yy,col=rainbow(5)[i])
              #  }
              # 

