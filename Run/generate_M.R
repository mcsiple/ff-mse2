rw.M <- function(Mbar, rho.m, sigma.m, n){
  M <- rep(Mbar,times=n)
  nu <- vector(length=length(M))
  nu[1] <- 0.5
  for(y in 2:n){
    tau.y <- rnorm(1,0,sigma.m)
    nu[y] <- rho.m * nu[y-1] + sqrt(1-rho.m^2) * tau.y #anuual deviations in M
    M[y] <- Mbar * exp(nu[y]- (sigma.m^2 /2))
  }
return(M)
}

regime.M <- function(Mbar, cutoff.yr, n){
  # n is the length of the total simulation time series (not just the years used to calculate metrics)
  # likewise, cutoff.yr is marked in the entire simulation (so if you want it to be halfway thru the measured years, have to calc accordingly)
   M <- rep(Mbar, times=n)
   M[cutoff.yr:n] <- Mbar*1.25 
   return(M)
}

# rw.M(Mbar=0.8,rho.m = 0.5,sigma.m = 0.2,n=50)
# regime.M(Mbar=0.8,cutoff.yr = 25, n=50)
