#Anchovy fishery control file

# Years, fishing, initial age dist
tot <- sum(c(6.221,4.170,2.795,1.873,1.256,0.842,0.564))
init.prop <- c(6.221,4.170,2.795,1.873,1.256,0.842,0.564)/tot # These proportions at age are the same as the first few age classes of sardine, just because I don't have any other info right now! The initial number at age should be based on K and the age distribution.
init.B <- init.prop*20000 # Multiplier should be total biomass in the first year!) - this value matches the biomass of CA sardine ~1960 (MacCall et al. 2016)
init.test <- init.B / lh.test$w.at.age # These are all the numbers at age (spawning and non-spawning biomass)
#init.test is initial population size... 
F0.test <- 0.3 # This is arbitrary
#cr.test <- list(Blim = 100, Btarget = 200, Fmax = 0.6) 

