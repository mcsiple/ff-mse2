#Menhaden fishery control file (NEED TO FIX BECAUSE I DON'T ACTUALLY HAVE ALL THE RIGHT VALUES IN HERE!)
tot <- sum(c(6.221,4.170,2.795,1.873,1.256,0.842,0.564))
init.prop <- c(6.221,4.170,2.795,1.873,1.256,0.842,0.564)/tot 
init.B <- init.prop*200000 # Biomass at age -- 110.54 is the B0 when R0 = 1000
init.test <- init.B / lh.test$w.at.age # These are all the numbers at age (spawning and non-spawning biomass)
F0.test <- 0.3 # This is arbitrary
cr.test <- list(Blim = 100, Btarget = 200, Fmax = 0.6) 