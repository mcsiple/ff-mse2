#Menhaden fishery control file (if possible, fix initial nums at age - can't find in stock assessment but get similar long-term results w different values)
# Menhaden parameters come from 
tot <- sum(c(6.221,4.170,2.795,1.873,1.256,0.842,0.564))
init.prop <- c(6.221,4.170,2.795,1.873,1.256,0.842,0.564)/tot # Initial proportion at age
init.B <- init.prop*200000 # Biomass at age -- 110.54 is the B0 when R0 = 1000
init.test <- init.B / lh.test$w.at.age # These are all the numbers at age (spawning and non-spawning biomass)
F0.test <- 0.3 # Random F0 for first year
# cr.test <- list(Blim = 100, Btarget = 200, Fmax = 0.6) 
