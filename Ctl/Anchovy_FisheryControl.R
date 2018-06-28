#Anchovy fishery control file

# Years, fishing, initial age dist
nums <- c(6.221,
           4.170,
           2.795,
           1.873,
           1.256,
           0.842,
           0.564)
tot <- sum(nums)
init.prop <- nums/tot # These proportions at age are the same as the first few age classes of sardine,  because numbers at age aren't available for anchovy yet. Initial number at age should be based on K and the age distribution.
init.B <- init.prop*20000 # Multiplier = total biomass in the first year
init.test <- init.B / lh.test$w.at.age # Initial population size
F0.test <- 0.3 # FO for testing

#cr.test <- list(Blim = 100, Btarget = 200, Fmax = 0.6) 

