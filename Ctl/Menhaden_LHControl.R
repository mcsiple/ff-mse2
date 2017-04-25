# Menhaden_LHControl
# Life history params for the menhaden model, to be used with menhaden-like recruitment patterns

ages.test <- 0:6
nages.test <- length(ages.test)
selectivity.test <- cbind(age=ages.test,selectivity = c(0.004,0.143,0.994,0.84,0.191,0.024,0))  # "domed" commercial selectivity from reduction fishery. Used Butterworth 2012 as an example
lh.test <- list(M = 0.5,   #from menhaden assessment (http://www.asmfc.org/uploads/file/55089931S40_AtlMenhadenSAR_CombinedFINAL_1.15.2015-reduced.pdf) - mortality is averaged over all ages, from all years (Boudreau & Dickie)
                selectivity = selectivity.test,
                ages = ages.test,
                l.at.age = seq(4,28,length.out=nages.test),
                w.at.age = c(0.0569, 0.1281,0.2317,0.3285,0.3711,0.5371,0.4481)*0.001, #weights at age are in kg  I think, multiply by 0.001 to get mt
                maturity = c(0,0.13,0.53,0.83,0.98,1,1), 
                R0=1e9) # R0 is flexible, and can be changed. I made it up.

