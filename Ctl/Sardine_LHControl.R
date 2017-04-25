# Sardine LH control file
# This file contains all the life history and fishery info to use in the MSE. Normally this would be a .ctl file but since I'm using R I'm just making it a source file for now.


# Sardine LH parameters

ages.test <- 0:15
nages.test <- length(ages.test)
selectivity.test <- cbind(age=ages.test,selectivity = c(0.263,1.00,1.000,0.669,0.471,0.390,0.358,0.345,0.339,0.335,0.333,0.332,0.332,0.331,0.331,0.331))  #fishery selectivity from sardine assessment; see Table App.A.1. in Hurtado-Ferro & Punt 2014 
#Appendix A.1 caption: "Fleet-averaged selectivity (computed using the output of model X6e of Hill et al. [2012]). Results are shown for 2011, 2007-2011, and 2002-2011."
lh.test <- list(M = 0.4,   #from sardine assessment (Hurtado-Ferro & Punt 2014)
                selectivity = selectivity.test,
                ages = ages.test,
                l.at.age = seq(4,28,length.out=nages.test), # Lengths are in mm...? Doesn't actually matter because lengths aren't used in the most basic model
                w.at.age = seq(0.014,0.180,length.out=nages.test)*0.001, #weights at age are in kg; multiplying by 0.001 gives mt
                maturity = c(0,0.4,0.85,0.99,rep(1,times=12)), # Roughly based on maturity ogives from Silva et al. 2006 ICES JMS.
                R0=1e9) # R0 is flexible, and can be changed. I kind of made it up.

