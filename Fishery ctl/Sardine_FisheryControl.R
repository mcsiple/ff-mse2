#Sardine fishery control file

# Years, fishing, initial age dist

tot <- sum(c(6.221,4.170,2.795,1.873,1.256,0.842,0.564,0.378,0.253,0.169,rep(0.0575,times=6)))
init.prop <- c(6.221,4.170,2.795,1.873,1.256,0.842,0.564,0.378,0.253,0.169,rep(0.0575,times=6))/tot #From Hill et al. 2012 (http://www.pcouncil.org/wp-content/uploads/MAIN_DOC_G3b_ASSMNT_RPT2_WEB_ONLY_NOV2012BB.pdf). The plus group (10+) is split up into 4, approximately. The initial number at age should be based on K and the age distribution.
init.B <- init.prop*200000 # Biomass at age -- 110.54 is the B0 when R0 = 1000
init.test <- init.B / lh.test$w.at.age # These are all the numbers at age (spawning and non-spawning biomass)
F0.test <- 0.3 # This is arbitrary
cr.test <- list(Blim = 100, Btarget = 200, Fmax = 0.6) 