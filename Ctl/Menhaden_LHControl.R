# Menhaden_LHControl
# Life history params for the menhaden model, to be used with menhaden-like recruitment patterns
# These are based on the Atlantic menhaden stock assessment
ages.test <- 0:6
nages.test <- length(ages.test)
selectivity.test <- cbind(age=ages.test,selectivity = c(0.004,0.143,0.994,0.84,0.191,0.024,0))  # "domed" commercial selectivity from reduction fishery. Used Butterworth 2012 as an example
lh.test <- list(type = "Menhaden",
                M = 0.5,   #from menhaden assessment (http://www.asmfc.org/uploads/file/55089931S40_AtlMenhadenSAR_CombinedFINAL_1.15.2015-reduced.pdf) - mortality is averaged over all ages, from all years (Boudreau & Dickie) - table 3.6.1
                selectivity = selectivity.test,
                ages = ages.test,
                l.at.age = seq(4,28,length.out=nages.test), #http://www.asmfc.org/uploads/file/55089931S40_AtlMenhadenSAR_CombinedFINAL_1.15.2015-reduced.pdf (Table 3.3.2)
                w.at.age = c(0.0569, 0.1281,0.2317,0.3285,0.3711,0.5371,0.4481)*0.001, #weights at age are in kg  I think, multiply by 0.001 to get mt
                maturity = c(0,0.13,0.53,0.83,0.98,1,1), # Table 3.4.1
                R0=1e9) # R0 is flexible, and can be changed. I made it up.

    #  Plot mat and sel to make sure FMSY values aren’t weird -----------------
    # plot(ages.test,selectivity.test[,2],type='l',col='red')
    # lines(ages.test,lh.test$maturity,col='blue')
    # legend("topleft",lwd=c(1,1),col=c('red','blue'),legend = c("selectivity","maturity"))
    # 
    # equilib = getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=steepness)
    # plot(equilib$x.vec,equilib$y.vec)
