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
            # setwd("/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Figures")
            # pdf(file = "Menhaden_Sel_Mat_MSY.pdf",width = 8,height = 8,useDingbats = FALSE)
            # par(mfrow=c(2,2))
            # plot(ages.test,selectivity.test[,2],type='l',col='red',xlab="Age",ylab = "Selectivity / Maturity") #main="Menhaden"
            # lines(ages.test,lh.test$maturity,col='blue')
            # 
            # ## Empty panel with legends
            # plot(0,type='n',axes=FALSE,ann=FALSE)
            # legend("top",lwd=c(1,1),col=c('red','blue'),legend = c("selectivity","maturity"))
            # legend("center",legend = c("Fmsy","0.5*Fmsy"),pch=c(19,19),col=c("black","blue"))
            # 
            # equilib = getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=0.6)
            # plot(equilib$x.vec,equilib$y.vec,type='l',xlab="F",ylab="Surplus Production")
            # text(x = equilib$Fmsy,y=4000,"h = 0.6")
            # points(x=equilib$Fmsy,y=max(equilib$y.vec), pch=19)
            # points(x=equilib$Fmsy*0.5,y=equilib$y.vec[which.min(abs(equilib$x.vec - 0.5*equilib$Fmsy))], pch=19,col="blue")
            # 
            # equilib = getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=0.9)
            # plot(equilib$x.vec,equilib$y.vec,type='l',xlab="F",ylab="Surplus Production")
            # text(x = equilib$Fmsy,y=4000,"h = 0.9")
            # points(x=equilib$Fmsy,y=max(equilib$y.vec), pch=19)
            # points(x=equilib$Fmsy*0.5,y=equilib$y.vec[which.min(abs(equilib$x.vec - 0.5*equilib$Fmsy))], pch=19,col="blue")
            # 
            # dev.off()

