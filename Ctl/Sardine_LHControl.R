# Sardine life history control file
# This file contains all the life history and fishery info to use in the MSE. Normally this would be a .ctl file but since I'm using R I'm just making it a source file for now.


# Sardine LH parameters
ages.test <- 0:15
nages.test <- length(ages.test)
selectivity.test <- cbind(age=ages.test,selectivity = c(0.263,1.00,1.000,0.669,0.471,0.390,0.358,0.345,0.339,0.335,0.333,0.332,0.332,0.331,0.331,0.331))  #fishery selectivity from sardine assessment; Table App.A.1. in Hurtado-Ferro & Punt 2014 
#Appendix A.1 caption: "Fleet-averaged selectivity (computed using the output of model X6e of Hill et al. [2012]). Results are shown for 2011, 2007-2011, and 2002-2011."
lh.test <- list(M = 0.4,   #from sardine assessment (Hurtado-Ferro & Punt 2014 which I think uses values from Hill et al. 2014)
                selectivity = selectivity.test,
                ages = ages.test,
                l.at.age = seq(4,28,length.out=nages.test), # (mm) These lengths are not real; not in final model
                w.at.age = seq(0.014,0.180,length.out=nages.test)*0.001, # Weights at age are in kg; multiplying by 0.001 gives mt
                #maturity = c(0,0.4,0.85,0.99,rep(1,times=12)), # Roughly based on maturity ogives from Silva et al. 2006 ICES JMS
                maturity = c(0,0.277,0.682,0.906,0.973,0.988,0.996,rep(1,times=9)),
                # Age(y) Prop. mature, from Hill et al. 2014 Figure 5b
                # 0.993	0.277
                # 2.005	0.682
                # 3.017	0.906
                # 4.003	0.973
                # 5.001	0.988
                # 6.048	0.996
                R0=5e9) # R0 is flexible, and can be changed. It was estimated in Hill et al. to be approx. 4.8 billion fish (https://swfsc.noaa.gov/uploadedFiles/Events/Meetings/Meeting_2014/H1b_2014_FULL_Electric_PacificSardine_StockAssmnt_APR2014BB.pdf)


#  Plot mat and sel to make sure FMSY values aren’t weird -----------------
        # setwd("/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Figures")
        # pdf("Sardine_Sel_Mat_MSY.pdf",width = 8,height = 8,useDingbats = FALSE)
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
        # plot(equilib$x.vec,equilib$y.vec,type='l',xlab="F",ylab="Surplus Production",xlim=c(0,4))
        # text(x = equilib$Fmsy+2,y=4000,"h = 0.6")
        # points(x=equilib$Fmsy,y=max(equilib$y.vec), pch=19)
        # points(x=equilib$Fmsy*0.5,y=equilib$y.vec[which.min(abs(equilib$x.vec - 0.5*equilib$Fmsy))], pch=19,col="blue")
        # 
        # equilib = getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=0.9)
        # plot(equilib$x.vec,equilib$y.vec,type='l',xlab="F",ylab="Surplus Production",xlim=c(0,4))
        # text(x = equilib$Fmsy,y=4000,"h = 0.9")
        # points(x=equilib$Fmsy,y=max(equilib$y.vec), pch=19)
        # points(x=equilib$Fmsy*0.5,y=equilib$y.vec[which.min(abs(equilib$x.vec - 0.5*equilib$Fmsy))], pch=19,col="blue")
        # dev.off()
        
