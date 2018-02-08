# Anchovy LH control file
# This file contains all the life history and fishery info to use in the MSE. Normally this would be a .ctl file but since I'm using R I'm just making it a source file for now.

# NOTE: THESE PARAMS ARE PIECED TOGETHER FROM SEVERAL REPORTS/PAPERS - THEY MAY BE IMPROVED IF I CAN GET SOUTH AFRICA ANCHOVY PARAMS OR ANOTHER STOCK W A GOOD ASSESSMENT

# Anchovy LH parameters

ages.test <- 0:6 # Mais CalCOFI report
nages.test <- length(ages.test)

#selectivity.test <- cbind(age=ages.test,selectivity = c(0,0.6,1.000,1.00,1.00,1,1))  # Fishery selectivity-- this is also made up except the part where age 0 
#aren't in survey or fishery-- that's a fact! 
selectivity.test <- cbind(age=ages.test,selectivity = c(0.1,0.6,1,1,1,1,1))
lh.test <- list(M = 1.05,   #mean between CA value (M = 1.2, from Jacobson et al. 2001, which cites Jacobson et al. (1994), Jacobson et al. (1995), and Methot (1989)), and the South Africa value (M = 0.9, Cunningham & Butterworth 2007) mean(c(1.2,0.9)
                selectivity = selectivity.test,
                ages = ages.test,
                l.at.age = c(8.0,10.8,14.0,16.3,17.8,18.9,19.6), # lengths in cm - from Huppert  et al. 1980
                w.at.age = c(0.0051,0.01278,0.02785,0.04395,0.05724,0.06852,0.07642)*0.001, #weights at age are in kg - from length-weight conversion in Huppert et al. -- so multiply by 0.001 to get mt
                maturity = c(0,0.4,0.85,0.99,1,1,1), # This is made up, based on anecdotal info from Alec (see email, end of January 2017).
                R0=1e9) # R0 is flexible, and can be changed. This one is based on the unfished behavior of the population...


#  Code for testing FMSY and other ref points -----------------------------
# Ideally the fishery parameters (like selectivity) could be adjusted until FMSY is similar to the estimate in the stock assessment. Howevever, for these stocks they aren't manage using FMSY, so I just looked for clues in the assessment about FMSY or proxy.
# selectivity.test <- cbind(ages.test, c(0.1,0.8,0.9,1,1,1,1))    
# 
        # setwd("/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Figures")
        # pdf("Anchovy_Sel_Mat_MSY.pdf",width = 8,height = 8,useDingbats = FALSE)
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
        # plot(equilib$x.vec,equilib$y.vec,type='l',xlab="F",ylab="Surplus Production",xlim=c(0,5))
        # text(x = equilib$Fmsy,y=250,"h = 0.6")
        # points(x=equilib$Fmsy,y=max(equilib$y.vec), pch=19)
        # points(x=equilib$Fmsy*0.5,y=equilib$y.vec[which.min(abs(equilib$x.vec - 0.5*equilib$Fmsy))], pch=19,col="blue")
        # 
        # equilib = getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=0.9)
        # plot(equilib$x.vec,equilib$y.vec,type='l',xlab="F",ylab="Surplus Production",xlim=c(0,5))
        # text(x = equilib$Fmsy,y=250,"h = 0.9")
        # points(x=equilib$Fmsy,y=max(equilib$y.vec), pch=19)
        # points(x=equilib$Fmsy*0.5,y=equilib$y.vec[which.min(abs(equilib$x.vec - 0.5*equilib$Fmsy))], pch=19,col="blue")
        # dev.off()