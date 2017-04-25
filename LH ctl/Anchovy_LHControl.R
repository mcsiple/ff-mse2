# Anchovy LH control file
# This file contains all the life history and fishery info to use in the MSE. Normally this would be a .ctl file but since I'm using R I'm just making it a source file for now.

# NOTE: THESE PARAMS ARE PIECED TOGETHER FROM SEVERAL REPORTS/PAPERS - THEY MAY BE IMPROVED IF I CAN GET SOUTH AFRICA ANCHOVY PARAMS OR ANOTHER STOCK W A GOOD ASSESSMENT
# Anchovy LH parameters

ages.test <- 0:6 # Mais CalCOFI report
nages.test <- length(ages.test)

selectivity.test <- cbind(age=ages.test,selectivity = c(0,0.6,1.000,1.00,1.00,1,1))  # Fishery selectivity-- this is also made up except the part where age 0 aren't in survey or fishery-- that's a fact! 
lh.test <- list(M = 1.2,   #from Jacobson et al. 2001, which cites Jacobson et al. (1994), Jacobson et al. (1995), and Methot (1989)
                selectivity = selectivity.test,
                ages = ages.test,
                l.at.age = c(8.0,10.8,14.0,16.3,17.8,18.9,19.6), # lengths in cm - from Huppert  et al. 1980
                w.at.age = c(0.0051,0.01278,0.02785,0.04395,0.05724,0.06852,0.07642)*0.001, #weights at age are in kg - from length-weight conversion in Huppert et al. -- so multiply by 0.001 to get mt
                maturity = c(0,0.4,0.85,0.99,1,1,1), # This is made up, based on anecdotal info from Alec (see email, end of January 2017).
                R0=1e9) # R0 is flexible, and can be changed. I made it up here, based on the unfished behavior of the population...

