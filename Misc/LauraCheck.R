
# Step 1: Get results -----------------------------------------------------


library(gtools)
library(RcppRoll)
library(plyr); 
library(dplyr)
library(scales)
library(gridExtra)
library(extrafont)
library(RColorBrewer)
library(reshape2)
library(ggplot2)

source("/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Code/ff-mse2/Plots/Megsieggradar.R")
source("/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Code/ff-mse2/Plots/SummaryFxns.R")
source("/Users/mcsiple/Dropbox/ChapterX-synthesis/Theme_Black.R")
Type = "Anchovy" #FF type to summarize
Date <- "2020-11-04"

# Set path to where the simulation results are:
path <- paste("/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Results/",Type,Date,"/",sep="")

setwd(path)    # ONLY RDATA FILES and TXT FILES should be in this dir

# Read all files into giant list
files <- list.files(path=path)
rm <- grep(files,pattern = ".txt") # Don't load text summary
files <- files[-rm]
files <- files[grep(files, pattern = ".RData")] # load rdata files
files <- mixedsort(files) # order in same order as scenario table
results <- sapply(files, function(x) mget(load(x)),simplify=TRUE) # giant list


# Step 2: Check the results and look at what they say ---------------------

str(results)
raw.table
results2 <- results[[2]] # h = 0.6
str(results2) 

sim = 4

plot(1:250,results2$biomass.total.true[sim,],type = 'l')
lines(1:250,results2$no.fishing.tb[sim,],col='red')
abline(h = 0.4*mean(results2$no.fishing.tb[sim,]),lty=2)

collapse.index(x = results2$biomass.total.true[sim,],F0.x = results2$no.fishing.tb[sim,])

