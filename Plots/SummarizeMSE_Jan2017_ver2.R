# SummariseMSE.R

library(gtools)
library(RcppRoll)
library(plyr); library(dplyr)
library(scales)
library(gridExtra)
library(extrafont)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
source("/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Code/ff-mse2/Plots/Megsieggradar.R")
source("/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Code/ff-mse2/Plots/SummaryFxns.R")
Type = "Menhaden" #FF type to summarize



# Set path to wherever the simulation results are, load them into a giant dataframe
#path <- paste("/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Results/",Type,"/",sep="")
path <- "/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Results/Menhaden_SavedOutputs/Tau_06"
  files <- list.files(path=path)
  rm <- grep(files,pattern = ".txt") # Don't load the text summary
  files <- files[-rm]
  files <- mixedsort(files) # IMPORTANT: need this line to order in same order as scenario table!
  setwd(path)
  results <- sapply(files, function(x) mget(load(x)),simplify=TRUE) # This is a giant list of all the results - ONLY RDATA FILES and TXT FILES should be in this dir, otherwise you'll get an error

nscenarios <- length(results)
scen.table <- read.table("Scenario_Table.txt")
nsims <- nrow(results[[1]]$biomass) # just count nrows to know how many sims there are
years.test <- ncol(results[[1]]$biomass) # just count cols to know how many years there are
nyrs.to.use <- 100 # How many years you want to use to calculate all your metrics - There are no big differences btwn 50 and 100 yrs
calc.ind <- tail(1:years.test, nyrs.to.use) # Which years to calculate median depletion over (length = nyrs.to.use)

# Add performance measure columns to table
performance.measures <- c("LTmeancatch","LTnonzeromeancatch","SDcatch","n.5yrclose","n.10yrclose","nyrs0catch","meanbiomass","good4preds","SDbiomass","p.bad4preds","meanDepl")
# Still haven't added : prob(catch falls below a threshold bc what should the threshold be?)
                #   mean interannual change in catches
                #   min closure length - bc it's a tough coding thing (but do it later)
                #   Others (check notes!)

scen.table[,performance.measures] <- NA
raw.table <- scen.table # This will contain raw info about performance measures, when scen.table includes things that are scaled so that higher numbers = good


# ------------------------------------------------------------------------
# Summarize everything in one giant "outputs" table - this is ugly, sorry
# ------------------------------------------------------------------------

for (s in 1:nscenarios){
  #**N** indicate metrics for which higher values mean worse performance (like SD(catch)) - these metrics are in scen.table as 1/x
  result.to.use <- results[[s]]
  scen.table[s,performance.measures[1]] <- raw.table[s,performance.measures[1]] <- median(rowMeans(result.to.use$total.catch[,calc.ind])) #calculate mean B over years to use in the index - the final number is the median (across all simulations) mean B
  nonzero.catch <- result.to.use$total.catch[,calc.ind]
  nonzero.catch <- ifelse(nonzero.catch<0.1,NA,nonzero.catch)
  mnz.catches <- rowMeans(nonzero.catch,na.rm=TRUE)
  #mnz.catches[ is.na(mnz.catches) ] <- NA
  scen.table[s,performance.measures[2]] <- raw.table[s,performance.measures[2]] <- median(mnz.catches,na.rm = TRUE) #"LTnonzeromeancatch"
  
  scen.table[s,performance.measures[3]] <- 1 / median(apply(X = result.to.use$total.catch[,calc.ind],FUN = sd,MARGIN = 1)) #"SDcatch" **N**
  raw.table[s,performance.measures[3]] <- median(apply(X = result.to.use$total.catch[,calc.ind],FUN = sd,MARGIN = 1)) 
  
  mean.num5yrclosures <- mean(n.multiyr.closures(result.to.use$total.catch[,calc.ind],threshold=0)$count5, na.rm=TRUE)
  scen.table[s,performance.measures[4]] <- ifelse(mean.num5yrclosures ==0, 1,
                                                  1 / mean.num5yrclosures) # "n.5yrclose" **N**
  raw.table[s,performance.measures[4]] <- mean.num5yrclosures # Mean number of 5-yr closures
  
  mean.num10yrclosures <- mean(n.multiyr.closures(result.to.use$total.catch[,calc.ind],threshold=0)$count10, na.rm=TRUE)
  scen.table[s,performance.measures[5]] <- ifelse(mean.num10yrclosures==0, 1,
                                                  1 / mean.num10yrclosures) # "n.10yrclose" **N**
  raw.table[s,performance.measures[5]] <- mean.num10yrclosures
  
  scen.table[s,performance.measures[6]] <- ifelse(median(apply(X = result.to.use$total.catch[,calc.ind],FUN = nzeroes,MARGIN = 1)) ==0, 1,
                                                  1 / median(apply(X = result.to.use$total.catch[,calc.ind],FUN = nzeroes,MARGIN = 1)) ) # "nyrs0catch" **N**
  raw.table[s,performance.measures[6]] <- median(apply(X = result.to.use$total.catch[,calc.ind],FUN = nzeroes,MARGIN = 1))
  
  scen.table[s,performance.measures[7]] <- 
    raw.table[s,performance.measures[7]] <- median(rowMeans(result.to.use$total.true.biomass[,calc.ind])) # "meanbiomass"
  
  scen.table[s,performance.measures[8]] <- 
    raw.table[s,performance.measures[8]] <- median(apply(X = result.to.use$total.true.biomass[,calc.ind],FUN = good4pred,MARGIN = 1, 
                                                         F0.x = result.to.use$no.fishing.tb[,calc.ind])) # "good4preds" - this is calculated from TOTAL biomass (including age 0)
  
  scen.table[s,performance.measures[9]] <- 1 / median(apply(X = result.to.use$total.true.biomass[,calc.ind],FUN = sd,MARGIN = 1)) # "SDbiomass" **N**
  raw.table[s,performance.measures[9]] <- median(apply(X = result.to.use$total.true.biomass[,calc.ind],FUN = sd,MARGIN = 1))  #Actual raw SD of Biomass
  
  yrs.bad <- apply(X = result.to.use$total.true.biomass[,calc.ind],FUN = bad4pred,MARGIN = 1, F0.x = result.to.use$no.fishing.tb[,calc.ind] ) # length of vector is nsims 
  scen.table[s,performance.measures[10]] <- ifelse(length(which(yrs.bad>0))/nsims == 0, 0,
                                               1 / length(which(yrs.bad>0))/nsims ) #Prob that biomass ever dipped below 10% of LT mean "p.bad4preds" **N**
  raw.table[s,performance.measures[10]] <- length(which(yrs.bad>0))/nsims
  scen.table[s,performance.measures[11]] <- 
    raw.table[s,performance.measures[11]] <- median(rowMeans(result.to.use$depl[,calc.ind])) 
}

write.csv(raw.table, file=paste(Type,"_outputs.csv",sep=""))


############################################################################
# PLOTS AND METRICS TO SHOW OUTPUTS ----------------------------
############################################################################

# Fxns for summarizing and plotting ---------------------------------------

# Colour palette for time series plots - some of these are from iWantHue and some are ColorBrewer
#palette <- brewer_pal(type="qual",palette=2)
#palette <- c("#d94313","#3097ff","#f5bd4e","#e259db","#009a3b","#da0b96","#38e096","#ff4471","#007733","#ff90f5","#588400","#feaedc","#a1d665","#42c7ff","#6f5500","#01b1be") 
palette <- brewer.pal(6,"Spectral")
show_col(palette)
hcr.colors <- palette[c(6,5,3,1,2)]
show_col(hcr.colors) # C1 (Oc), C2 (Len), constF, stability-favoring, trend-based (this is the order of the colors)

all.summaries <- lapply(results,FUN = summ.tab)
all.summaries <- do.call(rbind.data.frame, all.summaries)
all.summaries$scenario <- rep(1:20,each=length(performance.measures))

# Match the scenarios to type of error, etc.
all.summaries <- merge(all.summaries,scen.table[,1:6],by="scenario")
all.summaries <- mutate(all.summaries, obs.error.type = recode(obs.error.type, 
                                                               'Tim'='Delayed change detection',
                                                               'AC' = "Autocorrelated"),
                                        HCR = recode(HCR, 'cfp' = 'Stability-favoring',
                                                    'constF' = 'Constant F',
                                                    'lenfest' = 'C2',
                                                    'oceana' = 'C1',
                                                    'trend' = "Trend-based"))
all.summaries$HCR <- factor(all.summaries$HCR, levels = c("C1","C2","Constant F","Stability-favoring","Trend-based")) # Reorder factors so they plot in alphabetical order, the way they were intended to be!
# Change labels of things in the table! --------------------------
scen.table <- mutate(scen.table, obs.error.type = recode(obs.error.type, 
                                                        'Tim'='Delayed change detection',
                                                        'AC' = "Autocorrelated"),
                                           HCR = recode(HCR, 'cfp' = 'Stability-favoring',
                                                        'constF' = 'Constant F',
                                                        'lenfest' = 'C2',
                                                        'oceana' = 'C1',
                                                        'trend' = "Trend-based"))


######################################################################
###### MAKE A PDF WITH ALL THE OUTPUT FIGURES! #######################
######################################################################

pdf(paste(Type,"_May16.pdf",sep=""),width = 10,height = 9,onefile = TRUE)
# Put control rules in order so they plot right
scen.table$HCR <- factor(scen.table$HCR, levels = c("C1","C2","Constant F","Stability-favoring","Trend-based"))
# Compare each of the CRs together? It would be like pairs()



# Kite plots showing tradeoffs --------------------------------------------
mat <- matrix(1:4,nrow=2,byrow = TRUE)
plotnames <- list()
steepnesses <- unique(scen.table$h)
obs.error.types <- unique(scen.table$obs.error.type)
nice.pms <- data.frame(original = colnames(scen.table[-(1:6)]),
                       polished = c("LT mean catch","LT mean nonzero catch",
                                    "SD(Catch)","Number of \n 5-yr closures",
                                    "Number of \n 10-yr closures","Number of yrs \n w/ zero catch",
                                    "LT mean biomass","Number of yrs \n above pred threshold",
                                    "SD(Biomass)","Prob. biomass falls \n below pred threshold",
                                    "Mean depletion"))

for(steep in 1:2){
  for(obs in 1:2){
    tab <- subset(scen.table, obs.error.type == obs.error.types[obs] & h == steepnesses[steep])
    tab.metrics <- tab[,7:17]
    props <- tab.metrics
    maxes <- apply(X = tab.metrics,MARGIN = 2,FUN = max)
    for(i in 1:nrow(props)){
      props[i,] <- tab.metrics[i,] / maxes
    }
    final.tab <- cbind(tab[,3],props)
    colnames(final.tab)[1] <- "group"
    test.nas <- apply(X = final.tab,FUN = anyNA,MARGIN = 2)
    na.metrics <- names(which(test.nas))
    axis.labels <- nice.pms[,'polished']
    
    if(length(na.metrics>0)){
      print(paste("The following performance metrics had NAs and was removed from the figure: ",na.metrics))
      final.tab <- final.tab[,-which(test.nas)]
      rm.metrics <- which(nice.pms$original == na.metrics)
      axis.labels <- nice.pms[-rm.metrics,'polished']
    }
    
    legend.presence <- ifelse(mat[steep,obs] != 1,FALSE,TRUE)
    plotnames[[mat[steep,obs]]] <- ggradar(final.tab,font.radar = "Helvetica",grid.label.size=3,axis.label.size=2.5,
                                           legend.text.size = 2.5,
                                           axis.labels = axis.labels,
                                           plot.legend=legend.presence,palette.vec = hcr.colors)
  }}

grid.arrange(plotnames[[1]],plotnames[[2]],plotnames[[3]], plotnames[[4]])


# Zeh plots (e.g. Punt 2015) -----------------------------------------------

ggplot(all.summaries,aes(x=HCR,y=med,colour=HCR)) +
  geom_point(size=3,shape=20) +
  geom_pointrange(aes(ymin=loCI,ymax=hiCI)) +
  facet_grid(PM~h+obs.error.type, scales = "free_y") +
  scale_colour_manual(values = hcr.colors) +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5)) +
  ylab("Median (+/- 90% quantile)")

# Tradeoff figures ---------------------------------------------------
tradeoff.plot <- subset(scen.table,h==0.9)
# Scale all performance measures to max
for(i in 7:ncol(tradeoff.plot)){
  tradeoff.plot[,i] <- tradeoff.plot[,i] / max(tradeoff.plot[,i],na.rm=T) *100
}
melt.tradeoff <- melt(tradeoff.plot,id.vars=c("h","recruit.sd","recruit.rho","obs.error.type","HCR","scenario"))
ggplot(melt.tradeoff,aes(x=HCR,y=value,colour=HCR,shape=obs.error.type,label=scenario,alpha=obs.error.type)) +
  xlab("Control rule") + ylab("% of best performance, h = 0.9") +
  scale_colour_manual(values = hcr.colors) +
  scale_alpha_manual(values = c(0.6,1)) +
  geom_point(size=5)  + theme_bw(base_size = 18) + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5)) +
  facet_wrap(~variable)


tradeoff.plot <- subset(scen.table,h==0.6)
# Scale all performance measures to max
for(i in 7:ncol(tradeoff.plot)){
  tradeoff.plot[,i] <- tradeoff.plot[,i] / max(tradeoff.plot[,i],na.rm=T) *100
}
melt.tradeoff <- melt(tradeoff.plot,id.vars=c("h","recruit.sd","recruit.rho","obs.error.type","HCR","scenario"))
ggplot(melt.tradeoff,aes(x=HCR,y=value,colour=HCR,shape=obs.error.type,label=scenario,alpha=obs.error.type)) +
  xlab("Control rule") + ylab("% of best performance, h = 0.6") +
  scale_colour_manual(values = hcr.colors) +
  scale_alpha_manual(values = c(0.6,1)) +
  geom_point(size=5)  + theme_bw(base_size = 18) + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5)) +
  facet_wrap(~variable)

# Plot comparing biomass, catch, etc for one run each scenario ------------
attributes <- c("Biomass","Catches","Recruitment","Depletion (B/B0)")
att.ind <- c(1,2,4,5)

# Plot first run for each scenario
# Color order of these is funky because order of results list() is not the same as the order of  scen.table w/ the summary
for (scenario.index in 1:4){
  par(mfrow=c(2,2),mar=c(5,4,3,2)+0.1)
  for(i in 1:length(att.ind)){
    which.metric <- att.ind[i]
    plot(results[[scenario.index]][[which.metric]][1,calc.ind],col=hcr.colors[4],
         ylab = attributes[i], type="l", lwd=2,xlab="Year",
         ylim = range(c(results[[scenario.index]][[which.metric]][1,calc.ind],
                        results[[scenario.index+4]][[which.metric]][1,calc.ind],
                        results[[scenario.index+8]][[which.metric]][1,calc.ind],
                        results[[scenario.index+12]][[which.metric]][1,calc.ind],
                        results[[scenario.index+16]][[which.metric]][1,calc.ind])))
                        # The above is so messy but it's just to get the correct range for all the lines. 
    lines(results[[scenario.index+4]][[which.metric]][1,calc.ind],lwd=2,col=hcr.colors[3]) 
    lines(results[[scenario.index+8]][[which.metric]][1,calc.ind],lwd=2,col=hcr.colors[1]) 
    lines(results[[scenario.index+12]][[which.metric]][1,calc.ind],lwd=2,col=hcr.colors[2])
    lines(results[[scenario.index+16]][[which.metric]][1,calc.ind],lwd=2,col=hcr.colors[5])
  }
  
  mtext(paste("h=",scen.table[scenario.index,'h'],"-","ObsErrorType =",scen.table[scenario.index,'obs.error.type'],sep= " "),outer=TRUE,cex=1.2,side=3,line=-3)
  legend("topright",col = hcr.colors[1:5],lwd=rep(2,times=5),legend = sort(unique(scen.table$HCR)))
}

# Plot observed vs. true one-plus biomass, to show differences in true vs. observed B
par(mfrow=c(4,5))
for (scenario.index in 1:4){
  #par(mfrow=c(2,2),mar=c(5,4,3,2)+0.1)
    ts <- results[[scenario.index]][["biomass"]][1,calc.ind]
    plot(ts,col=hcr.colors[4],
         ylab = "Biomass", type="l", lwd=2,xlab="Year")
        lines(results[[scenario.index]][["obs.biomass"]][1,calc.ind],col = hcr.colors[4],lty=2)
        text(x = 40,y=max(ts)*0.95,labels = paste(scen.table[scenario.index,'obs.error.type'],sep= " "))
        text(x = 40,y=max(ts)*0.85,labels = paste(scen.table[scenario.index,'h'],sep= " "))
    # The above is so messy but it's just to get the correct range for all the lines. 
    ts <- results[[scenario.index+4]][["biomass"]][1,calc.ind]
    plot(ts,col=hcr.colors[3],
         ylab = "Biomass", type="l", lwd=2,xlab="Year")
    lines(results[[scenario.index+4]][["obs.biomass"]][1,calc.ind],col=hcr.colors[3],lty=2) 
    text(x = 40,y=max(ts)*0.95,labels = paste(scen.table[scenario.index,'obs.error.type'],sep= " "))
    
    ts <- results[[scenario.index+8]][["biomass"]][1,calc.ind]
    plot(ts,col=hcr.colors[1],
         ylab = "Biomass", type="l", lwd=2,xlab="Year")
    lines(results[[scenario.index+8]][["obs.biomass"]][1,calc.ind],col=hcr.colors[1],lty=2) 
    text(x = 40,y=max(ts)*0.95,labels = paste(scen.table[scenario.index,'obs.error.type'],sep= " "))
    
    ts <- results[[scenario.index+12]][["biomass"]][1,calc.ind]
    plot(ts,col=hcr.colors[2],
         ylab = "Biomass", type="l", lwd=2,xlab="Year")
    lines(results[[scenario.index+12]][["obs.biomass"]][1,calc.ind],col=hcr.colors[2],lty=2) 
    text(x = 40,y=max(ts)*0.95,labels = paste(scen.table[scenario.index,'obs.error.type'],sep= " "))
    
    ts <- results[[scenario.index+16]][["biomass"]][1,calc.ind]
    plot(ts,col=hcr.colors[5],
         ylab = "Biomass", type="l", lwd=2,xlab="Year")
    lines(results[[scenario.index+16]][["obs.biomass"]][1,calc.ind],col=hcr.colors[5],lty=2)
    text(x = 40,y=max(ts)*0.95,labels = paste(scen.table[scenario.index,'obs.error.type'],sep= " "))
    
  legend("topright",col = c("black","black"),lwd=c(2,2),lty=c(1,2),legend = c("true B1+","obs B1+"))
}

# Medians and 95% intervals for biomass, catch, etc ------------------


for (scenario in 1:length(results)){
  par(mfrow=c(2,2),mar=c(5,4,3,2)+0.1)
  for(i in 1:length(attributes)){
    plotintervals(results[[scenario]][[att.ind[i]]],ylab = attributes[i])
  }
  mtext(paste(scen.table[scenario,'HCR']), outer = TRUE, cex = 1.5,side = 3,line = -1.8)
  mtext(paste("h=",scen.table[scenario,'h'],"-","ObsErrorType =",scen.table[scenario,'obs.error.type'],sep= " "),outer=TRUE,cex=1.2,side=3,line=-3)
}

# Sample time series (just plot a few draws from each scenario) -----------

linecol <- add.alpha ("grey",alpha=0.4)
maxsim = 5 #Max number of example time series to plot

for (scenario in 1:length(results)){
  par(mfrow=c(2,2),mar=c(5,4,3,2)+0.1)
  for(i in 1:length(attributes)){
    plot(results[[scenario]][[att.ind[i]]][1,calc.ind],
         ylab = attributes[i], type="l", lwd=2,col=linecol,xlab="Year",
         ylim = range(results[[scenario]][[att.ind[i]]][1:maxsim,calc.ind]))
    for (sim in 2:maxsim){
      lines(results[[scenario]][[att.ind[i]]][sim,calc.ind],
            lwd=2,col=linecol)
    }
  }
  mtext(paste(scen.table[scenario,'HCR']), outer = TRUE, cex = 1.5,side = 3,line = -1.8)
  mtext(paste("h=",scen.table[scenario,'h'],"-","ObsErrorType =",scen.table[scenario,'obs.error.type'],sep= " "),outer=TRUE,cex=1.2,side=3,line=-3)
}


dev.off()