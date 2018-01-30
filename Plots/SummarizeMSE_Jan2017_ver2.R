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
Date <- "2017-10-05"


# Set path to wherever the simulation results are, load them into a giant dataframe
path <- paste("/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Results/",Type,Date,"/",sep="")
  files <- list.files(path=path)
  rm <- grep(files,pattern = ".txt") # Don't load the text table summary
  files <- files[-rm]
  files <- mixedsort(files) # IMPORTANT: need this line to order in same order as scenario table!
  setwd(path)
  results <- sapply(files, function(x) mget(load(x)),simplify=TRUE) # This is a giant list of all the results - ONLY RDATA FILES and TXT FILES should be in this dir, otherwise you'll get an error

nscenarios <- length(results)
scen.table <- read.table("Scenario_Table.txt")  #This empty table is generated when the simulations are run - then you fill it in after everything runs
nsims <- nrow(results[[1]]$biomass.oneplus.true) # just count nrows to know how many sims there are
years.test <- ncol(results[[1]]$biomass.oneplus.true) # just count cols to know how many years there are
nyrs.to.use <- 100 # How many years you want to use to calculate all your metrics - There are no big differences btwn 50 and 100 yrs
calc.ind <- tail(1:years.test, nyrs.to.use) # Which years to calculate PMs over (length = nyrs.to.use)

# Add performance measure columns to table
performance.measures <- c("LTmeancatch","LTnonzeromeancatch","SDcatch","n.5yrclose","n.10yrclose","nyrs0catch","meanbiomass","good4preds","SDbiomass","very.bad4preds","meanDepl","overallMaxCollapseLength","overallMaxBonanzaLength","BonanzaLength","CollapseLength","Prob.Collapse","Collapse.Severity","CV.Catch","Sustained.collapse")
pm.type <- c(rep("Fishery",times=6),rep("Ecosystem",times=11)) # for distinguishing types of PMs (mostly for plotting...)

#overall.max.coll.len,overall.max.bon.len,bon.length,coll.length
# Still haven't added : prob(catch falls below a threshold bc what should the threshold be?)
                #   mean interannual change in catches
                #   min closure length 
                #   Others (check notes!)

scen.table[,performance.measures] <- NA
raw.table <- scen.table # This will contain raw info about performance measures, whereas scen.table includes things that are scaled so that higher numbers = good


# ------------------------------------------------------------------------
# Summarize everything in one giant "outputs" table - this is ugly, sorry
# ------------------------------------------------------------------------

# Fxns for summarizing and plotting ---------------------------------------

# Colour palette for time series plots - some of these are from iWantHue and some are ColorBrewer
        #palette <- brewer_pal(type="qual",palette=2)
        #palette <- c("#d94313","#3097ff","#f5bd4e","#e259db","#009a3b","#da0b96","#38e096","#ff4471","#007733","#ff90f5","#588400","#feaedc","#a1d665","#42c7ff","#6f5500","#01b1be") 
palette <- brewer.pal(6,"Spectral")
hcr.colors <- palette[c(6,5,4,3,1,2)]
#show_col(hcr.colors) # C1 (Oc), C2 (Len), C3, constF, stability-favoring, trend-based (this is the order of the colors)
all.summaries <- NA
all.summaries <- lapply(results,FUN = summ.tab)   # This will take a while
all.summaries <- do.call(rbind.data.frame, all.summaries)
all.summaries$scenario <- rep(1:nscenarios,each=length(performance.measures))

# Match the scenarios to type of error, etc.
all.summaries <- merge(all.summaries,scen.table[,1:7],by="scenario")    # all.summaries is a giant table with 1080 rows = 72 scenarios * 15 PMs
all.summaries <- mutate(all.summaries, obs.error.type = recode(obs.error.type, 
                                                               'Tim'='Delayed change detection',
                                                               'AC' = "Autocorrelated"),
                        HCR = recode(HCR, 'cfp' = 'Stability-favoring',
                                     'constF' = 'Constant F',
                                     'C1' = 'C1',
                                     'C2' = 'C2',
                                     'C3' = 'C3',
                                     'trend' = "Trend-based"))

all.summaries$HCR <- factor(all.summaries$HCR, levels = c("C1","C2","C3","Constant F","Stability-favoring","Trend-based")) # Reorder factors so they plot in alphabetical order, the way they were intended to be!

write.csv(all.summaries, file = paste(Type,"_AllSummaries.csv",sep=""))



for (s in 1:nscenarios){
  #**N** indicate metrics for which higher values mean worse performance (like SD(catch)) - these metrics are in scen.table as 1/x
  result.to.use <- results[[s]]
  scen.table[s,performance.measures[1]] <- raw.table[s,performance.measures[1]] <- median(rowMeans(result.to.use$total.catch[,calc.ind])) # calculate mean B over years to use in the index - the final number is the median (across all simulations) of the mean B
  nonzero.catch <- result.to.use$total.catch[,calc.ind]
  nonzero.catch <- ifelse(nonzero.catch<0.1,NA,nonzero.catch)
  mnz.catches <- rowMeans(nonzero.catch,na.rm=TRUE)
  #mnz.catches[ is.na(mnz.catches) ] <- NA
  scen.table[s,performance.measures[2]] <- raw.table[s,performance.measures[2]] <- median(mnz.catches,na.rm = TRUE) #"LTnonzeromeancatch"
  
  scen.table[s,performance.measures[3]] <- 1 / median(apply(X = result.to.use$total.catch[,calc.ind],FUN = sd,MARGIN = 1)) #"SDcatch" **N**
  raw.table[s,performance.measures[3]] <- median(apply(X = result.to.use$total.catch[,calc.ind],FUN = sd,MARGIN = 1)) 
  
  ######
  five.yr.closure.given1 <- n.multiyr.closures(result.to.use$total.catch[,calc.ind],threshold=0)$count5 / 
    n.multiyr.closures(result.to.use$total.catch[,calc.ind],threshold=0)$count1
  median.P5 <- median(five.yr.closure.given1,na.rm=TRUE) # This is now the # of 5-year closures given a 1-year closure
  if(all(n.multiyr.closures(result.to.use$total.catch[,calc.ind],threshold=0)$count5==0) &
     all(n.multiyr.closures(result.to.use$total.catch[,calc.ind],threshold=0)$count1==0)){
    median.P5 = 0
  }
  scen.table[s,performance.measures[4]] <- ifelse(median.P5==0, 1, 1/median.P5)
  raw.table[s,performance.measures[4]] <- median.P5 # Mean number of 5-yr closures given a 1-yr closure
  
  ######
  ten.yr.closure.given5 <- n.multiyr.closures(result.to.use$total.catch[,calc.ind],threshold=0)$count10 / 
    n.multiyr.closures(result.to.use$total.catch[,calc.ind],threshold=0)$count5
  median.P10 <- median(ten.yr.closure.given5,na.rm=TRUE) 
  # If both P(5 yr closure|closure) and P(10 yr closure|closure) are 0, then set to zero:
  if(all(n.multiyr.closures(result.to.use$total.catch[,calc.ind],threshold=0)$count10==0) &
     all(n.multiyr.closures(result.to.use$total.catch[,calc.ind],threshold=0)$count5==0)){
    median.P10 = 0
  }
  scen.table[s,performance.measures[5]] <- ifelse(median.P10==0, 1, 1/median.P10)
  raw.table[s,performance.measures[5]] <- median.P10 # Mean number of 10-yr closures given a 5-yr closure
  
  ######
  scen.table[s,performance.measures[6]] <- ifelse(median(apply(X = result.to.use$total.catch[,calc.ind],FUN = nzeroes,MARGIN = 1)) ==0, 1,
                                                  1 / median(apply(X = result.to.use$total.catch[,calc.ind],FUN = nzeroes,MARGIN = 1)) ) # "nyrs0catch" **N**
  raw.table[s,performance.measures[6]] <- median(apply(X = result.to.use$total.catch[,calc.ind],FUN = nzeroes,MARGIN = 1))
  
  scen.table[s,performance.measures[7]] <- 
    raw.table[s,performance.measures[7]] <- median(rowMeans(result.to.use$biomass.total.true[,calc.ind])) # "meanbiomass"
  
  scen.table[s,performance.measures[8]] <- 
    raw.table[s,performance.measures[8]] <- median(apply(X = result.to.use$biomass.total.true[,calc.ind],FUN = good4pred,MARGIN = 1, 
                                                         F0.x = result.to.use$no.fishing.tb[,calc.ind])) # "good4preds" - this is calculated from TOTAL biomass (including age 0)
  
  scen.table[s,performance.measures[9]] <- 1 / median(apply(X = result.to.use$biomass.total.true[,calc.ind],FUN = sd,MARGIN = 1)) # "SDbiomass" **N**
  raw.table[s,performance.measures[9]] <- median(apply(X = result.to.use$biomass.total.true[,calc.ind],FUN = sd,MARGIN = 1))  #Actual raw SD of Biomass
  
  yrs.bad <- apply(X = result.to.use$biomass.total.true[,calc.ind],FUN = bad4pred,MARGIN = 1, F0.x = result.to.use$no.fishing.tb[,calc.ind] ) # Number of years that are very bad for preds - length of vector is nsims 

  scen.table[s,performance.measures[10]] <-  ifelse (median(yrs.bad)==0, 1, 1 / median(yrs.bad) ) #  **N** Number of years below a very low threshold (<10% of long term unfished biomass) "very.bad4preds" **N**
  raw.table[s,performance.measures[10]] <- median(yrs.bad)
  scen.table[s,performance.measures[11]] <- 
  raw.table[s,performance.measures[11]] <- NA #median(rowMeans(result.to.use$depl[,calc.ind]))  # Median depletion - depracated now because results don't include depletion
  
  sw <- subset(all.summaries,scenario==s)
  raw.table[s, performance.measures[12]] <- sw[12,'med']
  scen.table[s, performance.measures[12]] <- ifelse(is.na(sw[12,'med']),NA, 1 / sw[12,'med'])    # max collapse length
  
  raw.table[s, performance.measures[13]] <- sw[13,'med'] # **N**
  scen.table[s, performance.measures[13]] <- ifelse(is.na(sw[13,'med']),NA, 1 / sw[13,'med'])   # max bonanza length
  
  raw.table[s,performance.measures[14]] <- scen.table[s,performance.measures[14]] <- sw[14,'med']   # mean bonanza length
  
  raw.table[s,performance.measures[15]] <- sw[15,'med']
  scen.table[s,performance.measures[15]] <- ifelse(is.na(sw[15,'med']),NA, 1 / sw[15,'med'])    # mean collapse length
  
  raw.table[]
  
} # This loop is a hot mess and needs to be optimized - can also take out "scen.table" assignments because they're all done below for the kite plots

write.csv(raw.table, file=paste(Type,"_outputs.csv",sep=""))


############################################################################
# PLOTS AND METRICS TO SHOW OUTPUTS ----------------------------
############################################################################


# Change labels of things in the table! --------------------------
scen.table <- mutate(scen.table, obs.error.type = recode(obs.error.type, 
                                                        'Tim'='Delayed change detection',
                                                        'AC' = "Autocorrelated"),
                                             HCR = recode(HCR, 'cfp' = 'Stability-favoring',
                                                          'constF' = 'Constant F',
                                                          'C1' = 'C1',
                                                          'C2' = 'C2',
                                                          'C3' = 'C3',
                                                          'trend' = "Trend-based"))
raw.table <- mutate(raw.table, obs.error.type = recode(obs.error.type, 
                                                         'Tim'='Delayed change detection',
                                                         'AC' = "Autocorrelated"),
                     HCR = recode(HCR, 'cfp' = 'Stability-favoring',
                                  'constF' = 'Constant F',
                                  'C1' = 'C1',
                                  'C2' = 'C2',
                                  'C3' = 'C3',
                                  'trend' = "Trend-based"))


######################################################################
###### MAKE A PDF WITH ALL THE OUTPUT FIGURES! #######################
######################################################################

pdf(paste(Type,"AllPlots",Sys.Date(),".pdf",sep=""),width = 10,height = 9,onefile = TRUE)
# Put control rules in order so they plot right
scen.table$HCR <- raw.table$HCR <- factor(scen.table$HCR, levels = c("C1","C2","C3","Constant F","Stability-favoring","Trend-based"))




# Kite plots showing tradeoffs --------------------------------------------
mat <- matrix(1:4,nrow=2,byrow = TRUE)
plotnames <- list()
steepnesses <- unique(scen.table$h)
obs.error.types <- unique(scen.table$obs.error.type)
nice.pms <- data.frame(original = colnames(scen.table[-(1:7)]),
                       polished = c("LT mean catch","LT mean nonzero catch",
                                    "SD(Catch)","Probability of \n 5-yr closure",
                                    "Number of \n 10-yr closures","Number of yrs \n w/ zero catch",
                                    "LT mean biomass","Number of yrs \n above pred threshold",
                                    "SD(Biomass)","Number of yrs \n below low pred threshold",
                                    "Mean depletion","Max collapse length","Max bonanza length",
                                    "Bonanza length","Collapse length","Probability of collapse","Collapse severity"))

#For final version of paper, want to just show one scenario - this is is the basic one!

#for(steep in 1:2){
#for(obs in 1:2){
steep=2
obs=1
    
    tab <- subset(raw.table,obs.error.type == obs.error.types[obs] & h == steepnesses[steep] & M.type == "constant")
    tab.metrics <- tab[,-(1:7)]
    #Here are the PMs that are NEGATIVES (i.e., a high value for these is bad news)
    bad.pms <- c("SDcatch","n.5yrclose","n.10yrclose","nyrs0catch","SDbiomass","very.bad4preds","overallMaxCollapseLength","CollapseLength")
    which.bad <- which(colnames(tab.metrics) %in% bad.pms)
    tab.metrics[,which.bad] <- apply(tab.metrics[,which.bad],MARGIN = 2,FUN = function(x) ifelse(x==0,1,1/x)) # Turn all the "bad" PMs to their inverse
    props <- tab.metrics
    maxes <- apply(X = tab.metrics,MARGIN = 2,FUN = max,na.rm  = T )
    for(i in 1:nrow(props)){
      props[i,] <- tab.metrics[i,] / maxes
    }
    final.tab <- cbind(tab[,'HCR'],props)
    colnames(final.tab)[1] <- "group"
    test.nas <- apply(X = final.tab,FUN = anyNA,MARGIN = 2)
    na.metrics <- names(which(test.nas))
    
    if(length(na.metrics>0)){
      print(paste("The following performance metrics had NAs and was removed from the figure: ",na.metrics))
      final.tab <- final.tab[,-which(test.nas)]
      rm.metrics <- which(nice.pms$original %in% na.metrics)
      axis.labels <- nice.pms[-rm.metrics,'polished']
    }
    
    #legend.presence <- ifelse(mat[steep,obs] != 1,FALSE,TRUE)
    legend.presence = TRUE
    remove.these <- c("n.10yrclose","SDbiomass","meanDepl","LTnonzeromeancatch","good4preds","very.bad4preds")
    remove.ind <- which(colnames(final.tab) %in% remove.these)
    final.tab <- final.tab[-remove.ind]
    axis.labels <- nice.pms$polished[-(remove.ind-1)]

    
  pdf(file=paste(Type, "_Kite.pdf",sep=""),width = 10,height=7,useDingbats = FALSE)  
  ggradar(final.tab,font.radar = "Helvetica",grid.label.size=3,axis.label.size=4,
                                           legend.text.size = 4,
                                           axis.labels = axis.labels,
                                           plot.legend=legend.presence,palette.vec = hcr.colors)
dev.off()
    #}}

#grid.arrange(plotnames[[1]],plotnames[[2]],plotnames[[3]], plotnames[[4]])


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


# Plot just the predator ones ---------------------------------------------
            pred.df <- subset(all.summaries,PM %in% c("good4preds","very.bad4preds","overallMaxCollapseLength","overallMaxBonanzaLength","BonanzaLength","CollapseLength"))
            dodge <- position_dodge(.6)
            ggplot(pred.df,aes(x=HCR,y=med,colour=HCR,shape=obs.error.type,alpha=obs.error.type,label=med)) +
              scale_colour_manual(values = hcr.colors) +
              scale_alpha_manual(values = c(0.6,1)) +
              geom_point(size=5,position=dodge)  + 
              geom_errorbar(aes(ymin = loCI, ymax = hiCI), position = dodge,width=0.1) +
              #geom_text(nudge_x = 0.3) +
              theme_bw(base_size = 18) +
              theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5)) +
              facet_wrap(h~PM,scales="free_y")
            
            
# Plot comparing biomass, catch, etc for one run each scenario ------------
attributes <- c("Biomass","Catches","Recruitment","Depletion (B/B0)")
att.ind <- c(1,2,4,5)



# Plot time series of fishing rates leading up to collapses ---------------



# Plot first run for each scenario - DEPRACATED needs to be fixed for generic # scenarios
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
    ts <- results[[scenario.index]][["biomass.oneplus.true"]][1,calc.ind]
    plot(ts,col=hcr.colors[4],
         ylab = "Biomass", type="l", lwd=2,xlab="Year")
        lines(results[[scenario.index]][["biomass.oneplus.obs"]][1,calc.ind],col = hcr.colors[4],lty=2)
        text(x = 40,y=max(ts)*0.95,labels = paste(scen.table[scenario.index,'obs.error.type'],sep= " "))
        text(x = 40,y=max(ts)*0.85,labels = paste(scen.table[scenario.index,'h'],sep= " "))
    # The above is so messy but it's just to get the correct range for all the lines. 
    ts <- results[[scenario.index+4]][["biomass.oneplus.true"]][1,calc.ind]
    plot(ts,col=hcr.colors[3],
         ylab = "Biomass", type="l", lwd=2,xlab="Year")
    lines(results[[scenario.index+4]][["biomass.oneplus.obs"]][1,calc.ind],col=hcr.colors[3],lty=2) 
    text(x = 40,y=max(ts)*0.95,labels = paste(scen.table[scenario.index,'obs.error.type'],sep= " "))
    
    ts <- results[[scenario.index+8]][["biomass.oneplus.true"]][1,calc.ind]
    plot(ts,col=hcr.colors[1],
         ylab = "Biomass", type="l", lwd=2,xlab="Year")
    lines(results[[scenario.index+8]][["biomass.oneplus.obs"]][1,calc.ind],col=hcr.colors[1],lty=2) 
    text(x = 40,y=max(ts)*0.95,labels = paste(scen.table[scenario.index,'obs.error.type'],sep= " "))
    
    ts <- results[[scenario.index+12]][["biomass.oneplus.true"]][1,calc.ind]
    plot(ts,col=hcr.colors[2],
         ylab = "Biomass", type="l", lwd=2,xlab="Year")
    lines(results[[scenario.index+12]][["biomass.oneplus.obs"]][1,calc.ind],col=hcr.colors[2],lty=2) 
    text(x = 40,y=max(ts)*0.95,labels = paste(scen.table[scenario.index,'obs.error.type'],sep= " "))
    
    ts <- results[[scenario.index+16]][["biomass.oneplus.true"]][1,calc.ind]
    plot(ts,col=hcr.colors[5],
         ylab = "Biomass", type="l", lwd=2,xlab="Year")
    lines(results[[scenario.index+16]][["biomass.oneplus.obs"]][1,calc.ind],col=hcr.colors[5],lty=2)
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

pdf(paste(Type,"ErrorPlots",Sys.Date(),".pdf",sep=""),width=9.5,height=5)
# Plot just pred ones for paper -------------------------------------------
pred2.df <- subset(all.summaries, PM %in% c("good4preds","very.bad4preds") &
                     h ==0.9)
dodge <- position_dodge(.8)
ggplot(pred2.df,aes(x=HCR,y=med,colour=HCR,shape=obs.error.type,alpha=obs.error.type)) +
  scale_colour_manual(values = hcr.colors) +
  scale_alpha_manual(values = c(0.6,1)) +
  geom_point(size=5,position=dodge)  +
  geom_errorbar(aes(ymin = loCI, ymax = hiCI), position = dodge,width=0.1) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5)) +
  ylab("Median number of years") +
  facet_wrap(~PM,scales="free_y")
dev.off()


# Plot DD scenarios - for paper! ------------------------------------------
# For plotting the differences between AC error and delayed detection
# First, for menhaden, try scenario 2 (AC) vs. 4 (DD)-- biggest change in p(closure)
ac <- results[[2]]
dd <- results[[4]]
yrs <- 150:250
par(mfrow=c(2,1))
plot(ac$biomass.oneplus.true[1,yrs],type='l',ylab="Total biomass",xlab="Year",ylim=c(0,250000))
lines(ac$biomass.oneplus.obs[1,yrs],col='red')
lines(dd$biomass.oneplus.obs[1,yrs],col='blue')

lines(ac$total.catch[1,yrs],lty=2,col='red')
lines(dd$total.catch[1,yrs],lty=2,col='blue')

plot(ac$total.catch[1,yrs],type='l',col='red',ylab="Total catches",xlab="Year",ylim=c(0,250000))
lines(dd$total.catch[1,yrs],col='blue')

which(ac$total.catch[1,yrs] > ac$biomass.oneplus.true[1,yrs])
which(dd$total.catch[1,yrs] > dd$biomass.oneplus.true[1,yrs])
