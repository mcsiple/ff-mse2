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
#Testing
# Functions for calculating performance measures ---------------------------

nzeroes <- function(x){ # How many years with zero catch
  #x is a vector
  n <- length(which(x==0))
  return(n)
}

good4pred <- function(x){ # Number of years that are above a certain threshold (here, it's 25% of the long term mean biomass)
  ltm <- mean(x)
  thresh <- 0.75*ltm
  g4p <- length(which(x>thresh))
  return(g4p)
}

bad4pred <- function(x){ # Number of years that are below a certain threshold (here, it's 10% of the long term mean biomass)
  # x is a time series of biomass
  ltm <- mean(x)
  thresh <- 0.1*ltm
  b4p <- length(which(x<thresh))
  return(b4p)
}

n.multiyr.closures <- function(x, threshold = NA) { #where x is a matrix, rows are sims, cols are years
  count5 <- count10 <- vector(length=nrow(x))
  for(i in 1:nrow(x)){ #Either catch OR biomass
    ltm <- mean(x[i,])
    if(is.na(threshold)){ thresh <- 0.01*ltm } # threshold can be anything - CHANGE THIS TO BE DIFF FOR DIFF HCRS!!!
    else{thresh = threshold}
    badTorF <- x[i,] <= thresh
    fiveyr <- sum(roll_sum(badTorF, 5) == 5)
    tenyr <- sum(roll_sum(badTorF, 10) == 10)
    count5[i] <- fiveyr #number of five year closures
    count10[i] <- tenyr # number of 10 year closures
  }
  
  return(list(count5 = count5,count10 = tenyr)) # Mean number of 5- and 10-yr closures
}

# Rcpproll demo:
# set.seed(1); x <- sample(c(T, F), 100, replace = T); sum(RcppRoll::roll_sum(x, 3) == 3)

# Set path to wherever the simulation results are, load them into a giant dataframe
path <- "/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Results/Anchovy/"
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
nyrs.to.use <- 100 # How many years you want to use to calculate all your metrics
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
  scen.table[s,performance.measures[2]] <- raw.table[s,performance.measures[2]] <- median(rowMeans(nonzero.catch,na.rm=TRUE)) #"LTnonzeromeancatch"
  
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
    raw.table[s,performance.measures[7]] <- median(rowMeans(result.to.use$biomass[,calc.ind])) # "meanbiomass"
  
  scen.table[s,performance.measures[8]] <- 
    raw.table[s,performance.measures[8]] <- median(apply(X = result.to.use$biomass[,calc.ind],FUN = good4pred,MARGIN = 1)) # "good4preds"
  
  scen.table[s,performance.measures[9]] <- 1 / median(apply(X = result.to.use$biomass[,calc.ind],FUN = sd,MARGIN = 1)) # "SDbiomass" **N**
  raw.table[s,performance.measures[9]] <- median(apply(X = result.to.use$biomass[,calc.ind],FUN = sd,MARGIN = 1))  #Actual raw SD of Biomass
  
  yrs.bad <- apply(X = result.to.use$biomass[,calc.ind],FUN = bad4pred,MARGIN = 1) # length of vector is nsims 
  scen.table[s,performance.measures[10]] <- ifelse(length(which(yrs.bad>0))/nsims == 0, 0,
                                               1 / length(which(yrs.bad>0))/nsims ) #Prob that biomass ever dipped below 10% of LT mean "p.bad4preds" **N**
  raw.table[s,performance.measures[10]] <- length(which(yrs.bad>0))/nsims
  scen.table[s,performance.measures[11]] <- 
    raw.table[s,performance.measures[11]] <- median(rowMeans(result.to.use$depl[,calc.ind])) 
}

write.csv(raw.table, file="Anchovy_outputs.csv")


############################################################################
# PLOTS AND METRICS TO SHOW OUTPUTS ----------------------------
############################################################################

# Fxns for summarizing and plotting ---------------------------------------
## Add an alpha value to a colour (from Mages' blog, http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

# Colour palette for time series plots - some of these are from iWantHue and some are ColorBrewer
#palette <- brewer_pal(type="qual",palette=2)
#palette <- c("#d94313","#3097ff","#f5bd4e","#e259db","#009a3b","#da0b96","#38e096","#ff4471","#007733","#ff90f5","#588400","#feaedc","#a1d665","#42c7ff","#6f5500","#01b1be") 
palette <- brewer.pal(6,"Spectral")
show_col(palette)
hcr.colors <- palette[c(6,5,3,1,2)]
show_col(hcr.colors) # C1 (Oc), C2 (Len), constF, stability-favoring, trend-based (this is the order of the colors)

# Function to plot medians and certainty intervals from simulations:
plotintervals <- function(result.mat,ylab){ #result.mat is a matrix (e.g., biomass for results[[1]])
  median.vec <- apply(result.mat,MARGIN = 2,FUN = median)
  ints <- apply(result.mat,MARGIN = 2,FUN = quantile, probs = c(0.025,0.25,0.75,0.975))
  lo95 <- ints[1,calc.ind]
  hi95 <- ints[4,calc.ind]
  lo50 <- ints[2,calc.ind]
  hi50 <- ints[3,calc.ind]
  
  plot(1:nyrs.to.use,median.vec[calc.ind],
       type='l',lty=1,lwd=2,
       ylim=range(lo95,hi95),
       ylab = ylab,
       xlab = "Year")
  
  zz <- c(1:nyrs.to.use,tail(nyrs.to.use,n=1),rev(1:nyrs.to.use)) # for polygons
  aa <- c(hi95,0,rev(lo95))
  bb <-  c(hi50,0,rev(lo50))
  polygon(zz,aa,col=adjustcolor( "black", alpha.f = 0.2),border="NA")
  polygon(zz,bb,col=adjustcolor( "black", alpha.f = 0.2),border="NA")
}

#Function to get summary metrics so you can make Key plots (or whatever you want)
summ.tab <- function(result.list){ #result.list is one of the results (=1 harvest rule, 1000 time series of biomass, catch, fishing, rec, depl)
  for(i in 1:length(result.list)){
    result.list[[i]] <- result.list[[i]][,calc.ind]
  } # Trim results to the years we're using
  performance.measures
  LTmeans.list <- lapply(result.list,FUN = rowMeans) 
  # median and quantiles of LTM of all PMs
  ltm <- lapply(LTmeans.list,FUN = quantile, probs = c(0.05,0.5,0.95)) 
  # mean nonzero catch
  catch <- result.list$total.catch
  nonzero.catch <- ifelse(catch<0.1,NA,catch)
  ltm.nzc1 <- rowMeans(nonzero.catch,na.rm=TRUE)
  ltm.nzc2 <- quantile(ltm.nzc1,probs = c(0.05,0.5,0.95),na.rm = TRUE)
  # SDcatch
  sd.catches <- apply(catch, MARGIN = 1,FUN = sd, na.rm = TRUE)
  sd.catch <- quantile(sd.catches,probs = c(0.05,0.5,0.95))
  
  #5- and 10-yr closures
  n.5yrclose <- n.multiyr.closures(catch)$count5
  n.10yrclose <- n.multiyr.closures(catch)$count10
  
  #Number of years w zero catch
  nz1 <- apply(catch,MARGIN = 1,FUN = nzeroes)
  
  # SDbiomass
  biomass <- result.list$biomass
  sd.Bs <- apply(biomass, MARGIN = 1,FUN = sd, na.rm = TRUE)
  sd.B <- quantile(sd.Bs,probs = c(0.05,0.5,0.95))
  
  #Years that are "good for predators"
  g4p.vec <- apply(X = biomass,FUN = good4pred,MARGIN = 1)
  
  # Number of years that are below a predator threshold
  yrs.bad <- apply(X = biomass,FUN = bad4pred,MARGIN = 1) # length of vector is nsims 

  
  # Performance metrics
  interval <- c(0.05,0.5,0.95)
  
  ltm.c <- ltm$total.catch  # ltmcatch
  ltm.nzc <- ltm.nzc2       #ltmnonzerocatch
  SDcatch <- sd.catch       #SD(Catch)
  n5yr <- quantile(n.5yrclose,probs = interval)     #n.5yrclose
  n10yr <- quantile(n.10yrclose,probs = interval)   #n.10yrclose
  nz <- quantile(nz1,probs = interval)              #nyrs.0catch
  ltm.b <- ltm$biomass      #LTMBiomass
  g4p <- quantile(g4p.vec,probs = interval)         #Nyears "good for predator"
  sdB <- sd.B                               #SD(Biomass)
  b4p <- quantile(yrs.bad,probs = interval) #p(bad4preds)
  ltm.depl <- ltm$depl                      # Mean depletion
  
  output <- data.frame(PM = performance.measures, loCI = NA, med = NA, hiCI = NA)
  output[,-1] <- rbind(ltm.c,ltm.nzc,SDcatch,n5yr,n10yr,nz,ltm.b,g4p,sdB,b4p,ltm.depl)
  return(output)
}

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
all.summaries$HCR <- factor(all.summaries$HCR, levels = c("C1","C2","Constant F","Stability-favoring","Trend-based"))
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

pdf("Anchovy_May2.pdf",width = 10,height = 9,onefile = TRUE)
# Put control rules in order so they plot right
scen.table$HCR <- factor(scen.table$HCR, levels = c("C1","C2","Constant F","Stability-favoring","Trend-based"))
# Compare each of the CRs together? It would be like pairs()



# Kite plots showing tradeoffs --------------------------------------------
mat <- matrix(1:4,nrow=2,byrow = TRUE)
plotnames <- list()
steepnesses <- unique(scen.table$h)
obs.error.types <- unique(scen.table$obs.error.type)

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
    if(length(na.metrics>0)){
      print(paste("The following performance metrics had NAs and was removed from the figure: ",na.metrics))
      final.tab <- final.tab[,-which(test.nas)]
    }
    legend.presence <- ifelse(mat[steep,obs] != 1,FALSE,TRUE)
    plotnames[[mat[steep,obs]]] <- ggradar(final.tab,font.radar = "Helvetica",grid.label.size=4,axis.label.size=4,plot.legend=legend.presence,palette.vec = hcr.colors)
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

# Plot first run for each scenario
for (scenario.index in 1:4){
  par(mfrow=c(2,2),mar=c(5,4,3,2)+0.1)
  for(i in 1:length(att.ind)){
    which.metric <- att.ind[i]
    plot(results[[scenario.index]][[which.metric]][1,calc.ind],col=hcr.colors.lines[1],
         ylab = attributes[i], type="l", lwd=2,xlab="Year",
         ylim = range(results[[scenario.index]][[which.metric]][,calc.ind]))
    lines(results[[scenario.index+4]][[which.metric]][1,calc.ind],lwd=2,col=hcr.colors[2]) 
    lines(results[[scenario.index+8]][[which.metric]][1,calc.ind],lwd=2,col=hcr.colors[3]) 
    lines(results[[scenario.index+12]][[which.metric]][1,calc.ind],lwd=2,col=hcr.colors[4])
    lines(results[[scenario.index+16]][[which.metric]][1,calc.ind],lwd=2,col=hcr.colors[5])
  }
  
  mtext(paste("h=",scen.table[scenario.index,'h'],"-","ObsErrorType =",scen.table[scenario.index,'obs.error.type'],sep= " "),outer=TRUE,cex=1.2,side=3,line=-3)
  legend("topright",col = hcr.colors[1:5],lwd=rep(2,times=5),legend = sort(unique(scen.table$HCR)))
}


# Medians and 95% intervals for biomass, catch, etc ------------------
attributes <- c("Biomass","Catches","Recruitment","Depletion (B/B0)")
att.ind <- c(1,2,4,5)

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