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
source("/Users/mcsiple/Dropbox/ChapterX-synthesis/Theme_Black.R")
Type = "Anchovy" #FF type to summarize
Date <- "2018-03-09"

# Set path to wherever the simulation results are:
path <- paste("/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Results/",Type,Date,"/",sep="")
# Anchovy: "/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Results/",Type,"2017-07-19","/",sep=""
# Menhaden: "/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Results/",Type,"2017-07-20","/",sep=""
# Sardine: "/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Results/",Type,"2017-07-20","/",sep=""

setwd(path)    # ONLY RDATA FILES and TXT FILES should be in this dir, otherwise you'll get an error

# Read all files into giant list
  files <- list.files(path=path)
  rm <- grep(files,pattern = ".txt") # Don't load the text table summary
  files <- files[-rm]
  files <- files[grep(files, pattern = ".RData")] #only load rdata files
  files <- mixedsort(files) # IMPORTANT: need this line to order in same order as scenario table!
  results <- sapply(files, function(x) mget(load(x)),simplify=TRUE) # This is a giant list of all results
  
  
  # NOTE: SKIP TO ~LINE 79 IF RESULTS HAVE ALREADY BEEN SUMMARIZED
  
  nscenarios <- length(results)
  raw.table <- read.table("Scenario_Table.txt")  #This empty table is generated when the simulations are run - then you fill it in after everything runs
  nsims <- nrow(results[[1]]$biomass.oneplus.true) # count nrows to know how many sims there are
  years.test <- ncol(results[[1]]$biomass.oneplus.true) # count cols to know how many years there are
  nyrs.to.use <- 100 # How many years you want to use to calculate all your metrics - There are no big differences btwn 50 and 100 yrs
  calc.ind <- tail(1:years.test, nyrs.to.use) # Which years to calculate median depletion over (length = nyrs.to.use)
  
  # Add performance measure columns to table
  performance.measures <- c("LTmeancatch","LTnonzeromeancatch","SDcatch","n.5yrclose","n.10yrclose","nyrs0catch","meanbiomass","good4preds","SDbiomass","very.bad4preds","meanDepl","overallMaxCollapseLength","overallMaxBonanzaLength","BonanzaLength","CollapseLength","Prob.Collapse","Collapse.Severity","CV.Catch","Bonafide.Collapse")
  pm.type <- c(rep("Fishery",times=6),rep("Ecosystem",times=11),"Fishery","Ecosystem") # for distinguishing types of PMs (mostly for plotting...)
  raw.table[,performance.measures] <- NA


# ------------------------------------------------------------------------
# Summarize everything in one giant "outputs" table - this is ugly, sorry
# ------------------------------------------------------------------------

# Fxns for summarizing and plotting ---------------------------------------

all.summaries <- NA
all.summaries <- lapply(results,FUN = summ.tab, calc.ind = calc.ind, ou.ind = NA)   # This will take a while
all.summaries <- do.call(rbind.data.frame, all.summaries)
all.summaries$scenario <- rep(1:nscenarios,each=length(performance.measures))

# Match the scenarios to type of error, etc.
all.summaries <- merge(all.summaries,raw.table[,1:7],by="scenario")    # all.summaries is a giant table with 1080 rows = 72 scenarios * 15 PMs
all.summaries2 <- all.summaries %>% mutate(obs.error.type = recode_factor(obs.error.type, 
                                                                   "Tim"="Delayed change detection",
                                                                   "AC" = "Autocorrelated"),
                                           HCR = recode_factor(HCR, 'cfp' = 'Stability-favoring',
                                                        'constF' = 'Low F',
                                                        'C1' = 'Basic hockey stick',
                                                        'C2' = 'Low Blim',
                                                        'C3' = 'High Fmax',
                                                        'trend' = "Trend-based",
                                                        'constF_HI' = "High F"))
                                           
subset(all.summaries2, h==0.6 & obs.error.type=="Autocorrelated" & PM == "CollapseLength") 
# 14 is basic hockey, 20 is low Blim, 26 is high Fmax
write.csv(all.summaries2, file = paste(Type,"_AllSummaries.csv",sep=""))



# Just load the raw table, if you have already run the code below  --------
    opfile <- grep("outputs",x = list.files()) # Find outputs file
     raw.table <- read.csv(list.files()[opfile])
    if(colnames(raw.table)[1] == "X"){raw.table <- raw.table[,-1] } # if you use read.csv you need this

# See if performance metrics are correlated -------------------------------
# sims.all <- lapply(results,FUN = summ.tab, individual.sim = TRUE)
# pairs(sims.all[[33]],pch=19,col=rgb(0,0,0,0.2))

# Table of summary stats! -----------------------------------------------------------
 
for (s in 1:nscenarios){
  #**N** indicate metrics for which higher values mean worse performance (like SD(catch)) - these metrics are in scen.table as 1/x
  result.to.use <- results[[s]]
  raw.table[s,performance.measures[1]] <- median(rowMeans(result.to.use$total.catch[,calc.ind],na.rm = TRUE)) 
  #calculate mean B over years to use in the index - the final number is the median (across all simulations) mean B
  nonzero.catch <- result.to.use$total.catch[,calc.ind]
  nonzero.catch <- ifelse(nonzero.catch<0.1,NA,nonzero.catch)
  mnz.catches <- rowMeans(nonzero.catch,na.rm=TRUE)
  raw.table[s,performance.measures[2]] <- median(mnz.catches,na.rm = TRUE) #"LTnonzeromeancatch"
  raw.table[s,performance.measures[3]] <- median(apply(X = result.to.use$total.catch[,calc.ind],FUN = sd,MARGIN = 1)) 
  ######
  five.yr.closure.given1 <- n.multiyr.closures(result.to.use$total.catch[,calc.ind],threshold=0)$count5 / 
                            n.multiyr.closures(result.to.use$total.catch[,calc.ind],threshold=0)$count1
  median.P5 <- median(five.yr.closure.given1,na.rm=TRUE) # This is now the # of 5-year closures given a 1-year closure
  if(all(n.multiyr.closures(result.to.use$total.catch[,calc.ind],threshold=0)$count5==0) &
     all(n.multiyr.closures(result.to.use$total.catch[,calc.ind],threshold=0)$count1==0)){
    median.P5 = 0
  }

  raw.table[s,performance.measures[4]] <- median.P5 # Mean number of 5-yr closures given a 1-yr closure
  
  ######
  ten.yr.closure.given5 <- n.multiyr.closures(result.to.use$total.catch[,calc.ind],threshold=0)$count10 / 
    n.multiyr.closures(result.to.use$total.catch[,calc.ind],threshold=0)$count5
  median.P10 <- median(ten.yr.closure.given5,na.rm=TRUE) 
  if(all(n.multiyr.closures(result.to.use$total.catch[,calc.ind],threshold=0)$count10==0) &
     all(n.multiyr.closures(result.to.use$total.catch[,calc.ind],threshold=0)$count5==0)){
    median.P10 = 0
  }
  
  raw.table[s,performance.measures[5]] <- median.P10 # Mean number of 10-yr closures given a 5-yr closure
  ######
  raw.table[s,performance.measures[6]] <- median(apply(X = result.to.use$total.catch[,calc.ind],FUN = nzeroes,MARGIN = 1))
  
    raw.table[s,performance.measures[7]] <- median(rowMeans(result.to.use$biomass.total.true[,calc.ind])) # "meanbiomass"

    raw.table[s,performance.measures[8]] <- median(apply(X = result.to.use$biomass.total.true[,calc.ind],
                                                         FUN = good4pred,MARGIN = 1, 
                                                         F0.x = result.to.use$no.fishing.tb[,calc.ind])) # "good4preds" - this is calculated from TOTAL biomass (including age 0)
    
  raw.table[s,performance.measures[9]] <- median(apply(X = result.to.use$biomass.total.true[,calc.ind],FUN = sd,MARGIN = 1))  #Actual raw SD of Biomass
  
  yrs.bad <- apply(X = result.to.use$biomass.total.true[,calc.ind],FUN = bad4pred,MARGIN = 1, F0.x = result.to.use$no.fishing.tb[,calc.ind] ) # Number of years that are very bad for preds - length of vector is nsims 

  raw.table[s,performance.measures[10]] <- median(yrs.bad)
  
  sw <- subset(all.summaries,scenario==s)
  for(pm in 11:19){
    select.ind <- which(sw$PM == performance.measures[pm]) 
    raw.table[s,performance.measures[pm]] <- sw[select.ind,'med'] 
  }
    # 1. Depletion
    # 2. max collapse length
    # 3. **N**max bonanza length
    # 4. Mean bonanza length
    # 5. Mean collapse length
    # 6. # probability of collapse 
    # 7. Collapse severity
    # 8. CV.Catch
    # 9. Bonafide.collapse
} 

write.csv(raw.table, file=paste(Type,Sys.Date(),"_outputs.csv",sep=""))

str(results)
subset(raw.table,h==0.6 & obs.error.type=="AC")


# Set colors! :) ----------------------------------------------------------
# Colour palette options for plots - some of these are from iWantHue and some are ColorBrewer
#palette <- brewer_pal(type="qual",palette=2)
#palette <- c("#d94313","#3097ff","#f5bd4e","#e259db","#009a3b","#da0b96","#38e096","#ff4471","#007733","#ff90f5","#588400","#feaedc","#a1d665","#42c7ff","#6f5500","#01b1be") 
palette <- brewer.pal(6,"Spectral")
hcr.colors <- palette[c(6,5,4,1,3,2)]
#show_col(hcr.colors) # C1, C2, C3, constF, stability-favoring (this is the order of the colors)
# Plot a few sample time series (can modify to look at specific issues)
par(mfrow=c(2,1))
plot(results[[1]]$intended.f[2,],type='l',ylab="Fishing rate",ylim=c(0,5))
lines(results[[1]]$fishing[2,],col='red')
plot(results[[1]]$intended.f[3,],type='l',ylab="Fishing rate",ylim=c(0,5))
lines(results[[1]]$fishing[3,],col='red')
par(mfrow=c(4,1))
plot(results[[2]]$biomass.oneplus.true[1,],type='l',ylab="Biomass")
lines(results[[2]]$total.catch[1,],col='red')
plot(results[[2]]$biomass.oneplus.true[2,],type='l',ylab="Biomass")
lines(results[[2]]$total.catch[2,],col='red')
plot(results[[2]]$biomass.oneplus.true[3,],type='l',ylab="Biomass")
lines(results[[2]]$total.catch[3,],col='red')
plot(results[[2]]$biomass.oneplus.true[4,],type='l',ylab="Biomass")
lines(results[[2]]$total.catch[4,],col='red')
par(mfrow=c(4,1))
plot(results[[5]]$biomass.oneplus.true[1,],type='l',ylab="Biomass")
lines(results[[5]]$total.catch[1,],col='red')
plot(results[[5]]$biomass.oneplus.true[2,],type='l',ylab="Biomass")
lines(results[[5]]$total.catch[2,],col='red')
plot(results[[5]]$biomass.oneplus.true[3,],type='l',ylab="Biomass")
lines(results[[5]]$total.catch[3,],col='red')
plot(results[[5]]$biomass.oneplus.true[4,],type='l',ylab="Biomass")
lines(results[[5]]$total.catch[4,],col='red')
# plot(results[[2]]$total.catch[4,],col='red',type='l')

# Look at catches to make sure they're sort of following the rule
par(mfrow=c(4,1))
for(i in 1:4){
plot(results[[2]]$fishing[i,],col='red',type='l')
}

# Look at scenarios with no obs error and make sure intended f = f
sim = 22
plot(results[[5]]$intended.f[sim,],type='l')
lines(results[[5]]$fishing[sim,],col='red')
plot(results[[6]]$intended.f[sim,],type='l')
lines(results[[6]]$fishing[sim,],col='red')


# Compare AC and DD
plot(results[[14]]$total.catch[1,],type='l',ylab="Catches") # 14 is autocorrelated error
lines(results[[16]]$total.catch[1,],col='red')              # 16 is delayed detection
plot(results[[14]]$total.catch[2,],type='l',ylab="Catches")
lines(results[[16]]$total.catch[2,],col='red')




# Plot example time series of 1+ biomass for each of the control rules-- together!
    #autcorrelated errors: 2,8,14,20,26,32
    #dd errors: 4,10,16,22,28,34
plot.these <- subset(raw.table, h==0.6 & obs.error.type =="AC")$scenario
calc.ind=1:250 #for when you want to look at the whole simulation!
par(mfrow=c(2,2))
for(sim in 21:24){
plot(results[[plot.these[5]]]$biomass.oneplus.true[sim,calc.ind],type='n',lwd=1.5,ylab="True 1+ Biomass")
lines(results[[plot.these[1]]]$biomass.oneplus.true[sim,calc.ind],col=hcr.colors[6],lwd=1.5)
lines(results[[plot.these[2]]]$biomass.oneplus.true[sim,calc.ind],col=hcr.colors[5],lwd=1.5)
lines(results[[plot.these[3]]]$biomass.oneplus.true[sim,calc.ind],col=hcr.colors[1],lwd=1.5)
lines(results[[plot.these[4]]]$biomass.oneplus.true[sim,calc.ind],col=hcr.colors[2],lwd=1.5)
lines(results[[plot.these[5]]]$biomass.oneplus.true[sim,calc.ind],col=hcr.colors[3],lwd=1.5)
lines(results[[plot.these[6]]]$biomass.oneplus.true[sim,calc.ind],col=hcr.colors[4],lwd=1.5)
}

for(sim in 21:24){
  plot(results[[plot.these[5]]]$total.catch[sim,calc.ind],col=hcr.colors[5],type='n',lwd=1.5,ylab="Total catches")
  lines(results[[plot.these[1]]]$total.catch[sim,calc.ind],col=hcr.colors[6],lwd=1.5)
  lines(results[[plot.these[2]]]$total.catch[sim,calc.ind],col=hcr.colors[5],lwd=1.5)
  lines(results[[plot.these[3]]]$total.catch[sim,calc.ind],col=hcr.colors[1],lwd=1.5)
  lines(results[[plot.these[4]]]$total.catch[sim,calc.ind],col=hcr.colors[2],lwd=1.5)
  lines(results[[plot.these[5]]]$total.catch[sim,calc.ind],col=hcr.colors[3],lwd=1.5)
  lines(results[[plot.these[6]]]$total.catch[sim,calc.ind],col=hcr.colors[4],lwd=1.5)
}


#  Why is collapse severity worse for basic hockey than for Low Bl --------
par(mfrow=c(1,1))
plot(results[[plot.these[1]]]$no.fishing.tb[sim,],col='black',type='l',lwd=1.5,ylab="1+ Biomass")
lines(results[[plot.these[3]]]$biomass.total.true[sim,],col=hcr.colors[1],lwd=1.5)
lines(results[[plot.these[4]]]$biomass.total.true[sim,],col=hcr.colors[2],lwd=1.5)
lines(results[[plot.these[5]]]$biomass.total.true[sim,],col=hcr.colors[3],lwd=1.5)

nexamples = 16
scenario = 6
dim = sqrt(nexamples)
par(mfrow=c(dim,dim))
for(i in (1:nexamples)){
  plot(results[[scenario]]$biomass.total.true[i,],type='l',ylab="B_total true") #Check 34 too... should have some crashes
}


############################################################################
# PLOTS AND METRICS TO SHOW OUTPUTS ----------------------------
############################################################################

#pdf(paste(Type,"AllPlots",Sys.Date(),".pdf",sep=""),width = 10,height = 9,onefile = TRUE)

raw.table <- mutate(raw.table, HCR = recode_factor(HCR, 'cfp' = 'Stability-favoring',
                                            'constF' = 'Low F',
                                            'C1' = 'Basic hockey stick',
                                            'C2' = 'Low Blim',
                                            'C3' = 'High Fmax',
                                            'constF_HI' = "High F"))

# Kite plots showing tradeoffs --------------------------------------------

nice.pms <- data.frame(original = colnames(raw.table[-(1:7)]),
                       polished = c("Mean catch","Mean nonzero catch",
                                    "Catch stability","Minimize \n P(5-yr closure)",
                                    "Minimize \n P(10-yr closure)","Minimize years \n w/ zero catch",
                                    "Mean biomass","Number of yrs \n above pred threshold",
                                    "SD(Biomass)","Number of yrs \n below low pred threshold",
                                    "Mean depletion","Max collapse length","Max bonanza length",
                                    "Bonanza length","Minimize \n collapse length", "Minimize \n P(collapse)",
                                    "Minimize \n collapse severity","CV(Catch)","Minimize \n long collapses"))

plotnames <- list()
tileplots <- list()
all.scaled <- vector()
plots <- data.frame("steepness"=c(0.6,0.9,0.6),"obs.error.type" = c("AC","AC","Tim"))

for(p in 1:3){
    tab <- subset(raw.table,obs.error.type == as.character(plots$obs.error.type[p]) & h == plots$steepness[p] & M.type == "constant")
    tab.metrics <- tab[,-(1:7)]
    crs <- tab[,"HCR"]

    #tab.metrics[is.na(tab.metrics)] <- 1
    remove.these <- c("n.10yrclose","SDbiomass","meanDepl","LTnonzeromeancatch","good4preds","very.bad4preds","CV.Catch","overallMaxCollapseLength","overallMaxBonanzaLength","Bonafide.Collapse")
    # Removed "Bonafide collapse" metric bc all CRs were performing similarly on it (in the paper this is called an "extended collapse")
    remove.ind <- which(colnames(tab.metrics) %in% remove.these)
    tab.metrics <- tab.metrics[-remove.ind]
    
    # Sometimes there won't be collapses, so for those turn NA's into zeroes (collapse length and severity)
    tab.metrics[,c("CollapseLength","Collapse.Severity")][is.na(tab.metrics[,c("CollapseLength","Collapse.Severity")])] <- 0
    #Here we deal with the PMs that are NEGATIVES (i.e., a high value for these is bad news) - I chose Option 3 below
    bad.pms <- c("SDcatch","n.5yrclose","n.10yrclose","nyrs0catch","SDbiomass","very.bad4preds","overallMaxCollapseLength","CollapseLength","Prob.Collapse","Collapse.Severity","CV(Catch)","Bonafide.Collapse")
    which.bad <- which(colnames(tab.metrics) %in% bad.pms)
    
    # For PMs with non-decimal values:
    tab.metrics[,c("SDcatch","CollapseLength")] <- 1 / tab.metrics[,c("SDcatch","CollapseLength")]
    tab.metrics$CollapseLength[is.infinite(tab.metrics$CollapseLength)] <- 1 # for when there are no collapses (get rid of infinite values)
    tab.metrics$SDcatch[is.infinite(tab.metrics$SDcatch)] <- 1 # there is also a case where SDcatch is zero because catches crash before the final 100 yrs
    
    # PMs with decimal values:
    tab.metrics[,c("Prob.Collapse","Collapse.Severity","n.5yrclose")] <- 1 - tab.metrics[,c("Prob.Collapse","Collapse.Severity","n.5yrclose")]
    tab.metrics[,c("nyrs0catch")] <- 1-(tab$nyrs0catch/nyrs.to.use) # This is essentially now the proportion of years when there *wasn't* 0 catch
    
    #tab.metrics[,-which.bad]
    props <- tab.metrics
    props <- apply(props, MARGIN = 2,FUN = function(x) x/max(x,na.rm=T))
    all(props<=1) # check to make sure everything worked
    
    # OPTION #1: this scales everything between 0 and 1... which is kind of weird because it can make little differences look HUGE
    #tab.metrics[,which.bad] <- apply(tab.metrics[,which.bad],MARGIN = 2,FUN = function(x) 1-(x/ifelse(all(x==0),1,max(x,na.rm=T)))) # Turn all the "bad" PMs to their inverse
    #tab.metrics[,-which.bad] <- apply(tab.metrics[,-which.bad],MARGIN = 2,FUN = function(x) (x-min(x))/(max(x)-min(x)))
    
    # OPTION #2: if you want to scale to the best performance but not necessarily between 0 (worst) and 1 (best)
    # PRscaled = 1- min(1, PR / Z)
    
    # OPTION #3: scale differently for different metrics. 
    # 1/x for all the "bad" performance measures that come in numbers >= 1: PRscaled = 1/PR
            # SD catch
            # Collapse length
    # 1-x for all the "bad" performance measures that come in numbers <1: PRscaled = 1-PR
            # Probability of collapse
            # Probability of a 5-yr closure given that there was a closure
            # Proportion of years with zero catch
            # Collapse severity
    # if the whole column is zero for any of these, make all PRscaled=1

    final.tab <- data.frame(group = crs,props)
    test.nas <- apply(X = final.tab,FUN = anyNA,MARGIN = 2)
    na.metrics <- names(which(test.nas))
    
    if(length(na.metrics>0)){
      print(paste("The following performance metrics had NAs and were removed from the figure: ",na.metrics))
      final.tab <- final.tab[,-which(test.nas)]
      rm.metrics <- which(nice.pms$original %in% na.metrics)
      axis.labels <- nice.pms[-rm.metrics,'polished']
    }
    
    legend.presence <- ifelse(p == 1,TRUE,FALSE)
    axis.ind <- match(x = colnames(final.tab)[-1],table = nice.pms$original)
    axis.labels <- nice.pms$polished[axis.ind]
    final.tab$group <- factor(final.tab$group,levels=c("Basic hockey stick","Low Blim","High Fmax","High F","Low F","Stability-favoring"))
  plotnames[[p]] <- ggradar(final.tab,font.radar = "Helvetica",   # Add "_b" to fxn name if making w black background
                            grid.label.size=3,axis.label.size=8, 
                                           legend.text.size = 4,
                                            axis.labels = axis.labels,
                                           plot.legend=legend.presence,palette.vec = hcr.colors,
                            manual.levels = levels(final.tab$group))
  
  ftm <- melt(final.tab,id.vars="group")
  ftm$name <- ftm$variable
  ftm$name <- factor(ftm$name, levels=c("LTmeancatch", "meanbiomass", "BonanzaLength","SDcatch","nyrs0catch","n.5yrclose","CollapseLength","Prob.Collapse","Collapse.Severity"))
  ftm$group <- factor(ftm$group,c("Basic hockey stick","Low Blim","High Fmax","Low F","High F","Stability-favoring"))
  ftm <- mutate(ftm,name = recode_factor(name, 'LTmeancatch' = "Mean catch",
                                     "meanbiomass" = "Mean biomass",
                                     "BonanzaLength" = "Bonanza Length",
                                     'SDcatch' = "Minimize \n catch variation",
                                     'nyrs0catch' = "Minimize \n years with 0 catch",
                                     'n.5yrclose' = "Minimize \n P(5 yr closure|closure)",
                                     "CollapseLength" = "Minimize \n collapse length",
                                     "Prob.Collapse" = "Minimize P(collapse)",
                                     "Collapse.Severity" = "Minimize collapse severity"))
  tileplots[[p]] <- ggplot(ftm,aes(x=name,y=group)) + geom_tile(aes(fill=value)) + scale_fill_distiller(palette="RdYlBu",trans="reverse") + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5)) + geom_text(aes(label=round(value,digits = 1)))
  all.scaled <- rbind(all.scaled,final.tab)
}

pdf(file = paste(Type,Sys.Date(),"2_KitePlots.pdf",sep=""),width = 15,height=27,useDingbats = FALSE)
grid.arrange(plotnames[[1]],plotnames[[2]],plotnames[[3]],ncol=1)
dev.off()

pdf(file = paste(Type,Sys.Date(),"_TilePlots.pdf",sep=""),width = 8,height = 14,useDingbats = FALSE)
grid.arrange(tileplots[[1]],tileplots[[2]],tileplots[[3]],ncol=1)
dev.off()


# Try a pairs plot --------------------------------------------------------

all.scaled$scen <- rep(c(1,2,3),each=length(unique(raw.table$HCR))) # This is different from "scenario" in tables above
all.scaled$cols <- rep(hcr.colors[c(6,5,1,2,3,4)],times=length(unique(all.scaled$scen)))
pairsnames <- c("Base case","High steepness","Delayed detection")
pdf(paste("PAIRS_",Type,".pdf",sep=""), width=11,height=8.5,onefile = TRUE)
for(i in 1:3){
p1 <- subset(all.scaled,scen==i)
par(las=1)
rm <- which(colnames(p1) %in% c("group","scen","cols")) # take out cols without performance in them
pairs(p1[,-rm],col=p1$cols,pch=19,xlim=c(0,1),ylim=c(0,1),labels=axis.labels,lower.panel = NULL)
#title(paste(pairsnames[i]),line = 0)
}
dev.off()

# Change labels of things in the table! --------------------------
raw.table <- mutate(raw.table, obs.error.type = recode(obs.error.type, 
                                                       'Tim'='Delayed change detection',
                                                       'AC' = "Autocorrelated"),
                    HCR = recode(HCR, 'cfp' = 'Stability-favoring',
                                 'constF' = 'Constant F',
                                 'C1' = 'C1',
                                 'C2' = 'C2',
                                 'C3' = 'C3',
                                 'trend' = "Trend-based"))

# Plot DD scenarios & check random seeds to make sure they match ------------
# For plotting the differences between AC error and delayed detection
# First, for menhaden, try scenario 2 (AC) vs. 4 (DD)-- biggest change in p(closure)
ac <- results[[2]]
dd <- results[[4]]
yrs <- 1:250
par(mfrow=c(1,1))
sim <- 5

# True biomass
yrange <- range(c(0,ac$biomass.total.true[sim,yrs],dd$biomass.total.true[sim,yrs]))
plot(ac$biomass.oneplus.true[sim,yrs],type='l',ylab="Total biomass",xlab="Year",ylim=yrange,axes=FALSE)
lines(dd$biomass.oneplus.true[sim,yrs],col='#2b83ba')

# Observations
lines(ac$biomass.oneplus.obs[sim,yrs],lty=2)
lines(dd$biomass.oneplus.obs[sim,yrs],lty=2,col='#2b83ba')

legend("topright",c("Autocorrelated","Delayed detection"),col=c('black','#2b83ba'),lwd=c(1,1),border = "white")
axis(1)

plot(ac$biomass.oneplus.true[5,yrs],type='l')
lines(ac$no.fishing.tb[5,yrs],col='green')
lines(dd$no.fishing.tb[5,yrs],col='green',lty=2)


lines(ac$total.catch[1,yrs],lty=2,col='red')
lines(dd$total.catch[1,yrs],lty=2,col='blue')

plot(ac$total.catch[1,yrs],type='l',col='red',ylab="Total catches",xlab="Year",ylim=c(0,5000))
lines(dd$total.catch[1,yrs],col='blue')


# Why the increase in catches with the stability-favoring rule ------------
#... when there is delayed detection?
ac <- results[[2]]
dd <- results[[4]]
# OR (high steepness):
# ac <- results[[1]]
# dd <- results[[3]]

yrs <- 1:250
par(mfrow=c(1,1))
sim <- 14
plot(ac$total.catch[sim,yrs],col=hcr.colors[6],type='l',lwd=2)
lines(dd$total.catch[sim,yrs],col=hcr.colors[6],lwd=2,lty=2)

ac.means <- apply(ac$total.catch[,yrs],MARGIN = 1,FUN = mean)
dd.means <- apply(dd$total.catch[,yrs],MARGIN = 1,FUN = mean)

median(ac.means)
median(dd.means)


###########################################################################
###########################################################################
# Everything below this line can be used but isn’t really necessary -------
###########################################################################
###########################################################################

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
# *** NOTE: Need to scale these to best (as in kite plots); scen.table depracated
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



#  Steepness tests --------------------------------------------------------

hiS <- results[[7]]
loS <- results[[8]]
yrs <- 150:250
par(mfrow=c(2,1))
plot(hiS$biomass.total.true[1,yrs],type='l',ylab="Total biomass",xlab="Year")
lines(loS$biomass.total.true[1,yrs],col='red')


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
