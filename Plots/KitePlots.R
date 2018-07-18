# Figure 5

library(gtools)
library(plyr); library(dplyr)
library(scales)
library(gridExtra)
library(extrafont)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
source("~/Dropbox/Chapter4-HarvestControlRules/Code/ff-mse2/Plots/Megsieggradar.R")
figwd <- "~/Dropbox/Chapter4-HarvestControlRules/Figures"

# Define directories where each forage fish type is found
s.dir <- "~/Dropbox/Chapter4-HarvestControlRules/Results/Sardine2018-07-16"
a.dir <- "~/Dropbox/Chapter4-HarvestControlRules/Results/Anchovy2018-07-16"
m.dir <- "~/Dropbox/Chapter4-HarvestControlRules/Results/Menhaden2018-07-16"
dirs <- c(s.dir,a.dir,m.dir)
Types <- c("Sardine","Anchovy","Menhaden")

# Table for cleaning up axis labels
nice.pms <- data.frame(original = c("LTmeancatch", "LTnonzeromeancatch", "SDcatch", "n.5yrclose", 
                                    "n.10yrclose", "nyrs0catch", "meanbiomass", "good4preds", "SDbiomass", 
                                    "very.bad4preds", "meanDepl", "overallMaxCollapseLength", "overallMaxBonanzaLength", 
                                    "BonanzaLength", "CollapseLength", "Prob.Collapse", "Collapse.Severity", 
                                    "CV.Catch", "Bonafide.Collapse"),
                       polished = c("Mean catch","Mean nonzero catch",
                                    "Catch stability","Minimize \n P(5-yr \n closure)",
                                    "Minimize \n P(10-yr closure)","Minimize years \n w/ zero catch",
                                    "Mean biomass","Number of yrs \n above pred threshold",
                                    "SD(Biomass)","Number of yrs \n below low pred threshold",
                                    "Mean depletion","Max collapse length","Max bonanza length",
                                    "Bonanza length","Minimize \n collapse length", "Minimize \n P(collapse)",
                                    "Minimize \n collapse severity","CV(Catch)","Minimize \n long collapses"))

# Get all the data, compile it
all.dat <- data.frame()
for( ff in 1:3){
  setwd(dirs[ff])
  opfile <- grep("outputs.csv",x = list.files()) # Find outputs file - csv only
  raw.table <- read.csv(list.files()[opfile])
  if(colnames(raw.table)[1] == "X"){raw.table <- raw.table[,-1] } # if you use read.csv you need this
  raw.table$Type <- Types[ff]
  all.dat <- rbind(all.dat,raw.table)
  }

all.dat <- mutate(all.dat, HCR = recode_factor(HCR, 'cfp' = 'Stability-favoring',
                                                   'constF' = 'Low F',
                                                   'C1' = 'Basic hockey stick',
                                                   'C2' = 'Low Blim',
                                                   'C3' = 'High Fmax',
                                                   'constF_HI' = "High F"))

plotnames <- list() #list of subfig names

Figure5a <- function(steepness=0.6,obs="AC",M.type=constant){
    for(ff in 1:3){
      tab <- subset(all.dat,obs.error.type == obs & h == steepness & M.type == M.type & Type==Types[ff])
      tab.metrics <- tab[,-(c(1:7,27))]
      crs <- tab[,"HCR"]
      
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
      tab.metrics[,c("nyrs0catch")] <- 1-(tab$nyrs0catch/nyrs.to.use) # Proportion of years when there *wasn't* 0 catch
      props <- tab.metrics
      props <- apply(props, MARGIN = 2,FUN = function(x) x/max(x,na.rm=T))
      all(props<=1) # check to make sure everything worked
      
      final.tab <- data.frame(group = crs,props)
      test.nas <- apply(X = final.tab,FUN = anyNA,MARGIN = 2)
      na.metrics <- names(which(test.nas))
      
      if(length(na.metrics>0)){
        print(paste("The following performance metrics had NAs and were removed from the figure: ",na.metrics))
        final.tab <- final.tab[,-which(test.nas)]
        rm.metrics <- which(nice.pms$original %in% na.metrics)
        axis.labels <- nice.pms[-rm.metrics,'polished']
      }
      
      legend.presence <- FALSE#
      axis.ind <- match(x = colnames(final.tab)[-1],table = nice.pms$original)
      axis.labels <- nice.pms$polished[axis.ind]
      final.tab$group <- factor(final.tab$group,levels=c("Basic hockey stick","Low Blim","High Fmax","High F","Low F","Stability-favoring"))
      plotnames[[ff]] <- ggradar(final.tab,font.radar = "Helvetica",
                                grid.label.size=3,axis.label.size=7, 
                                legend.text.size = 4,
                                axis.labels = axis.labels,
                                plot.legend=legend.presence,palette.vec = hcr.colors,
                                manual.levels = levels(final.tab$group),
                                axis.label.offset=1.1)
    }

  p1 <- plotnames[[1]]
  p2 <- plotnames[[2]]
  p3 <- plotnames[[3]]
  grid.arrange(p1,p2,p3,ncol=3)
  } # end of Fig 5a function



setwd(figwd)
pdf(file = "Figure5a.pdf",width=27,height=15,useDingbats = FALSE)
Figure5a()
dev.off()

setwd(figwd)
pdf(file = "Figure5b.pdf",width=27,height=15,useDingbats = FALSE)
Figure5a(steepness = 0.9)
dev.off()

setwd(figwd)
pdf(file = "Figure5c.pdf",width=27,height=15,useDingbats = FALSE)
Figure5a(obs = "Tim")
dev.off()


