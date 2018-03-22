library(car)
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

resultsdir <- "~/Dropbox/Chapter4-HarvestControlRules/Results"

# *** Before starting: Check to make sure these directories are results you want
aa <- read.csv(file.path(resultsdir,"Anchovy2018-03-09/Anchovy2018-03-19_outputs.csv"))
aa$Type = "Anchovy"
mm <-read.csv(file.path(resultsdir,"Menhaden2018-03-09/Menhaden2018-03-19_outputs.csv"))
mm$Type = "Menhaden"
ss <- read.csv(file.path(resultsdir,"Sardine2018-03-09/Sardine2018-03-19_outputs.csv"))
ss$Type = "Sardine"
dat <- rbind.fill(aa,mm,ss)     # This is an amalgamation of all the raw.tables from the results files

dat <- dat[,-1] # Remove row numbers
head(dat)
dat <- subset(dat, h==0.6,obs.error.type != "noerror")
remove <- c("scenario","M.type","h","recruit.rho","recruit.sd")
dat <- dat[,-which(colnames(dat) %in% remove)]
mdat <- melt(dat,id.vars=c("Type","obs.error.type","HCR"))


dat2 <- dcast(mdat,Type+HCR+variable ~ obs.error.type)
dat2$percentdiff <- ((dat2$Tim - dat2$AC) / dat2$AC) * 100 # Calc % diff between delay detection and AC 

dat3 <- subset(dat2,!variable %in% c("overallMaxCollapseLength","LTnonzeromeancatch","n.10yrclose","meanDepl","good4preds","very.bad4preds","Bonafide.Collapse","overallMaxBonanzaLength","SDbiomass","CV.Catch"))

dat3 <- mutate(dat3,HCR = recode_factor(HCR, 'cfp' = 'Stability-favoring',
                                 'constF' = 'Low F',
                                 'C1' = 'Basic hockey stick',
                                 'C2' = 'Low Blim',
                                 'C3' = 'High Fmax',
                                 'trend' = "Trend-based",
                                 'constF_HI' = "High F"))
setwd("~/Dropbox/Chapter4-HarvestControlRules/Tables")
for(i in 1:nrow(dat3)){
  if(dat3$AC[i]==0 & dat3$Tim[i]==0){
    dat3$percentdiff[i] <- 0    
  }
  if(is.infinite(dat3$percentdiff[i])){
    dat3$percentdiff[i] <- 100
  }
}
write.csv(dat3,"PercentDiffs_032118.csv")

# Start here if you’ve already calculated percent differences -------------
setwd("~/Dropbox/Chapter4-HarvestControlRules/Figures/")
#dat3 <- read.csv("PercentDiffs.csv",header=T)
palette <- brewer.pal(6,"Spectral")
hcr.colors <- palette[c(6,5,4,1,3,2)]


# Multiply all the "negative" performance measures by -1 so they're all scaled to make positive = good
negativos <- c("SDcatch","nyrs0catch","n.5yrclose","CollapseLength","Prob.Collapse","Collapse.Severity")
dat3.new <- dat3
dat3.new$percentdiff[which(dat3$variable %in% negativos)] <- -1 * dat3.new$percentdiff[which(dat3$variable %in% negativos)]
unique(dat3$variable)
dat3.new$name <- dat3.new$variable
# Change names to nice labels so they make sense
dat3.new <- mutate(dat3.new,name = recode_factor(name, 'LTmeancatch' = "Mean catch",
                                       "meanbiomass" = "Mean biomass",
                                       "BonanzaLength" = "Bonanza length",
                                       'SDcatch' = "Minimize catch variation",
                                       'nyrs0catch' = "Minimize years with 0 catch",
                                       'n.5yrclose' = "Minimize P(5 yr closure|closure)",
                  "CollapseLength" = "Minimize collapse length",
                  "Prob.Collapse" = "Minimize P(collapse)",
                  "Collapse.Severity" = "Minimize collapse severity"))

# Order control rules and performance measures for plotting
dat3.new$name <- factor(dat3.new$name, levels=rev(c("Mean biomass","Bonanza length","Minimize collapse length","Minimize P(collapse)","Minimize collapse severity","Mean catch","Minimize catch variation","Minimize P(5 yr closure|closure)","Minimize years with 0 catch")))
dat3.new$HCR <- factor(dat3.new$HCR, levels=c("Basic hockey stick","Low Blim","High Fmax","High F","Low F","Stability-favoring"))
dat3.new$Type <- factor(dat3.new$Type, levels = c("Sardine","Anchovy","Menhaden"))

# Some metrics have a >100% change and won't show up if you don't set them w/in plot limits. 
# make a new column with labels for these 
dat3.new$percentlabel <- NA

# Fill in the labels from percentdiffs column
dat3.new$percentlabel[which(dat3.new$percentdiff > 100)] <- dat3.new$percentdiff[which(dat3.new$percentdiff > 100)]
dat3.new$percentlabel[which(dat3.new$percentdiff < -100)] <- dat3.new$percentdiff[which(dat3.new$percentdiff < -100)]
dat3.new$percentlabel <- round(dat3.new$percentlabel,digits = 0)
dat3.new$location <- 100
dat3.new$location[which(dat3.new$percentdiff< (-100))] <- -100

dat3.new$percentdiff[dat3.new$percentdiff > 100 & dat3.new$percentdiff < 1e6] = 100
dat3.new$percentdiff[dat3.new$percentdiff < (-100) & dat3.new$percentdiff > (-1e6)] = -100
dat3.new$percentdiff[is.infinite(dat3.new$percentdiff)] <- NA
dat3.new$percentlabel[which(is.infinite(dat3.new$percentlabel))] <- NA # replace the two "inf" values with NA

setwd("/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Figures")
pdf(file = "PercentDiffsErrors_032118.pdf",width = 11,height = 9,useDingbats = FALSE)
ggplot(dat3.new, aes(x=name,y=percentdiff)) +
  geom_bar(colour='black',aes(fill=HCR),stat = "identity") + 
  scale_fill_manual(values = hcr.colors) +
  geom_text(aes(x=name,y=location, label=percentlabel,hjust=ifelse(sign(percentlabel)>0, -0.2, 1.1)), 
            position = position_dodge(width=1),size=2.5) +
  geom_hline(yintercept=0)+facet_grid(HCR~Type) +
  theme_classic(base_size = 14) + 
  theme(strip.background = element_blank(),strip.text.y = element_blank()) +
  ylab("% change from delayed detection model") +
  xlab("Performance metric") + 
  ylim(c(-130,130)) + coord_flip()
dev.off()



# Maybe add line plots to the top with time series? -----------------------
#sardine.eg <- # just pick out the scenarios you want, load results




#  Make a dark-background one for presentations ---------------------------
pdf(file = "PercentDiffsErrors_black.pdf",width = 11,height = 9,useDingbats = FALSE)
ggplot(dat3.new, aes(x=name,y=percentdiff)) +
  geom_bar(colour='black',aes(fill=HCR),stat = "identity") + 
  scale_fill_manual(values = hcr.colors) +
  geom_hline(yintercept=0,colour='white')+facet_grid(HCR~Type) +
  theme_black(base_size = 14) +
  theme(strip.background = element_blank(),strip.text.y = element_blank()) +
  ylab("% change from delayed detection model") +
  xlab("Performance metric") + 
  ylim(c(-100,100)) + coord_flip() 
dev.off()
