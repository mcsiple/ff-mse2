library(car)
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

resultsdir <- "~/Dropbox/Chapter4-HarvestControlRules/Results"

# Before starting: Check to make sure these directories are results you want
aa <- read.csv(file.path(resultsdir,"Anchovy2017-10-05/Anchovy2017-10-05_outputs.csv"))
aa$Type = "Anchovy"
mm <-read.csv(file.path(resultsdir,"Menhaden2017-10-05/Menhaden2017-10-09_outputs.csv"))
mm$Type = "Menhaden"
ss <- read.csv(file.path(resultsdir,"Sardine2017-10-05/Sardine2017-10-09_outputs.csv"))
ss$Type = "Sardine"
dat <- rbind.fill(aa,mm,ss) # This is an amalgamation of all the raw.tables from the results files

dat <- dat[,-1] # Remove row numbers
head(dat)
dat <- subset(dat, h==0.6,obs.error.type != "noerror")
remove <- c("scenario","M.type","h","recruit.rho","recruit.sd")
dat <- dat[,-which(colnames(dat) %in% remove)]
mdat <- melt(dat,id.vars=c("Type","obs.error.type","HCR"))


#dat2 <- dcast(mdat,variable+value ~ Type+obs.error.type+HCR+recruit.sd+recruit.rho)
dat2 <- dcast(mdat,Type+HCR+variable ~ obs.error.type)
#dat2 <- dcast(mdat,value ~ Type+variable+obs.error.type+HCR+recruit.sd+recruit.rho)
dat2$percentdiff <- ((dat2$Tim - dat2$AC) / dat2$Tim) * 100
dat3 <- subset(dat2,!variable %in% c("overallMaxCollapseLength","LTnonzeromeancatch","n.10yrclose","meanDepl","good4preds","very.bad4preds","Bonafide.Collapse","overallMaxBonanzaLength","SDbiomass","CV.Catch"))

dat3 <- mutate(dat3,HCR = recode_factor(HCR, 'cfp' = 'Stability-favoring',
                            'C1' = 'C1',
                            'C2' = 'C2',
                            'C3' = 'C3',
                            "constF" = "Constant F - low",
                            'constF_HI' = "Constant F - high"))

write.csv(dat3,"PercentDiffs_101017.csv")

# Start here if you’ve already calculated percent differences -------------
setwd("~/Dropbox/Chapter4-HarvestControlRules/Figures/")
#dat3 <- read.csv("PercentDiffs.csv",header=T)
palette <- brewer.pal(6,"Spectral")
hcr.colors <- palette[c(6,5,4,3,1,2)]


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

# Order factors for plotting
dat3.new$name <- factor(dat3.new$name, levels=rev(c("Bonanza length","Mean biomass","Mean catch", "Minimize collapse severity","Minimize P(collapse)","Minimize collapse length","Minimize P(5 yr closure|closure)","Minimize years with 0 catch","Minimize catch variation")))
dat3.new$HCR <- factor(dat3.new$HCR, levels=c("C1","C2","C3","Constant F - low","Constant F - high","Stability-favoring"))

setwd("/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Figures")
pdf(file = "PercentDiffsErrors_101017.pdf",width = 11,height = 9,useDingbats = FALSE)
ggplot(dat3.new, aes(x=name,y=percentdiff)) +
  geom_bar(colour='black',aes(fill=HCR),stat = "identity") + 
  scale_fill_manual(values = hcr.colors[c(1,2,3,4,5,6)]) +
  geom_hline(yintercept=0)+facet_grid(HCR~Type) +
  theme_classic(base_size = 14) + 
  theme(strip.background = element_blank(),strip.text.y = element_blank()) +
  ylab("% change from delayed detection model") +
  xlab("Performance metric") + 
  ylim(c(-100,100)) + coord_flip()
dev.off()



#  Make a dark-background one for presentations ---------------------------
pdf(file = "PercentDiffsErrors_black.pdf",width = 11,height = 9,useDingbats = FALSE)
ggplot(dat3.new, aes(x=name,y=percentdiff)) +
  geom_bar(colour='black',aes(fill=HCR),stat = "identity") + 
  scale_fill_manual(values = hcr.colors[c(1,2,3,6,4,5)]) +
  geom_hline(yintercept=0,colour='white')+facet_grid(HCR~Type) +
  #geom_vline(xintercept = 6.5) +
  #theme_classic(base_size = 14) + 
  theme_black(base_size = 14) +
  theme(strip.background = element_blank(),strip.text.y = element_blank()) +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2)) +
  ylab("% change from delayed detection model") +
  xlab("Performance metric") + 
  ylim(c(-100,100)) + coord_flip() 
dev.off()
