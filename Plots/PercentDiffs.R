library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
dat <- read.csv("/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Results/FF_MSE_summary.csv")
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
#dat3 <- subset(dat3,HCR !="constF")
dat3 <- mutate(dat3,HCR = recode(HCR, 'cfp' = 'Stability-favoring',
                            'C1' = 'C1',
                            'C2' = 'C2',
                            'C3' = 'C3',
                            "constF" = "Constant F - low",
                            'constF_HI' = "Constant F - high"))

write.csv(dat3,"PercentDiffs.csv")
palette <- brewer.pal(6,"Spectral")
hcr.colors <- palette[c(6,5,4,3,1,2)]

# Reorder outputs so the good ones are separate from the bad ones
unique(dat3$variable)
dat3$name <- dat3$variable
dat3$name <- factor(dat3$name, levels=c("SDcatch","nyrs0catch","n.5yrclose","CollapseLength","Prob.Collapse","Collapse.Severity","LTmeancatch", "meanbiomass", "BonanzaLength"))

dat3 <- mutate(dat3,name = recode (name, 'LTmeancatch' = "Mean catch",
                                       "meanbiomass" = "Mean biomass",
                                       "BonanzaLength" = "Bonanza Length",
                                       'SDcatch' = "Catch variation",
                                       'nyrs0catch' = "Years with 0 catch",
                                       'n.5yrclose' = "P(5 yr closure|closure)",
                  "CollapseLength" = "Collapse length",
                  "Prob.Collapse" = "P(collapse)",
                  "Collapse.Severity" = "Collapse severity"))

setwd("/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Figures")
pdf(file = "PercentDiffsErrors.pdf",width = 10,height=9,useDingbats = FALSE)
ggplot(dat3, aes(x=name,y=percentdiff)) +
  geom_bar(colour='black',aes(fill=HCR),stat = "identity") + 
  scale_fill_manual(values = hcr.colors[c(1,2,3,6,4,5)]) +
  geom_hline(yintercept=0)+facet_grid(HCR~Type) +
  geom_vline(xintercept=6.5)+
  theme_classic(base_size = 14) + 
  theme(strip.background = element_blank(),strip.text.y = element_blank()) +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2)) +
  ylab("% change from delayed detection model") +
  xlab("Performance metric") + 
  ylim(c(-100,100)) + coord_flip()
dev.off()
#fill=percentdiff
