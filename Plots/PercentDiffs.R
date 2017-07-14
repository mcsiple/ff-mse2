dat <- read.csv("/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Results/FF_MSE_summary.csv")

dat <- subset(dat, h==0.6)
remove <- c("scenario","M.type","h")
dat <- dat[,-which(colnames(dat) %in% remove)]
mdat <- melt(dat,id.vars=c("Type","obs.error.type","HCR"))


#dat2 <- dcast(mdat,variable+value ~ Type+obs.error.type+HCR+recruit.sd+recruit.rho)
dat2 <- dcast(mdat,Type+HCR+variable ~ obs.error.type)
#dat2 <- dcast(mdat,value ~ Type+variable+obs.error.type+HCR+recruit.sd+recruit.rho)
dat2$percentdiff <- ((dat2$Tim - dat2$AC) / dat2$Tim) * 100
dat3 <- subset(dat2,!variable %in% c("n.10yrclose","meanDepl","good4preds","very.bad4preds"))
dat3 <- subset(dat3,HCR !="constF")
dat3 <- mutate(dat3,HCR = recode(HCR, 'cfp' = 'Stability-favoring',
                            #'constF' = 'Constant F',
                            'C1' = 'C1',
                            'C2' = 'C2',
                            'C3' = 'C3',
                            'trend' = "Trend-based"))
library(ggplot2)
setwd("/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Figures")
pdf(file = "PercentDiffsErrors.pdf",width = 12,height=10,useDingbats = FALSE)
ggplot(dat3, aes(x=variable,y=percentdiff),) +geom_bar(stat = "identity") + geom_hline(yintercept=0)+facet_grid(Type~HCR) +theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +ylab("% change when peaks hard to anticipate") +xlab("Performance metric") + ylim(c(-100,100))
dev.off()
#fill=percentdiff
