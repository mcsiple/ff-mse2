library(rfishbase)
library(ggplot2)
library(plyr)
library(dplyr)
# maturity
# length_weight
# length_frequency
# reproduction

length_freq(species_list = "Sardina pilchardus",limit = 200)
length_weight(species_list = "Sardina pilchardus",limit = 200)
species <- c("Sardina pilchardus","Sardinops sagax","Sardinella aurita")
mat.table <- maturity(species_list = species)
summ.table <- mat.table %>% group_by(sciname) %>% summarize(minagemat = min(AgeMatMin,na.rm=T),
                                                            maxagemat = max(AgeMatMax,na.rm=T))
par(mfrow=c(1,1))
ggplot(mat.table, aes(x=Locality, y=Lm)) + geom_bar(stat = "identity")
