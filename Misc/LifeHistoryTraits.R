library(rfishbase)
library(ggplot2)
library(plyr)
library(dplyr)
# maturity
# length_weight
# length_frequency
# reproduction

length_freq(species_list = "Sardina pilchardus",limit = 200)
lw <- length_weight(species_list = "Sardina pilchardus",limit = 200)

plot(0:30,lw$a[1]*(0:30)^lw$b[1],type="l",
     col = rainbow(n=15)[1],
     ylim=c(0,400),xlab="Length (mm)",ylab="Weight")
for(i in 2:nrow(lw)){
  lines(0:30,lw$a[i]*(0:30)^lw$b[i],col=rainbow(n=15)[i])
}


species <- c("Sardina pilchardus","Sardinops sagax","Sardinella aurita")
mat.table <- maturity(species_list = species)
summ.table <- mat.table %>% group_by(sciname) %>% summarize(minagemat = min(AgeMatMin,na.rm=T),
                                                            maxagemat = max(AgeMatMin,na.rm=T))
par(mfrow=c(1,1))
ggplot(mat.table, aes(x=Locality, y=Lm)) + geom_bar(stat = "identity")
