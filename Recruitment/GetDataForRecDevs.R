# This code is just for extracting time series in order to estimate rec residuals

load("~/Dropbox/Chapter 3 - Sardine Anchovy/All Datasets/RAM/RAM.RData")      #RAM
load("~/Dropbox/Chapter 3 - Sardine Anchovy/All Datasets/FAO/FAO.RData")      #FAO
load("~/Dropbox/Chapter 3 - Sardine Anchovy/All Datasets/Barange/Barange_mystocks.RData")   #barange_noNAs

colnames(barange_noNAs)[7] <- 'sp'

library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)


# RAM stocks -- European pilchard and Pacific sardine
RAM.new <- RAM
RAM.new$datasource <- "RAM"
europilchard <- subset(RAM.new,scientificname=="Sardina pilchardus")
# ** I calculated catches by multiplying ssb * fishing.mortality
# europilchard$totalcatch <- europilchard$ssb*europilchard$fishing.mortality
# min(europilchard$year)
# max(europilchard$year)
# mean(europilchard$rec)
# write.csv(europilchard,"europilchard.csv")

# Barange et al. stocks - 
Barange.new <-  data.frame(datasource=rep("Barange",times=nrow(barange_noNAs)),scientificname=rep(NA,times=nrow(barange_noNAs)),barange_noNAs)

Barange.sardines <- subset(Barange.new,sp=="Sardine")
barange.sardines <- unique(Barange.new$stock)

# Combine data from all the datasets
alldat <- rbind.fill(list(Barange.new,RAM.new))
alldat$year <- as.numeric(alldat$year)
sards <- subset(alldat, sp=="Sardine")
sardstocks <- unique(sards$stock)


JS <- subset(Barange.new,stock=="Japanese sardine") # No need to calculate catches for this one!

setwd("/Users/mcsiple/Dropbox/Chapter 4 - Harvest Control Rules/Datasets/Sardine")

pick_fill <- function(stockname,sourcename){
  df <- subset(alldat,stock==stockname & datasource==sourcename)
  if(nrow(df)==0) stop("no data for this stock in this dataset")
  df <- df %>%
    mutate(ssb=replace(ssb, is.na(ssb), -1),
           rec=replace(rec, is.na(rec), -1),
           landings=replace(landings, is.na(landings), -1)) %>% as.data.frame()
  write.csv(df,file=paste(stockname,sourcename,".csv",sep=""))
  print(head(df))
}



#setwd("/Users/mcsiple/Dropbox/Chapter 4 - Harvest Control Rules/Datasets/Sardine")
#write.csv(JS,"japanesesardine_Barange.csv")


# Sunset plots -----------------------------------------------------------
# These are the SSB/rec plots where lighter colors = later years. Should
# show whether SSB is a strong influence on recruitment
# test plot
sards$stock <- as.factor(sards$stock)
# get just RAM and Barange data
completes <- subset(sards,datasource %in% c("Barange","RAM"))
stockies=list()
nstocks=length(completes)
for(i in 1:nstocks){
  onestock <- subset(completes,stock==sardstocks[i])
  
  plot1 <- ggplot(onestock,aes(x=ssb,y=rec,colour=year)) + geom_point(size=2) + geom_path(size=0.75) +
    scale_colour_continuous(low = "blue",high="orange") +
    theme_bw() + ylab("Recruitment") + xlab("SSB") + ggtitle(paste(sardstocks[i])) + theme(legend.position="none")
  stockies[[i]] <- plot1
}
stockies[[11]] <- NULL
stockies[[7]] <- NULL
n <- length(stockies)
nCol <- floor(sqrt(n))
par=mfrow=c(1,1)
do.call("grid.arrange", c(stockies, ncol=nCol))


ep <- ggplot(europilchard,aes(x=ssb,y=rec,colour=year)) + geom_point(size=4) + geom_path(size=1) +
  scale_colour_continuous(low = "blue",high="orange") +
  theme_bw() + ylab("Recruitment") + xlab("SSB") 
#  + facet_wrap(~region) 
js <- ggplot(JS,aes(x=ssb,y=rec,colour=year)) + geom_point(size=4) + geom_path(size=1) +
  scale_colour_continuous(low = "blue",high="orange") +
  theme_bw() + ylab("Recruitment") + xlab("SSB") 

library(gridExtra)
stockies <- list()
stockies[[1]] <- ep
stockies[[2]] <- js
n <- length(stockies)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(stockies, ncol=nCol))
#grid.arrange(ep,js,ncol=2)
#scale_colour_brewer(type="seq",palette ="BuPu" )
