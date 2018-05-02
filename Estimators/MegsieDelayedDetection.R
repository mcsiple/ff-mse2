tau0 <-  0.1
sigma0 <- 0.2
tau1 <- (1/tau0^2+1/sigma0^2)^(-0.5)


amp <- 100
freq <- 15
# get a new random seed
rm(.Random.seed)

seed.index <- runif(1,1,100)
set.seed(seed.index)
timelist <- seq(1,60)
B <- amp + amp/1.3 * sin(pi*timelist/freq)
y <- rep(NA, length(timelist))
E <- rep (NA, length(timelist))
E[1]<-B[1]*exp(rnorm(1,0, sigma0))

for (t in 2:length(timelist)) {
  y[t] <- log (B[t]/E[t-1])
  mu1.tmp <- y[t]*(1-sigma0^2/(tau0^2+sigma0^2))
  E[t] <- E[t-1] * exp(rnorm(1,mu1.tmp,tau1))
}

setwd('/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Figures')
plotfilename <- "DelayedDetection_blackbkd.pdf"
pdf(file = plotfilename, height = 7, width = 10,useDingbats = FALSE)
par(mfrow=c(1,1))
lotau.col <- "#FB323B" # blues: "#0F2B5F" #"red" #539E59" 
hitau.col <-  "#C5E137" # blues: "#8EC1E7" #"orange" #BBD961" # light green
par(bg = 'black', fg = 'white') # set background to black, foreground white
plot(timelist,B, type = "l", lwd = 3, col = "white", xlab = "", ylab = " ", ylim = c(0, 275), axes = F)
box() #make col="black" for plotting in paper (or other white bkd) 
mtext(side = 2, text = "Population Size", line = 1)
points(timelist, B, pch=21, bg = "white") #change bg="black" if paper
lines(timelist, E, lwd = 3, col = lotau.col)
points(timelist, E, pch=21, bg = lotau.col, col=lotau.col)


tau0 <-  1000
sigma0 <- 0.2
tau1 <- (1/tau0^2+1/sigma0^2)^(-0.5)


set.seed(seed.index)
y <- rep(NA, length(timelist))
E <- rep (NA, length(timelist))
E[1]<-B[1]*exp(rnorm(1,0, sigma0))

for (t in 2:length(timelist)) {
  y[t] <- log (B[t]/E[t-1])
  mu1.tmp <- y[t]*(1-sigma0^2/(tau0^2+sigma0^2))
  E[t] <- E[t-1] * exp(rnorm(1,mu1.tmp,tau1))
}

lines(timelist, E, lwd = 3, col = hitau.col)
points(timelist, E, pch=21, bg = hitau.col,col = hitau.col)

legend('topright',pch=rep(21,times=3),pt.bg=c('white',lotau.col,hitau.col), # put 'black' here for paper
       col=c('white',lotau.col,hitau.col),
       legend=c("Biomass",
                expression(paste("Large changes not expected (",tau," = 0.1",")",sep="")),
                expression(paste("Large changes expected (",tau," = 1000",")",sep=""))))

dev.off()
#system2("open",args=c("-a Skim.app",plotfilename))
