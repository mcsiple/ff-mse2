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
# B <- c(247728.6, 243145.9, 249756.3, 270354.6, 299320.8, 315663.6, 
#        326716.2, 332712.2, 344883.1, 353495.8, 375949.5, 394098.1, 428431, 
#        444125.4, 436833.9, 482789.8, 524401.2, 565965.3, 615554.8, 644774.2, 
#        673944.3, 704280.4, 720881.7, 729185.4, 729201.7, 787458.8, 820629.6, 
#        811141.4, 787299.2, 764825.5, 752510.2, 731537.9, 690350.8, 671351.9, 
#        651918.2, 619961.5, 585162.9, 551323.8, 541180, 540778.3, 559508.9, 
#        569509.9, 579629.6, 580662.4, 579789.7, 562445.6, 526011.4, 518258, 
#        531159.7, 525685.4, 507961, 475416.6, 474533.3, 511293.3, 546476.8, 
#        589062.1, 615273, 621396.5, 605216.1, 574825.6, 571952.5, 574753.7, 
#        581205.3, 592003.9, 631733.3, 657159.7, 652735.8, 675975.8, 687891.6, 
#        677807.6, 641046.1, 585274.8, 525519.2, 479589.5, 460408.3, 465499.1, 
#        471731.3, 479816.2, 476590.2, 461649.8, 458247.4, 444747, 463652.2, 
#        486741.3, 512796.1, 521745.9, 515376.8, 488594.1, 457759.5, 435229.8, 
#        423622.6, 408471.6, 387875.4, 365308, 364256.1, 357080.6, 349669.6, 
#        372125.1, 399010.3, 405034.5, 397760.3)[1:60]
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
lotau.col <- "#bdbdbd"  
hitau.col <-  "#636363"
#par(bg = 'black', fg = 'white') # set background to black, foreground white - for presentations
plot(timelist,B, type = "l", lwd = 3, col = "black", xlab = "", ylab = " ",ylim=range(c(B,E)), axes = F)
box() #make col= "black" for plotting in paper (or other white bkd) 
mtext(side = 2, text = "Population Size", line = 1)
mtext(side = 1, text = "Time",line = 1)
points(timelist, B, pch=21, bg = "black") #change bg="black" if paper
lines(timelist, E, lwd = 3, col = lotau.col)
points(timelist, E, pch=21, bg = lotau.col, col=lotau.col)


tau0 <-  100
sigma0 <- 0.1
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
