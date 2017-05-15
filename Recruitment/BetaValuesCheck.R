# COMPARE RECRUITMENT TIME SERIES' SPECTRAL PROPERTIES TO THE OBSERVED ONES
# using calc. beta.fast to get beta of every time series generated using the ts function
# 


# Run these functions first -----------------------------------------------
## Add an alpha value to a colour (from Mages' blog, http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html)
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

# fast functon to calculate beta, does not use the multiple segment method
calc.beta.fast<-function(ts.data){
  
  tmp.spec<-spectrum(ts.data,plot=FALSE)
  freq<-tmp.spec$freq
  spec<-tmp.spec$spec
  tmp.data.frame<-data.frame(freq=log(freq),spec=log(spec))
  lm.out<-lm(spec~freq,data=tmp.data.frame)
  beta<--lm.out$coef[2]
  return(beta)
}


years.test <- ncol(results[[1]]$biomass) # just count cols to know how many years there are
nyrs.to.use <- 100 # How many years you want to use to calculate all your metrics
calc.ind <- tail(1:years.test, nyrs.to.use) # Which years to calculate median depletion over (length = nyrs.to.use)


# Read table of beta and SD for stocks
B.data.mat <- read.csv("/Users/mcsiple/Dropbox/Chapter 4 - Harvest Control Rules/Datasets/AnchHerring_Rec_Spectra.csv")

# Here is the spectral properties of the rec time series we have:
s <- B.data.mat[,c("X","Rec.Length","SD.R","Beta.r","lowerCI","upperCI","beta.SE")]
nscenarios <- nrow(s)
empty.mat <- matrix(NA,nrow=length(results),ncol=ncol(s))
colnames(empty.mat) <- names(s)
s <- rbind(s,empty.mat)
levels(s$X) <- c(levels(s$X),"Simulated")
s[(nscenarios+1):nrow(s),"X"] <- "Simulated"

# First, plot estimated betas of the stocks that we have time series for:
par(mfrow=c(1,1))
plot(1:nrow(s),s[,'Beta.r'],pch=19,xaxt="n", xlab="Stock",ylab="Estimated Beta",ylim=range(c(s[,'lowerCI'],s[,'upperCI']),na.rm = T))
axis(1, at=1:nrow(s), labels=rownames(s))
arrows(x0 = 1:nrow(s),y0=s[,'lowerCI'],x1 = 1:nscenarios,y1 = s[,'upperCI'],length=0.05, angle=90, code=3)

# Then, add some points from the simulated recruitment time series
palette <- add.alpha(c('#4ac7cb',
  '#e47647',
  '#57d190',
'#c173dd',
'#6ac350',
 '#d889c6',
  '#a6ba4c',
  '#799ae1',
  '#d7ac37',
  '#e66c85',
  '#76b37d',
 '#c5a061'),alpha = 0.4)

# Plot all the different scenarios on the same plot!
for(r in 1:length(results)){
  beta.vec <- vector()
  for(i in 1:100){
    beta.vec[i] <- calc.beta.fast(results[[r]]$rec[i,calc.ind])
  }
  points(x=jitter(rep(nscenarios+r,times=length(beta.vec))),y=beta.vec,col= palette[r],pch=19)
}


# Next, just plot one (because in this case, they’re all really similar) --------
par(mfrow=c(1,2),las=2,mar=c(15,5,1,3))
plot(1:(nscenarios+1),s[1:(nscenarios+1),'Beta.r'],pch=19,xaxt="n", xlab="",ylab="Estimated Beta",ylim=range(c(s[,'lowerCI'],s[,'upperCI']),na.rm = T))
axis(1, at=1:(nscenarios+1), labels=s[1:(nscenarios+1),"X"])
arrows(x0 = 1:(nscenarios+1),y0=s[,'lowerCI'],x1 = 1:(nscenarios+1),y1 = s[,'upperCI'],length=0.05, angle=90, code=3)

points(x=jitter(rep((nscenarios+1),times=length(beta.vec))),y=beta.vec,col= 'lightblue',pch=19)

# Same thing with SD
plot(1:(nscenarios+1),s[1:(nscenarios+1),'SD.R'],pch=19,xaxt="n", xlab="",ylab="SD(Standardized Recruitment)",ylim=c(0,2))
axis(1, at=1:(nscenarios+1), labels=s[1:(nscenarios+1),"X"])
sd.vec <- vector()
for(r in 1:length(results)){
  for(i in 1:100){
    stdized.rec <- results[[r]]$rec[i,calc.ind] / mean(results[[r]]$rec[i,calc.ind])
    sd.vec[i] <- sd(stdized.rec)
  }
  points(x=jitter(rep((nscenarios+1),times=length(sd.vec))),y=sd.vec,col= 'lightblue',pch=19)
}

