# COMPARE RECRUITMENT TIME SERIES' SPECTRAL PROPERTIES TO THE OBSERVED ONES
# using calc. beta.fast to get beta of every time series generated using the ts function
# 
Type = "Anchovy"
library(gtools)
library(RColorBrewer)
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

Date <- "2018-02-13"
# File path for both Macs:
#basepath <- "/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/"
# File path for PC laptop:
basepath <- "C:/Users/Megsie Siple/Dropbox/Chapter4-HarvestControlRules/"
path <- paste(basepath,"Results/",Type,Date,"/", sep="")
files <- list.files(path=path)
rm <- grep(files,pattern = ".txt") # Don't load the text summary
files <- files[-rm]
files <- files[grep(files, pattern = ".RData")] # only load RData files
files <- mixedsort(files) # IMPORTANT: need this line to order in same order as scenario table!
setwd(path)
results <- sapply(files, function(x) mget(load(x)),simplify=TRUE) 

years.test <- ncol(results[[1]]$biomass.oneplus.true) # just count cols to know how many years there are
nyrs.to.use <- 100 # How many years you want to use to calculate all your metrics
calc.ind <- tail(1:years.test, nyrs.to.use) # Which years to calculate median depletion over (length = nyrs.to.use)


# Read table of beta and SD for stocks
if(Type == "Anchovy"){B.data.mat <- read.csv(paste(basepath,"Datasets/AnchHerring_Rec_Spectra.csv",sep=""))}
if(Type == "Sardine"){B.data.mat <- read.csv(paste(basepath,"Datasets/Sardine_Rec_Spectra.csv",sep=""))}
if(Type == "Menhaden"){B.data.mat <- read.csv(paste(basepath,"Datasets/Datasets/Menhaden_Rec_Spectra.csv",sep=""))}


# Here is the spectral properties of the recruitment estimates we have:
s <- B.data.mat[,c("X","Rec.Length","SD.R","Beta.r","lowerCI","upperCI","beta.SE")]
nstocks <- nrow(s)
empty.mat <- matrix(NA,nrow=20,ncol=ncol(s)) #This is the number of "scenarios" there are in the MSE. *** this may change
colnames(empty.mat) <- names(s)
s <- rbind(s,empty.mat)
levels(s$X) <- c(levels(s$X),"Simulated")
s[(nstocks+1):nrow(s),"X"] <- "Simulated"

# Plot estimated betas of the stocks that we have time series for:
par(mfrow=c(1,1))
plot(1:nrow(s),s[,'Beta.r'],pch=19,xaxt="n", 
     xlab="Stock",ylab="Estimated Beta",
     ylim=range(c(s[,'lowerCI'],s[,'upperCI']),na.rm = T))
axis(1, at=1:nrow(s), labels=rownames(s))
arrows(x0 = 1:nrow(s),y0=s[,'lowerCI'],
       x1 = 1:nrow(s),y1 = s[,'upperCI'],
       length=0.05, angle=90, code=3)

# Then, add some points from the simulated recruitment time series
palette <- rainbow(length(results),alpha = 0.4)

# Plot all the different scenarios on the same plot!
for(r in 1:length(results)){
  beta.vec <- vector()
  for(i in 1:50){
    beta.vec[i] <- calc.beta.fast(results[[r]]$rec[i,calc.ind]) # Calculate beta for each of 50 random recruitment time series taken from the simulations
  }
  points(x=jitter(rep(nstocks+r,times=length(beta.vec))),y=beta.vec,col= palette[r],pch=19)
}


# Next, just plot one (because in this case, they’re all really similar) --------
pdf(paste("BetaCheck_",Type,".pdf",sep=""),width=12,height=7,useDingbats = FALSE)
          par(mfrow=c(1,2),las=2,mar=c(15,5,1,3))
          plot(1:(nstocks+1),s[1:(nstocks+1),'Beta.r'],pch=19,xaxt="n", xlab="",ylab="Estimated Beta",ylim=range(c(s[,'lowerCI'],s[,'upperCI']),na.rm = T))
          axis(1, at=1:(nstocks+1), labels=s[1:(nstocks+1),"X"])
          arrows(x0 = 1:(nstocks+1),y0=s[,'lowerCI'],x1 = 1:(nstocks+1),y1 = s[,'upperCI'],length=0.05, angle=90, code=3)
          
          points(x=jitter(rep((nstocks+1),times=length(beta.vec))),y=beta.vec,col= 'lightblue',pch=19)
          
          # Same thing with SD
          plot(1:(nstocks+1),s[1:(nstocks+1),'SD.R'],pch=19,xaxt="n", xlab="",ylab="SD(Standardized Recruitment)",ylim=c(0,2))
          axis(1, at=1:(nstocks+1), labels=s[1:(nstocks+1),"X"])
          sd.vec <- vector()
          for(r in 1:length(results)){
            for(i in 1:50){
              stdized.rec <- results[[r]]$rec[i,calc.ind] / mean(results[[r]]$rec[i,calc.ind])
              sd.vec[i] <- sd(stdized.rec)
            }
            points(x=jitter(rep((nstocks+1),times=length(sd.vec))),y=sd.vec,col= 'lightblue',pch=19)
          }
dev.off()
