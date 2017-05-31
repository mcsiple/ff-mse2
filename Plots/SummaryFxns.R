# Functions for summarizing and plotting performance measures
# Functions for calculating performance measures ---------------------------
library(RcppRoll)

nzeroes <- function(x){ # How many years with zero catch
  #x is a vector
  n <- length(which(x==0))
  return(n)
}

good4pred <- function(x, F0.x){ # Nyears above a "good for predators" threshold (here, it's 75% of the long term mean unfished biomass, given rec variation)
  B.bar <- mean(F0.x)
  thresh <- 0.8*B.bar
  g4p <- length(which(x>thresh))
  return(g4p)
}

bad4pred <- function(x, F0.x){ # Number of years that are below a certain threshold (here, it's 10% of the long term mean unfished biomass, given rec variation)
  # x is a time series of biomass
  B.bar <- mean(F0.x)
  thresh <- 0.2*B.bar
  b4p <- length(which(x<thresh))
  return(b4p)
}

collapse.index <- function(x,F0.x){ # Which years have "collapses"-- use this function to ID the collapse years and look at fishing rates leading up to those collapses.
  B.bar <- mean(F0.x)
  thresh <- 0.2*B.bar
  if(any(x<thresh)){
    collapse.ind <- which(x<thresh)
  }
  else collapse.ind <- 0
  return(collapse.ind)
}

bonanza.index <- function(x,F0.x){ # Which years have "bonanzas"-- use this function to ID the bonanza years and what conditions led up to them.
  B.bar <- mean(F0.x)
  thresh <- 0.8*B.bar
  if(any(x>thresh)){
    bonanza.ind <- which(x>thresh)
  }
  else bonanza.ind <- 0
  return(bonanza.ind)
}

mean.duration.bc <- function(x, F0.x){
  collapses <- collapse.index(x = x, F0.x = F0.x)
  collapse_tf <- vector(FALSE, length = length(F0.x))
  collapse_tf[collapses] <- TRUE # change to True/False vector so you can measure duration of collapses
  y <- rle(collapse_tf)
  collapse.lengths <- y$lengths[y$values==TRUE]
  avg.collapse.length <- mean(collapse.lengths)
  return(avg.collapse.length)
}

mean.duration.bonanza <- function(x, F0.x){
    bonzanzas <- bonanza.index(x = x, F0.x = F0.x)
    bonanza_tf <- vector(FALSE, length = length(F0.x))
    bonanza_tf[bonzanzas] <- TRUE # change to True/False vector so you can measure duration of bonanzas
    y <- rle(bonanza_tf)
    bonanza.lengths <- y$lengths[y$values==TRUE]
    avg.bonanza.length <- mean(bonanza.lengths)
    return(avg.bonanza.length)
  }
n.multiyr.closures <- function(x, threshold = NA) { #where x is a matrix, rows are sims, cols are years
  count5 <- count10 <- vector(length=nrow(x))
  for(i in 1:nrow(x)){ #Either catch OR biomass
    ltm <- mean(x[i,])
    if(is.na(threshold)){ thresh <- 0.01*ltm } # threshold can be anything - CHANGE THIS TO BE DIFF FOR DIFF HCRS!!!
    else{thresh = threshold}
    badTorF <- x[i,] <= thresh
    fiveyr <- sum(roll_sum(badTorF, 5) == 5)
    tenyr <- sum(roll_sum(badTorF, 10) == 10)
    count5[i] <- fiveyr #number of five year closures
    count10[i] <- tenyr # number of 10 year closures
  }
  
  return(list(count5 = count5,count10 = tenyr)) # Mean number of 5- and 10-yr closures
}

# Rcpproll demo:
# set.seed(1); x <- sample(c(T, F), 100, replace = T); sum(RcppRoll::roll_sum(x, 3) == 3)


# Get summary metrics so you can make Zeh plots (or whatever you want)
summ.tab <- function(result.list){ #result.list is one of the results (=1 harvest rule, 1000 time series of biomass, catch, fishing, rec, depl)
  for(i in 1:length(result.list)){
    result.list[[i]] <- result.list[[i]][,calc.ind]
  } # Trim results to the years we're using
  performance.measures
  LTmeans.list <- lapply(result.list,FUN = rowMeans) 
  # median and quantiles of LTM of all PMs
  ltm <- lapply(LTmeans.list,FUN = quantile, probs = c(0.05,0.5,0.95)) 
  # mean nonzero catch
  catch <- result.list$total.catch
  nonzero.catch <- ifelse(catch<0.1,NA,catch)
  ltm.nzc1 <- rowMeans(nonzero.catch,na.rm=TRUE)
  ltm.nzc2 <- quantile(ltm.nzc1,probs = c(0.05,0.5,0.95),na.rm = TRUE)
  # SDcatch
  sd.catches <- apply(catch, MARGIN = 1,FUN = sd, na.rm = TRUE)
  sd.catch <- quantile(sd.catches,probs = c(0.05,0.5,0.95))
  
  #5- and 10-yr closures
  n.5yrclose <- n.multiyr.closures(catch)$count5
  n.10yrclose <- n.multiyr.closures(catch)$count10
  
  #Number of years w zero catch
  nz1 <- apply(catch,MARGIN = 1,FUN = nzeroes)
  
  # SDbiomass
  true.biomass <- result.list$biomass.total.true
  sd.Bs <- apply(true.biomass, MARGIN = 1,FUN = sd, na.rm = TRUE)
  sd.B <- quantile(sd.Bs,probs = c(0.05,0.5,0.95))
  
  #Years that are "good for predators"
  g4p.vec <- vector()
  for(i in 1:nrow(true.biomass)){
    g4p.vec[i] <- good4pred(x = true.biomass[i,],F0.x = result.list$no.fishing.tb[i,])
  }

  # Number of years that are below a predator threshold
  yrs.bad <- vector()
  for(i in 1:nrow(true.biomass)){
    yrs.bad[i] <- bad4pred(x = true.biomass[i,],F0.x = result.list$no.fishing.tb[i,])
  }
 
  
  # Performance metrics
  interval <- c(0.05,0.5,0.95)
  
  ltm.c <- ltm$total.catch  # ltmcatch
  ltm.nzc <- ltm.nzc2       #ltmnonzerocatch
  SDcatch <- sd.catch       #SD(Catch)
  n5yr <- quantile(n.5yrclose,probs = interval)     #n.5yrclose
  n10yr <- quantile(n.10yrclose,probs = interval)   #n.10yrclose
  nz <- quantile(nz1,probs = interval)              #nyrs.0catch
  ltm.b <- ltm$biomass.total.true      #LTMBiomass
  g4p <- quantile(g4p.vec,probs = interval)         #Nyears "good for predator"
  sdB <- sd.B                               #SD(Biomass)
  b4p <- quantile(yrs.bad,probs = interval) #p(bad4preds)
  ltm.depl <- ltm$depl                      # Mean depletion
  
  output <- data.frame(PM = performance.measures, loCI = NA, med = NA, hiCI = NA)
  output[,-1] <- rbind(ltm.c,ltm.nzc,SDcatch,n5yr,n10yr,nz,ltm.b,g4p,sdB,b4p,ltm.depl)
  return(output)
}

# Function to plot medians and certainty intervals from simulations:
plotintervals <- function(result.mat,ylab){ #result.mat is a matrix (e.g., biomass for results[[1]])
  median.vec <- apply(result.mat,MARGIN = 2,FUN = median)
  ints <- apply(result.mat,MARGIN = 2,FUN = quantile, probs = c(0.025,0.25,0.75,0.975))
  lo95 <- ints[1,calc.ind]
  hi95 <- ints[4,calc.ind]
  lo50 <- ints[2,calc.ind]
  hi50 <- ints[3,calc.ind]
  
  plot(1:nyrs.to.use,median.vec[calc.ind],
       type='l',lty=1,lwd=2,
       ylim=range(lo95,hi95),
       ylab = ylab,
       xlab = "Year")
  
  zz <- c(1:nyrs.to.use,tail(nyrs.to.use,n=1),rev(1:nyrs.to.use)) # for polygons
  aa <- c(hi95,0,rev(lo95))
  bb <-  c(hi50,0,rev(lo50))
  polygon(zz,aa,col=adjustcolor( "black", alpha.f = 0.2),border="NA")
  polygon(zz,bb,col=adjustcolor( "black", alpha.f = 0.2),border="NA")
}

## Add an alpha value to a colour (from Mages' blog, http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

