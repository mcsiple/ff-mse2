# Functions for summarizing and plotting performance measures
# Functions for calculating performance measures ---------------------------
library(RcppRoll)

cv <- function(ts){
  mean = mean(ts)
  sd = sd(ts)
  CV <- sd/mean
  return(CV)
}

nzeroes <- function(x){         # How many years with zero catch
  #x is a vector
  n <- length(which(x==0))
  return(n)
}

good4pred <- function(x, F0.x){ # Nyears above a "good for predators" threshold (here, it's 80% of the long term mean unfished biomass, given rec variation). This is the same as "number of bonanza years"
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

duration.collapse <- function(x, F0.x){
  # This function returns the min, max, and mean duration of collapses
  collapses <- collapse.index(x = x, F0.x = F0.x)
  collapse_tf <- logical(length = length(F0.x))
  collapse_tf[collapses] <- TRUE # change to True/False vector so you can measure duration of collapses
  y <- rle(collapse_tf)
  collapse.lengths <- y$lengths[y$values==TRUE]
  
  if(length(collapse.lengths) == 0){
    shortest.collapse <- NA
    longest.collapse <- NA
    avg.collapse.length <- NA
  }  else{
    shortest.collapse <- min(collapse.lengths, na.rm = T)
    longest.collapse <- max(collapse.lengths, na.rm = T)
    avg.collapse.length <- mean(collapse.lengths, na.rm = T)
  }
  return(list("shortest.collapse" = shortest.collapse,
              "longest.collapse" = longest.collapse,
              "avg.collapse.length" = avg.collapse.length))
}

duration.bonanza <- function(x, F0.x){
  # This function returns the min, max, and mean duration of bonanzas
    bonanzas <- bonanza.index(x = x, F0.x = F0.x)
    bonanza_tf <- logical(length = length(F0.x))
    bonanza_tf[bonanzas] <- TRUE # change to True/False vector so you can measure duration of bonanzas
    y <- rle(bonanza_tf)
    bonanza.lengths <- y$lengths[y$values==TRUE]
    if(length(bonanza.lengths) == 0){
      shortest.bonanza <- NA
      longest.bonanza <- NA
      avg.bonanza.length <- NA
    }  else{
      shortest.bonanza <- min(bonanza.lengths, na.rm = T)
      longest.bonanza <- max(bonanza.lengths, na.rm = T)
      avg.bonanza.length <- mean(bonanza.lengths, na.rm = T)
    }
      return(list("shortest.bonanza" = shortest.bonanza,
                  "longest.bonanza" = longest.bonanza,
                  "avg.bonanza.length" = avg.bonanza.length))
}

#mean.duration.bonanza(x=testie$biomass.total.true,F0.x = F0.Type)
n.bonafide.collapses <- function(ts.length, coll.yrs){
  binary.coll.vec <- rep(0, times=ts.length)
  binary.coll.vec[coll.yrs] <- 1 # Turn into vector of 0s and 1s
  m <- rle(binary.coll.vec==1)
  c <- m$lengths
  d <- m$values
  counter <- 0
  if(any(d)){
  starts <- which(d) # indexes rle outputs to get locations where collapses start
  for(l in 1:length(starts)){
    ind.coll <- starts[l]
    if(c[ind.coll]>=2 & ind.coll != 1 & (ind.coll+2) < length(c)){
      #if the sequence of collapse is preceded and followed by at least 4 yrs non-collapse 
      if(c[(ind.coll-1)]>=4 & c[ind.coll+1]>=4){counter <- counter+1} else{ 
        if(c[ind.coll+1]<=4 & c[ind.coll+2]>=2){counter <- counter +1} else{
        counter=counter}}}else{counter=counter}
    # if the 4 years leading up to the collapse are not-collapsed, and same with the four years after, it's a "bonafide collapse"
  }
  # Last case: if the final years are collapsed, need to count those!
  if(tail(d,1) & tail(c,1) > 4){counter=counter+1}else{counter=counter}
  }
  
  return(counter)
}


n.multiyr.closures <- function(x, threshold = NA) { #where x is a matrix, rows are sims, cols are years
  count1 <- count5 <- count10 <- vector(length=nrow(x))
  for(i in 1:nrow(x)){ #Either catch OR biomass
    ltm <- mean(x[i,],na.rm = T)
    if(is.na(threshold)){ thresh <- 0.01*ltm } # threshold can be anything - default is 1% of long term mean C or B
    else{thresh = threshold}
    badTorF <- x[i,] <= thresh
    oneyr <- sum(roll_sum(badTorF, 1) == 1,na.rm = T)
    fiveyr <- sum(roll_sum(badTorF, 5) == 5,na.rm = T)
    tenyr <- sum(roll_sum(badTorF, 10) == 10,na.rm = T)
    count1[i] <- oneyr # number of one-year closures
    count5[i] <- fiveyr #number of five year closures
    count10[i] <- tenyr # number of 10 year closures
  }
  
  return(list(count1 = count1,count5 = count5,count10 = tenyr)) # Mean number of 5- and 10-yr closures
}

# Rcpproll demo:
# set.seed(1); x <- sample(c(T, F), 100, replace = T); sum(RcppRoll::roll_sum(x, 3) == 3)

# For testing summ.tab - run these then you can run the interior of the function:
# result.list = results[[14]] #C1=14; C2=20
# individual.sim = FALSE
# calc.ind = 151:250
# ou.ind = NA
# collapse.severity.C1[which(!is.na(collapse.severity.C1))]
# collapse.severity.C2[which(!is.na(collapse.severity.C2))]


# Get summary metrics so you can make Zeh plots (or whatever you want)
summ.tab <- function(result.list, individual.sim = FALSE,calc.ind = calc.ind, ou.ind = NA){ #result.list is one of the results (=1 harvest rule, 1000 time series of biomass, catch, fishing, rec, depl)
  # calc.ind is the index of which years to calculate over
  for(i in 1:length(result.list)){
    result.list[[i]] <- result.list[[i]][,calc.ind]
    if(length(ou.ind)>1){result.list[[i]] <- result.list[[i]][ou.ind,]}
    }  # Trim results to the years we're using

  LTmeans.list <- lapply(result.list,FUN = rowMeans,na.rm=TRUE) 
  # median and quantiles of LTM of all PMs
  ltm <- lapply(LTmeans.list,FUN = quantile, probs = c(0.05,0.5,0.95),na.rm=TRUE) 
  # mean nonzero catch
  catch <- result.list$total.catch
  nonzero.catch <- ifelse(catch<0.1,NA,catch)
  ltm.nzc1 <- rowMeans(nonzero.catch,na.rm=TRUE)
  ltm.nzc2 <- quantile(ltm.nzc1,probs = c(0.05,0.5,0.95),na.rm = TRUE)
  # SDcatch
  sd.catches <- apply(catch, MARGIN = 1,FUN = sd, na.rm = TRUE)
  sd.catch <- quantile(sd.catches,probs = c(0.05,0.5,0.95))
  
  #Number of years w zero catch
  nz1 <- apply(catch,MARGIN = 1,FUN = nzeroes)
  
  #5- and 10-yr closures
  n.5yrclose <- n.multiyr.closures(catch,threshold = 0)$count5 / nz1          # P(5yr closure | closure)
  n.10yrclose <- n.multiyr.closures(catch,threshold = 0)$count10 / n.5yrclose # P(10yr closure | 5yr closure) - this one is depracated
  
  # SDbiomass
  true.biomass <- result.list$biomass.total.true
  sd.Bs <- apply(true.biomass, MARGIN = 1,FUN = sd, na.rm = TRUE)
  sd.B <- quantile(sd.Bs,probs = c(0.05,0.5,0.95))
  
  #Years that are "good for predators"
  g4p.vec <- vector()
  for(sim in 1:nrow(true.biomass)){
    g4p.vec[sim] <- good4pred(x = true.biomass[sim,],F0.x = result.list$no.fishing.tb[sim])
  }

  # Number of years that are below a predator threshold
  yrs.bad <- vector()
  for(sim in 1:nrow(true.biomass)){
    yrs.bad[sim] <- bad4pred(x = true.biomass[sim,],F0.x = result.list$no.fishing.tb[sim,])
  }
 
  # Maximum duration of collapse
  max.duration.collapse = min.duration.collapse = avg.duration.collapse <- vector()
  max.duration.bonanza = min.duration.bonanza = avg.duration.bonanza <- vector()
  prob.collapse = collapse.severity = cv.vec = bonafide.collapse = vector()
  depl <- matrix(nrow=nrow(true.biomass),ncol=ncol(true.biomass))
  for(sim in 1:nrow(true.biomass)){ # Each one = 1 simulation
    max.duration.collapse[sim] <- duration.collapse(x = true.biomass[sim,],F0.x = result.list$no.fishing.tb[sim,])$longest.collapse
    min.duration.collapse[sim] <- duration.collapse(x = true.biomass[sim,],F0.x = result.list$no.fishing.tb[sim,])$shortest.collapse
    avg.duration.collapse[sim] <- duration.collapse(x = true.biomass[sim,],F0.x = result.list$no.fishing.tb[sim,])$avg.collapse.length
    
    max.duration.bonanza[sim] <- duration.bonanza(x = true.biomass[sim,],F0.x = result.list$no.fishing.tb[sim,])$longest.bonanza
    min.duration.bonanza[sim] <- duration.bonanza(x = true.biomass[sim,],F0.x = result.list$no.fishing.tb[sim,])$shortest.bonanza
    avg.duration.bonanza[sim] <- duration.bonanza(x = true.biomass[sim,],F0.x = result.list$no.fishing.tb[sim,])$avg.bonanza.length
    
    # Probability of collapse
    coll.yrs <- collapse.index(x = true.biomass[sim,],F0.x = result.list$no.fishing.tb[sim,]) 
    prob.collapse[sim] <- length(coll.yrs) / length(true.biomass[sim,])
    
    # Severity of collapse (as a percentage of mean unfished biomass)
    coll.thresh <- mean(result.list$no.fishing.tb[sim,]) * 0.2
    collapse.yrs <- true.biomass[sim,coll.yrs]
    collapse.severity[sim] <- 1 - (mean(collapse.yrs,na.rm=T) / coll.thresh ) 
    # median biomass at collapse, over the median unfished biomass for that simulation - a higher value of that is "better" (less severe) so substract it from 1 to get higher number == worse severity.
    
    #Depletion based on unfished B0
      depl[sim,] <-  true.biomass[sim,]/result.list$no.fishing.tb[sim,]
      
      #CV of catches
      cv.vec[sim] <- cv(catch[sim,])
      
      # N "bonafide" collapse periods (this is different than dipping below 0.2Bbar, you have to have a streak of non-collapse years before and after)
      bonafide.collapse[sim] <- n.bonafide.collapses(ts.length = length(true.biomass[sim,]),coll.yrs = coll.yrs)           
                  # 
    } #end of sim loop
  
  mean.depl <- apply(depl,MARGIN = 1,FUN = mean)
  
  # Performance metrics
  interval <- c(0.05,0.5,0.95)
  
  ltm.c <- ltm$total.catch  # ltmcatch
  ltm.nzc <- ltm.nzc2       #ltmnonzerocatch
  SDcatch <- sd.catch       #SD(Catch)
  n5yr <- quantile(n.5yrclose,probs = interval, na.rm=TRUE)     #n.5yrclose
  n10yr <- quantile(n.10yrclose,probs = interval, na.rm=TRUE)   #n.10yrclose
  nz <- quantile(nz1,probs = interval)              #nyrs.0catch
  ltm.b <- ltm$biomass.total.true                   #LTMBiomass
  g4p <- quantile(g4p.vec,probs = interval)         #Nyears "good for predator"
  sdB <- sd.B                               #SD(Biomass)
  b4p <- quantile(yrs.bad,probs = interval) #p(bad4preds)
  ltm.depl <- quantile(depl, probs = interval,na.rm=T)                     # Mean depletion
  prob.coll <- quantile(prob.collapse,probs = interval,na.rm=T)
  severity <- quantile(collapse.severity,probs = interval,na.rm=T)
  cv.catch <- quantile(cv.vec,probs = interval, na.rm=T)
  b.collapse <- quantile(bonafide.collapse,probs = interval, na.rm=T)
  # Awkward but necessary for when there are NAs in these vectors
  if(all(is.na(max.duration.collapse))){
    overall.max.coll.len <- rep(NA,times=3)
      }else{
  overall.max.coll.len <-  c(NA, max(max.duration.collapse,na.rm=T), NA)  #Fill quantiles w NAs because other metrics have quantiles and these don't!
      }
  
  # Ugh
  if(all(is.na(max.duration.bonanza))){
    overall.max.bon.len <- rep(NA,times=3)
  }else{
    overall.max.bon.len <- c(NA, max(max.duration.bonanza,na.rm=T), NA)
  }
  
  # Sigh 
  if(all(is.na(avg.duration.bonanza))){
    bon.length <- rep(NA,times=3)
  } else{
  bon.length <- quantile(avg.duration.bonanza,probs = interval,na.rm = T)
  }
  
  # whatever
  if(all(is.na(avg.duration.collapse))){
    coll.length <- rep(NA,times=3)
  } else{
    coll.length <- quantile(avg.duration.collapse,probs = interval, na.rm = T)
  }
  
  output <- data.frame(PM = performance.measures, loCI = NA, med = NA, hiCI = NA)
  output[,-1] <- rbind(ltm.c,ltm.nzc,SDcatch,n5yr,n10yr,nz,ltm.b,g4p,sdB,b4p,ltm.depl,overall.max.coll.len,overall.max.bon.len,bon.length,coll.length,prob.coll,severity,cv.catch,b.collapse)
  sim.output <- list()
  #LTmean catch
  #LTmean nonzero catch
  #SDcatch
  #n5yr
  #n10yr
  #yrs w zero catch
  #LTmean biomass
  #yrs good4preds
  #SDbiomass
  #yrs bad4preds
  #ltm depletion
  #max collapse length (across ALL SIMULATIONS) - do not use in correlations
  #max bonanza length (across ALL SIMULATIONS) - do not use in correlations
  #mean bonanza length
  #mean collapse length
  #probability of a collapse
  #collapse severity
  #cv(catch)
  #sustained collapse
  sim.output[[1]] <- LTmeans.list$total.catch
  sim.output[[2]] <- ltm.nzc1 
  sim.output[[3]] <- sd.catches
  sim.output[[4]] <- n.5yrclose
  sim.output[[5]] <- n.10yrclose
  sim.output[[6]] <- nz1
  sim.output[[7]] <- LTmeans.list$biomass.total.true
  sim.output[[8]] <- g4p.vec
  sim.output[[9]] <- sd.Bs
  sim.output[[10]] <- yrs.bad
  sim.output[[11]] <- LTmeans.list$depl
  sim.output[[12]] <- avg.duration.bonanza
  sim.output[[13]] <- avg.duration.collapse
  sim.output[[14]] <- prob.coll
  sim.output[[15]] <- severity
  sim.output[[16]] <- cv.catch
  sim.output[[17]] <- b.collapse
  if(individual.sim==TRUE){return(sim.output)}
  else{return(output)}
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


# Quickly plot time series for a given scenario and simulation # ----------

plot.scenario <- function(result.ind, nyrs.to.use = 100, sim = 1){
  # quickly plot basic stuff for one scenario (combo of LH traits and control rule)
  to.use <- results[[result.ind]]
  nsims <- nrow(to.use$biomass.oneplus.true) # just count nrows to know how many sims there are
  years.test <- ncol(to.use$biomass.oneplus.true) # just count cols to know how many years there are
  calc.ind <- tail(1:years.test, nyrs.to.use) # Which years to calculate median depletion over (length = nyrs.to.use)
  for(i in 1:length(to.use)){ # Trim results to yrs calculating metrics for
    to.use[[i]] <- to.use[[i]][,calc.ind]
  } 
  mtu <- melt(to.use)
  colnames(mtu) <- c("Sim","Year","value","variable")
  toplot <- subset(mtu, Sim==sim)
  ggplot(toplot,aes(x=Year,y=value)) + geom_line() + facet_wrap(~variable,scales="free_y")
}

