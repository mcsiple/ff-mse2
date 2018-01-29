
# Generate recruitment deviations ----------------------------------------------------
# This is used in Thorson et al. 2014, CJFAS, Eq. 16
# http://www.nrcresearchpress.com/doi/full/10.1139/cjfas-2013-0645#.Wmo95JM-dR4
# The CV function is included just for plotting

toplot = FALSE

cv <- function(ts){
  mean = mean(ts)
  sd = sd(ts)
  CV <- sd/mean
  return(CV)
}

generate.devs <- function(N, rho, sd.devs, plot=FALSE,method=1){
  #' @description Function to generate recruitment deviations using the "time series" method 
  #' @param N - number of years for which to generate recruitment deviations 
  #' @param rho - autocorrelation
  #' @param sd.devs - standard deviation of the recruitment devs
  # This used to have a burn-in, but I checked and it doesn't seem necessary
  
  #Generate noise
  innov <- rnorm(N,mean = 0,sd = sd.devs) # For bias correction, can use mean = -sd.devs^2)/2 (I didn't see a big difference in devs, so have left mean = 0)

  #Make a vector to hold simulated deviations
  dev.ts <- numeric(length=N)
  dev.ts[1] <- innov[1]
  
  for(yr in 2:(N)){
    dev.ts[yr] <- rho * dev.ts[yr-1] + sqrt(1-rho^2) * innov[yr] # sqrt(1-rho^2) is a correction so it doesn't just wander unidirectionally
  }
  
  if(method==1){
    dev.ts <- exp(dev.ts)         # Since these are log deviations
    dev.ts <- dev.ts/mean(dev.ts) # Standardize so that mean is 1
  } 
    if(method==2){
    dev.ts <- dev.ts - mean(dev.ts) - 0.5 * var(dev.ts)
    dev.ts <- exp(dev.ts)
    }

  if(plot==TRUE){ plot(1:N, dev.ts, type='l',yaxt="n",ylab='',xlab='Year')
                  title(ylab="Recruit deviations", line=0.4, cex.lab=1.2)
                  #CV <- round(cv(dev.ts),digits = 2)
                  #text(x = 46,y = 0.99*max(dev.ts),labels = paste("CV =",CV))
                  }
  return(dev.ts)
}


# Test and plot the devs! -------------------------------------------------
              if(toplot == TRUE){
                
                set.seed(123)
                x1 <- generate.devs(N = 1000,rho = 0.9,sd.devs = .3,plot = FALSE,method=1)
                set.seed(123)
                x2 <- generate.devs(N = 1000,rho = 0.9,sd.devs = .3,plot = FALSE,method=2)
                plot(1:1000,x1,type='l')
                lines(1:1000,x2,col='red')
                      par(mfcol=c(3,3))
                      # "Sardine"
                      for(i in 1:3){
                        generate.devs(N=50,rho = 0.9,sd.devs = 0.1,plot = TRUE)
                        if(i==1){mtext(side = 3,text = "Sardine-like",outer = FALSE)}
                        }

                      # "Anchovy"
                      for(i in 1:3){
                        generate.devs(N=50,rho = 0.5,sd.devs = 0.3,plot = TRUE)
                        if(i==1){mtext(side = 3,text = "Anchovy-like",outer = FALSE)}
                      }

                      # "Menhaden"
                      for(i in 1:3){
                        generate.devs(N=50,rho = 0,sd.devs = 0.7,plot = TRUE)
                        if(i==1){mtext(side = 3,text = "Menhaden-like",outer = FALSE)}
                      }

              # Plot rec devs for presentation ------------------------------------------

                    sardine <- generate.devs(N=50,rho = 0.9,sd.devs = 0.1,plot = F)
                    anchovy <- generate.devs(N=50,rho = 0.5,sd.devs = 0.3,plot = F)
                    menhaden <- generate.devs(N=50,rho = 0,sd.devs = 0.7,plot = F)
                    df <- data.frame(year=1:length(devs), s = sardine, a = anchovy, m = menhaden)
                    cols <- c("#0C66AC", "#FCEA1B", "#E74690")
                    mdf <- melt(df, id.vars="year")
                    pdf("~/Dropbox/ChapterX-synthesis/RecruitmentDevExamples.pdf",width=10,height=3,useDingbats = F)
                    ggplot(mdf,aes(x=year,y=value,colour=variable)) + 
                      geom_line(lwd=2) + 
                      scale_colour_manual(values = cols) +
                      theme_black(base_size=14) + 
                      ylab("") + xlab("") +
                      facet_wrap(~variable,scales = "free_y") +
                      theme(axis.text.x = element_blank(),
                      axis.text.y = element_blank(),legend.position = "none",strip.background = element_blank(),
                      strip.text.x = element_blank()) # order is: sardine, anchovy, menhaden
                    dev.off()
                
                # Now plot unfished biomass with the rec devs -- ****fix this part to show biomass over time
                        # Load all the contents of the control rules folder:
                        setwd("/Users/mcsiple/Dropbox/Chapter 4 - Harvest Control Rules/Code/Control Rules")
                        lapply(list.files(pattern = "[.]R$", recursive = TRUE), source)
                        par(mfcol=c(3,3))


                      # "Sardine"
                      for(i in 1:3){
                        devs <- generate.devs(N=200,rho = 0.9,sd.devs = 0.1,plot = FALSE)
                        biomass <- calc.trajectory(lh = lh.test, obs.cv = NULL, init = init.test, rec.dev = devs, F0, cr, years = years.test, hcr.type = "hockeystick",obs.type="LN",const.f.rate = 0.2,equilib = NULL, buffer = NULL,steepness = 0.9, R0.traj = 0,tim.params=NULL)$biomass

                        if(i==1){mtext(side = 3,text = "Sardine-like",outer = FALSE)}
                      }

                      # "Anchovy"
                      for(i in 1:3){
                        generate.devs(N=200,rho = 0.5,sd.devs = 0.3,plot = FALSE)
                        if(i==1){mtext(side = 3,text = "Herring/anchovy-like",outer = FALSE)}
                      }

                      # "Menhaden"
                      for(i in 1:3){
                        generate.devs(N=200,rho = 0,sd.devs = 0.7,plot = FALSE)
                        if(i==1){mtext(side = 3,text = "Menhaden-like",outer = FALSE)}
                      }
              }



##################################################################
##################################################################

# pdf("~/Dropbox/Chapter4-HarvestControlRules/Figures/RecruitmentDevsExample.pdf",height = 3, width = 10)
# par(mfcol=c(1,3))
# # "Sardine"
#   generate.devs(N=50,rho = 0.9,sd.devs = 0.1,plot = TRUE)
#   mtext(side = 3,text = "Sardine-like",outer = FALSE)
# 
# # "Herring/anchovy"
#   generate.devs(N=50,rho = 0.5,sd.devs = 0.3,plot = TRUE)
#   mtext(side = 3,text = "Anchovy-like",outer = FALSE)
# 
# # "Menhaden"
#   generate.devs(N=50,rho = 0,sd.devs = 0.7,plot = TRUE)
#   mtext(side = 3,text = "Menhaden-like",outer = FALSE)
# dev.off()
