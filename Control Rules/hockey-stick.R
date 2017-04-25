# "Hockey stick" function for a harvest control rule. 
# This function calculates f given a biomass level and the parameters for the h.s. function
example=FALSE

calc.F.stick <- function(Bt, Blim, Btarget, Fmax){
  # Blim is the biomass at which fishing is no longer allowed
  # B target is the biomass where you could catch @ rate Fmax - above Btarget you catch max rate, below Btarget you catch linearly less until you reach Blim.
  f <- NA
  slope = Btarget/(Btarget-Blim)      # slope of the diagonal part of the hockey stick
  adj.constant <- Btarget/Fmax     # scales y axis to max fishing mortality
  if (Bt <= Blim) {f = 0}
  if (Bt > Blim & Bt <= Btarget) {f <- slope * (Bt-Blim) / adj.constant}
  if (Bt > Btarget) {f = Fmax} #Blim / adj.constant
  return(f)
}


# calc.F.stick(Bt = 30,Blim = 20, Btarget = 50, Fmax = 0.4)
if(example==TRUE){
              B0 = 110.548
              BMSY = 31.49111
              FMSY = 1
              par(mfrow=c(3,2))
                  
              # Here is an example
              rule.table <- data.frame(Blim=c(0.5*BMSY,0.4*B0,0.1*B0,0.1*B0,0.5*BMSY),Btarget=c(BMSY,0.8*B0,BMSY,BMSY,BMSY),Fmax=c(FMSY,0.5*lh.test$M,FMSY,FMSY,FMSY))
              rule.names <- c("Classic US hockey stick","Oceana recommendations","40-10 rule","40-10 rule w/ buffer","Common Fisheries Policy","Constant F")
              for(n in 1:nrow(rule.table)){
              x<- seq(1,100,by=0.1)         # Biomass range
              y <- vector(length=length(x)) # Fishing mortality
              Blim = rule.table$Blim[n]
              Btarget = rule.table$Btarget[n]
              Fmax = rule.table$Fmax[n]
              
                  for(i in 1:length(x)){
                    y[i] <- calc.F.stick(Bt=x[i],Blim=Blim,Btarget=Btarget,Fmax = Fmax)
                  }
              
              plot(x,y,type='l',lwd=3, xlab='Biomass',ylab="F",axes=F,ylim=c(0,1.5),main=rule.names[n])
              text(x = mean(x)*1.1,y = 0.33*Fmax, labels = expression('B'[target]))
              text(labels = expression('F'[max]),x=0,y=Fmax, adj = c(-0.1,0))
              text(labels = expression('B'[lim]),x=Blim,y=0,pos=4)
              axis(side = 1,labels = FALSE)
              axis(side = 2,labels = FALSE)
              abline(v=Btarget,lty=2)
              abline(h=Fmax,lty=3)
              }
              
              
              plot(x,rep(0.6,times=length(x)),type='l',lwd=3, xlab='Biomass',ylab="F",axes=F,ylim=c(0,1.5),main="Constant F")
              # text(x = mean(x)*1.1,y = 0.33*Fmax, labels = expression('B'[target]))
              # text(labels = expression('F'[max]),x=0,y=Fmax, adj = c(-0.1,0))
              # text(labels = expression('B'[lim]),x=Blim,y=0,pos=4)
              axis(side = 1,labels = FALSE)
              axis(side = 2,labels = FALSE)
              text(labels = expression('F'[constant]),x=0,y=0.6, adj = c(-0.2,-0.2))
              
              
              # Here is a figure showing the function, for use in paper or proposal
              # setwd('~/Dropbox/Chapter 4 - Harvest Control Rules/Figures')
                      # tiff(file='HCR_Hockey.tiff',res=600,width=6,height=6,units='in',pointsize=12)
                      # par(mfrow=c(1,1))
                      # plot(x,y,type='l',lwd=3, xlab='Biomass',ylab="F",axes=F)
                      # text(x = mean(x)*1.1,y = 0.33*Fmax, labels = expression('B'[target]))
                      # text(labels = expression('F'[max]),x=0,y=Fmax, adj = c(-0.1,0))
                      # text(labels = expression('B'[lim]),x=Blim,y=0,pos=4)
                      # axis(side = 1,labels = FALSE)
                      # axis(side = 2,labels = FALSE)
                      # abline(v=Btarget,lty=2)
                      # abline(h=Fmax,lty=3)
                      # dev.off()
}