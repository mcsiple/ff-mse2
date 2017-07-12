# Demo figure for control rules
source(file.path(basedir,"Control Rules/smith_oceana.R"))
source(file.path(basedir,"Control Rules/cfp.R"))
source(file.path(basedir,"Control Rules/hockey-stick.R"))
source(file.path(basedir,"Control Rules/trend-based-rule.R"))

B <- 1:10000
B0 <- 7000
Bmsy <- 3000
Fmsy = 0.7
m = 0.8

C1.vec <- C2.vec <- C3.vec <- cfp.vec <- vector()
for(i in 1:length(B)){
C1.vec[i] <- calc.F.oceana(Bt = B[i],Blim = 0.4*B0,Btarget = 0.8*B0, M = m)
C2.vec[i] <- calc.F.oceana(Bt = B[i], Blim = 0.1*B0, Btarget = 0.8*B0, M = m)
C3.vec[i] <- calc.F.stick(Bt = B[i], Blim = 0.4*B0, Btarget = 0.8*B0, Fmax = Fmsy)
cfp.vec[i] <- calc.F.cfp(prevCatch = 10,Bt = B[i], Btarget = Bmsy, Fmax = Fmsy,Blim = 0.5*Bmsy )
}

pdf("ControlRules.pdf",width = 7,height = 7,useDingbats = FALSE)
lwdp = 3
par(las=2) # Rotate axis labels
plot(C1.vec,type='l',col=hcr.colors[1],lwd=lwdp, ylim=c(0,0.9),
     axes=FALSE,xlab="Biomass",ylab="Fishing rate (F)") #xaxt="n",yaxt="n"

ticks = c(0,0.1*B0,0.5*Bmsy,0.4*B0,max(B))
axis(side = 1, at = ticks,labels = c(" ","0.1*B0","0.5Bmsy","0.4*B0"," "))
axis(side = 2)

lines(C2.vec, col=hcr.colors[2],lwd=lwdp)
lines(C3.vec, col=hcr.colors[3],lwd=lwdp)
lines(cfp.vec,col=hcr.colors[5],lwd=lwdp)
#abline(h = 0.5*Fmsy,col=hcr.colors[4],lwd = lwdp)
legend("bottomright",legend = c("C1","C2","C3","Stability-favoring"),bty = "n",lwd=rep(lwdp,times=4),col=hcr.colors[c(1,2,3,5)])
dev.off()
