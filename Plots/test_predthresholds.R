result.list <- results[[17]] 
sim=1
par(mfrow=c(2,2))
for(sim in 1:25){
plot(1:years.test,result.list$total.true.biomass[sim,],type='l',ylim=c(0,max(c(result.list$biomass[sim,],
                                                                    result.list$no.fishing.tb[sim,]))))
lines(1:years.test,result.list$no.fishing.tb[sim,],col="blue")
recruits.biomass <- result.list$total.true.biomass[sim,] - result.list$biomass[sim,] #result.list$biomass is the true one-plus biomass

mean.unfished <- mean(result.list$no.fishing.tb[sim,])
abline(h=0.2*mean.unfished,col="blue",lty=2)
lines(1:years.test, recruits.biomass,col="red")
}



# Code for looking at collapses
collapses <- collapse.index(x = result.list$biomass[sim,])