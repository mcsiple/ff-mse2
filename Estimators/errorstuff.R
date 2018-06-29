
basedir <- "~/Dropbox/Chapter4-HarvestControlRules/Code/ff-mse2"
resultsdir <- "~/Dropbox/Chapter4-HarvestControlRules/Results"

source(file.path(basedir,"Recruitment/GenerateDevs.R")) 
source(file.path(basedir,"Estimators/CalcFTrue.R"))
source(file.path(basedir,"Run/HCR_Trajectory_NEW.R"))
source(file.path(basedir,"Estimators/Estimators.R"))
source(file.path(basedir,"Run/generate_M.R"))



test.ts <- c(129790, 210226, 268468, 390606, 398058, 397318, 360216, 358585, 
             374015, 362225, 371147, 404184, 389995, 351421, 303666, 310247, 
             267873, 239348, 233029, 222568, 199544, 332979, 333917, 288796, 
             251272, 273175, 275501, 256436, 274347, 245027, 249175, 268043, 
             275516, 390729, 404027, 353527, 315368, 279771, 252693, 328303, 
             419058, 401525, 338971, 272203, 236587, 424096, 428526, 915854, 
             893034, 792664, 600809, 535045, 406554, 336333, 314980, 263382, 
             343581, 314591, 275886, 269675, 285701, 252091, 251580, 245164, 
             299508, 275334, 270394, 305106, 297652, 232627, 221962, 211044, 
             255357, 229605, 202242, 195857, 340297, 494339, 457842, 422973, 
             457724, 500911, 440822, 594920, 572725, 480533, 360052, 338472, 
             257437, 272528, 372553, 370993, 305588, 270280, 277552, 239371, 
             212529, 217231, 313825, 390237, 361535, 339711, 313088, 267748, 
             209095, 289711, 260227, 211904, 161633, 155718, 120459, 105902, 
             100124, 96604, 122619, 120161, 145985, 250999, 362795, 404525, 
             375865, 329231, 279981, 289788, 261451, 224297, 205369, 280646, 
             515630, 527285, 497506, 535635, 537885, 516387, 509951, 474588, 
             371023, 305414, 299816, 257533, 292202, 266364, 236219, 546059, 
             578987, 489714, 383187, 358955, 340774, 336392, 389630, 359374, 
             337507, 346819, 307102, 256898, 328943, 310683, 296801, 238059, 
             230740, 338667, 449285, 447861, 395256, 378358, 301371, 228776, 
             166265, 154457, 201497, 240833, 299670, 357728, 430009, 405332, 
             438952, 381992, 343458, 320398, 335295, 291373, 253898, 251732, 
             230440, 224267, 220523, 207248, 181627, 182546, 169953, 155230, 
             190514, 269183, 292842, 278632, 316353, 323410, 318420, 314301, 
             406794, 533092, 557539, 771823, 766456, 696629, 560839, 454830, 
             339962, 279609, 273985, 274599, 265466, 326317, 411317, 377061, 
             401487, 343129, 297018, 246823, 225816, 197123, 165008, 248715, 
             287621, 342946, 331932, 334180, 402338, 372293, 331572, 267601, 
             253587, 214134, 226303, 238854, 266070, 351503, 331102, 275112, 
             225694, 209395, 316549, 299872, 252782, 189004, 166132, 120317, 
             115794, 109419)
# curly.phi <- rnorm(1,0,sig.s) # random deviations on top of the autocorrelation
# epsilon.curr <- rho * epsilon.prev + sqrt(1-(rho^2)) * curly.phi # error in the current year
sigma.test <- seq(0.25,0.35, length.out =50)
tim.params <- list(sigma0 = 0.2,tau0 = 0.1)
tau1 <-  (1/tim.params$tau0^2 + 1/tim.params$sigma0^2)^(-0.5)
sigma.save <- matrix(NA,2,length(sigma.test))
#tim.params$sigma0 <- sigma_est

for (i in 1:length(sigma.test)){

obs.ts <-eps.vec <- obs.ts.ac <-vector()
obs.ts[1] <- test.ts[1]*exp(rnorm(1,0,(1/tim.params$tau0^2 +1/tim.params$sigma0^2)^(-0.5)))


sig.s2 <- sigma.test[i] #(1/tim.params$tau0^2 +1/tim.params$sigma0^2)^(-0.5) # diff values
curly.phi2 <- rnorm(1,0,sig.s2) # random deviations on top of the autocorrelation

#epsilon.curr2 <- 1 # error in the current year
#err2 <- (epsilon.curr2-(0.5*sig.s2^2))
obs.ts.ac[1] <- test.ts[1] * exp(curly.phi2)

# sig.s <- (1/tim.params$tau0^2 +1/tim.params$sigma0^2)^(-0.5) # target sd 0.4
eps.vec[1] <- 1
#obs.ts[1] <- biomass.true[,yr] * exp(tim.rand.inits)
  
for(t in 2:250){
  obs.ts[t] <- tim.assessment(B = test.ts[t],Eprev = obs.ts[t-1],sigma0 = tim.params$sigma0, tau0 = tim.params$tau0,tau1 = tau1,err = ,err_a = )
  #Eprev = biomass[,yr-1],
  # B = biomass.true[,yr],
  # sigma0 = sigma0,
  # tau0 = tau0,
  # tau1 = tau1,
  # err_a = tim.rands[,yr]
  x<- c(0.1,0.2,50,40,30)
  x.prev <- c(0.2,0.8,60,50,40)
  
  eps.prev =0.5;sig.s=0.3;rho=0.5;curly.phi=-0.37
  by.age <- add.wied.error(biomass.true=x,epsilon.prev = eps.prev,sig.s = sig.s,rho = rho,curly.phi = curly.phi)
  sum(by.age$biomass.est)
  total <- add.wied.error(biomass.true=sum(x),epsilon.prev = eps.prev,sig.s = sig.s,rho = rho,curly.phi = curly.phi)
  sum(total$biomass.est)
  
  by.age2 <- tim.assessment(Eprev = x.prev,B = x,sigma0 = tim.params$sigma0, tau0 = tim.params$tau0,tau1 = tau1,err_a = rep(-0.08,times=length(x)))
  total.2 <- tim.assessment(B = sum(x),Eprev = sum(x.prev),sigma0 = tim.params$sigma0,tau0 = tim.params$tau0,tau1 = tau1,err = -0.08)
  
  sum(by.age2)
  total.2
    # tim.assessment(Eprev = biomass[,yr-1],
    #                B = biomass.true[,yr],
    #                sigma0 = sigma0,
    #                tau0 = tau0,
    #                tau1 = tau1,
    #                err_a = tim.rands[,yr])
  
  # if y==1, biomass.true[,yr] * exp(tim.rand.inits)
  est <- add.wied.error(biomass.true = test.ts[t],epsilon.prev = eps.vec[t-1],sig.s = sig.s2,rho=0.5)
  obs.ts.ac[t] <- est$biomass.est
  eps.vec[t] <- est$epsilon.curr
}
# par(mfrow=c(1,1))
# plot(test.ts,type='l')
# lines(obs.ts, col="red")
# lines(obs.ts.ac, col="blue")

 sigma.save[1,i] <- sd(log(test.ts/obs.ts))
 sigma.save[2,i] <-  sd(log(test.ts/obs.ts.ac))
 

}

plot(sigma.test,sigma.save[1,], ylim = c(0.05,0.5))
lines(sigma.test,sigma.save[2,])
# Figure out where these two sigma values intersect; this is where error(AC) ~~ error(DD)




hist(sigma.save[1,])

plot(test.ts, type = 'l')
lines(obs.ts.ac, col = 'blue')
lines(obs.ts, col = 'red')


epis.prev <- 1
rho <- 0

tau0 <- 0.5

sig.1 <- 0.4 # Specify the sd on 'add wied error'

sigma_est <- (1/(0.1^(-2)-1/tau0^2))^0.5 # Calculate the value for sigma2 given tau - max value is 0.5 (inf)
sig.2 <-sigma_est # specifiy the sd on 'tim.assessment'
sigma0 <- sig.2
ntest <- vector()
ntest2 <- vector()

for (i in 1:10000){ # Test it numerically
ntest[i] <-add.wied.error(epis.prev,sig.1,rho)
ntest2[i] <- tim.assessment(tau0, sig.2)
}

par(mfrow = c(2,1), mar = c(4,4,1,1))
nmean <- mean(ntest)
nsd <- sd(ntest)
x <- seq(-5,5, length.out = 100)
h <- hist(ntest, breaks = 20, density = 10,
          col = "gray",  main = paste('SD = ',format(nsd, digits = 2),', mean = ',format(nmean, digits = 2)))
xfit <- seq(min(ntest), max(ntest), length = 40) 
yfit <- dnorm(xfit, mean = mean(ntest), sd = sd(ntest)) 
yfit <- yfit * diff(h$mids[1:2]) * length(ntest) 

lines(xfit, yfit, col = "black", lwd = 2)

test.ts <- testie$biomass.oneplus.true
sd(ntest)
sd(ntest2)
sd(test.ts*ntest[1:250])
sd(test.ts*ntest2[1:250])

# Try the other error model and test the algebra
nmean <- mean(ntest2)
nsd <- sd(ntest2)
x <- seq(-5,5, length.out = 100)
h <- hist(ntest2, breaks = 20, density = 10,
          col = "gray", main = paste('SD = ',format(nsd, digits = 2),', mean = ',format(nmean, digits = 2))) 
xfit <- seq(min(ntest2), max(ntest2), length = 40) 
yfit <- dnorm(xfit, mean = mean(ntest2), sd = sd(ntest2)) 
yfit <- yfit * diff(h$mids[1:2]) * length(ntest2) 
lines(xfit, yfit, col = "black", lwd = 2)

sigma.crossval <- sqrt(mean(test.ts[2:250] / test.ts[1:249]) / 0.5)
