# Conservative harvest control rule #1
# Loosely based on 
# This is the "nice to fish" endpoint
# This info comes from the Lenfest Forage Fish guidelines and an email from Ben Enticknap at Oceana, who provided some guidance. 
calc.F.oceana <- function(Bt, Blim, Btarget, M){
  Fmax = 0.5*M #The max F is determined from natural mortality
  #Btarget = 0.75*B0
  if(Bt==0){f=0}else{
  f <- NA
  slope = Btarget/(Btarget-Blim)      # slope of the diagonal part of the hockey stick
  adj.constant <- Btarget/Fmax     # scales y axis to max fishing mortality
  if (Bt <= Blim) {f = 0}
  if (Bt > Blim & Bt <= Btarget) {f <- slope * (Bt-Blim) / adj.constant}
  if (Bt > Btarget) {f = Fmax }
  }
  return(f)
}

# 
# calc.F.oceana(Bt = 14,Blim = 0.25*B0, B0 = B0,FMSY = FMSY)
# 
# x <- 1:1000
# for (i in 1:length(x)){
#   y[i] <- calc.F.oceana(Bt = x[i],Blim = 0.25*B0, B0 = B0,FMSY = FMSY)
# }
# 
# 
# plot(x,y,type='l')