# Libraries
library(doParallel)
library(doSNOW)
library(plyr)

# Set directories
basedir <- "/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Code/ff-mse2"
resultsdir <- "/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Results"
types <- c("Anchovy"="Anchovy",
           "Menhaden" = "Menhaden",
           "Sardine" = "Sardine")

# Parallelize so each simulation is on a different core:
registerDoParallel(4)
llply(.data=types,.fun = run.mod,.parallel = TRUE)
stopCluster()

# for (f in 1:2){
#   fftype = types[f]

      #run.mod <- function(fftype){
        basedir <- "/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Code/ff-mse2"
        resultsdir <- "/Users/mcsiple/Dropbox/Chapter4-HarvestControlRules/Results"
        #subDir <- "Anchovy" # Name of ff type
        subDir <- fftype
        
        #Set up other simulation params
        years.test = 250
        nsims = 1000
        tim.params = list(sigma0 = 0.2,tau0 = 0.1)
        sig.s = 0.3 #0.30 (changed to 0.001 to check whether catches were limited at higher F by observation error)
        R0.sens = NA #NO DYNAMIC R0 anymore-- ignore
        
        # Load packages
        library(reshape2)
        library(ggplot2)
        library(plyr)
        library(dplyr)
        
        # Load rev devs generating fxn, MSE main model, estimator fxns
        #toplot=FALSE      # Don't plot examples of rec trajectories
        source(file.path(basedir,"Recruitment/GenerateDevs.R")) 
        source(file.path(basedir,"Estimators/CalcFTrue.R"))
        source(file.path(basedir,"Run/HCR_Trajectory_NEW.R"))
        source(file.path(basedir,"Estimators/Estimators.R"))
        source(file.path(basedir,"Run/generate_M.R"))

        # Load control files & set parameter values
        #     Type     TargetSD
        # 1  Sardine 0.5117551
        # 2  Anchovy 0.3221512
        # 3 Menhaden 0.3009440
        # Sardines
        if(subDir == "Sardine"){
            source(file.path(basedir,"Ctl/Sardine_LHControl.R"))
            source(file.path(basedir,"Ctl/Sardine_FisheryControl.R"))
            # Sardine params
            recruit.sd <- 0.6
            recruit.rho <- 0.9
            }
        # Anchovy/Herring
        if(subDir == "Anchovy"){
            source(file.path(basedir,"Ctl/Anchovy_LHControl.R"))
            source(file.path(basedir,"Ctl/Anchovy_FisheryControl.R"))
            #Anchovy recruitment dev params
            recruit.sd <- 0.6
            recruit.rho <- 0.5
        }
        
        # Menhaden
        if(subDir == "Menhaden"){
          source(file.path(basedir,"Ctl/Menhaden_LHControl.R"))
          source(file.path(basedir,"Ctl/Menhaden_FisheryControl.R"))
          #Menhaden recruitment dev params
          recruit.sd <- 0.8
          recruit.rho <- 0.2
        }
        
        
        # Load harvest rules
        source(file.path(basedir,"Control Rules/smith_oceana.R"))
        source(file.path(basedir,"Control Rules/cfp.R"))
        source(file.path(basedir,"Control Rules/hockey-stick.R"))
        source(file.path(basedir,"Control Rules/trend-based-rule.R"))
        
        # If a results folder doesn't exist already, create one!
        resultsfolder <- paste(subDir,Sys.Date(),sep="")
        dir.create(file.path(resultsdir, resultsfolder))
        setwd(file.path(resultsdir, resultsfolder))
        
        # Scenarios
        h = c(0.9, 0.6)
        obs.error.type = c("AC","Tim")
        #obs.error.type = "noerror"
        HCR = c("constF","constF_HI")
        #HCR = c("cfp","constF","C1","C2","C3","constF_HI") # Took out trend because it was unrealistic-- but using trend in CPUE as adjustment (data-poor method) might be a good idea!
        M.type = c("constant") # took out "regimeshift" and "time-varying" to save time but can be added back in for sensitivity
        
        scenarios <- expand.grid(h,obs.error.type,HCR,recruit.sd,recruit.rho,M.type)
        colnames(scenarios) <- c("h","obs.error.type","HCR","recruit.sd","recruit.rho","M.type")
        nscenarios <- nrow(scenarios)
        scenarios$scenario <- 1:nscenarios #Label them so it's easier to find/index em later
        
        
        write.table(scenarios,file = "Scenario_Table.txt")
        
        CFP <- C1 <- C2 <- C3 <- constF <- trend <- constF_HI <- 
                    list(biomass.oneplus.true=matrix(nrow = nsims,ncol = years.test), 
                    total.catch=matrix(nrow = nsims,ncol = years.test),
                    fishing= matrix(nrow = nsims,ncol = years.test),
                    intended.f=matrix(nrow = nsims,ncol = years.test),
                    rec= matrix(nrow = nsims,ncol = years.test),
                    biomass.oneplus.obs = matrix(nrow = nsims,ncol = years.test),
                    biomass.total.true = matrix(nrow = nsims,ncol = years.test),
                    no.fishing.tb = matrix(nrow = nsims,ncol = years.test))
        
        
        # Test params and runs to make sure they look good ------------------------
                steepness = scenarios$h[1]
                obs.type <- scenarios$obs.error.type[1]
                HCR <- scenarios$HCR[1]
                recruit.sd = .6 #scenarios$recruit.sd[1]
                recruit.rho = .9 #scenarios$recruit.rho[1]
                equilib = getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=steepness)
                rec.dev.test <- generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd)
                test.constF <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test,rec.ram = NA, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "constF", const.f.rate = 0, steepness = steepness,obs.type = obs.type,equilib=equilib,R0.traj = R0.sens, tim.params = tim.params,time.var.m = NA,sig.s = sig.s)
                nofish <- melt(test.constF[-c(1,4,9,10)])
                nofish$year <- rep(1:years.test,times=length(unique(nofish$L1)))
                ggplot(nofish,aes(x=year,y=value)) + geom_line() + facet_wrap(~L1,scales = "free_y") #+ xlim(c(150,250))
                
        # Test pop and see if it crashes ------------------------------------------
        # Check that population will still sometimes collapse even without fishing.
        nexamples <- 1 # How many time series do you want to plot
        var.to.plot <- "rec"
        set.seed(123)
              big.df <- data.frame()
              for(i in 1:nexamples){
              steepness = scenarios$h[2]
              obs.type <- scenarios$obs.error.type[2]
              HCR <- scenarios$HCR[1]
              recruit.sd = .6 #scenarios$recruit.sd[1]
              recruit.rho = .9 #scenarios$recruit.rho[1]
              equilib = getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=steepness)
              rec.dev.test <- generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd)
              test.constF <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "constF", const.f.rate = 0, steepness = steepness,obs.type = obs.type,equilib=equilib,R0.traj = R0.sens, tim.params = tim.params,time.var.m = NA,sig.s = sig.s)
              nofish <- melt(test.constF[-c(1,4,9,10)])
              nofish$year <- rep(1:years.test,times=length(unique(nofish$L1)))
              vars.to.plot <- subset(nofish, L1 %in% var.to.plot)
              vars.to.plot$rep <- paste(i)
              big.df <- rbind(big.df,vars.to.plot)
              }
        
              
        ndf <- big.df %>% group_by(rep) %>% mutate(low.B.thresh = 0.2*mean(value)) %>% as.data.frame()
        ggplot(ndf,aes(x=year,y=value)) + geom_line() + facet_wrap(~rep,scales = "free_y") + geom_line(aes(y=low.B.thresh),col="red") #xlim(c(150,250)
        
        #New, just for the recruitment time series
        # fff <- nofish %>% subset(L1=="biomass.oneplus.true") %>% mutate(low.B.thresh = 0.2*mean(value,na.rm=TRUE))
        # nofish %>% subset(L1=="biomass.oneplus.true") %>% as.data.frame() %>% ggplot(aes(x=year,y=value))+geom_line()+ geom_hline(yintercept = 15639.94,col = "red")
        # 
        ###########################################################################
        # SIMULATIONS -------------------------------------------------------------
        ###########################################################################
        
        #
        for(s in 1:nscenarios){  #
          steepness = scenarios$h[s]
          obs.type <- scenarios$obs.error.type[s]
          HCR <- scenarios$HCR[s]
          recruit.sd = scenarios$recruit.sd[s]
          recruit.rho = scenarios$recruit.rho[s]
          M.type = scenarios$M.type[s]
        
          equilib = getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=steepness) # NO recruitment devs used in the equilibrium calculations, so don't need to embed in the loop
          const.f.rate = 0.5*equilib$Fmsy 
                  no.fishing <- matrix(NA, nrow = nsims, ncol = years.test)
                  set.seed(123) # Start each round of sims at same random seed
                  time.var.m <- NA # Base case: M constant 
                  if(M.type == "timevar"){time.var.m <- rw.M(Mbar = lh.test$M, rho.m = 0.6, sigma.m = 0.2,n = years.test)}
                  if(M.type == "regimeshift"){time.var.m <- regime.M(Mbar = lh.test$M,cutoff.yr = 201,n = years.test)}
                  for (sim in 1:nsims){
                        rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
                        F0.Type <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "constF", const.f.rate = 0, steepness = steepness,obs.type = obs.type,equilib=equilib,R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s)$biomass.total.true # Only need to do this 1x for each simulation (not repeat for each CR) because the seed is the same and there is no fishing.
                        no.fishing[sim,] <- F0.Type
                          }
                
                
          if(HCR=="cfp"){    
          set.seed(123) # Start each round of sims at same random seed
          for (sim in 1:nsims){
                rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
                expt.cfp <- calc.trajectory(lh = lh.test,obs.cv = NA, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "cfp",equilib = equilib,steepness=steepness,obs.type = obs.type,R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s,rec.ram = NA)
              
            CFP[["biomass.oneplus.true"]][sim,] <- expt.cfp$biomass.oneplus.true # This is the true one-plus biomass
            CFP[["total.catch"]][sim,] <- expt.cfp$total.catch  
            CFP[["fishing"]][sim,] <- expt.cfp$fishing
            CFP[["intended.f"]][sim,] <- expt.cfp$intended.f    
            CFP[["rec"]][sim,] <- expt.cfp$rec
            CFP[["biomass.oneplus.obs"]][sim,] <- expt.cfp$biomass.oneplus.obs        # This is the observed one-plus biomass
            CFP[["biomass.total.true"]][sim,] <- expt.cfp$biomass.total.true          # This is the true total biomass
            CFP[["no.fishing.tb"]] <- no.fishing    # True total biomass with no fishing
            
          }
            save(CFP,file=paste("All",s,"CFP",".RData",sep="_"))
            }
            
          if(HCR=="constF"){
            set.seed(123) # same seed
            for (sim in 1:nsims){
                  rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
                  expt.constF <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "constF", const.f.rate = const.f.rate, steepness = steepness,obs.type = obs.type,equilib=equilib,R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s)
                  
            constF[["biomass.oneplus.true"]][sim,] <- expt.constF$biomass.oneplus.true
            constF[["total.catch"]][sim,] <- expt.constF$total.catch
            constF[["fishing"]][sim,] <- expt.constF$fishing
            constF[["intended.f"]][sim,] <- expt.constF$intended.f   
            constF[["rec"]][sim,] <- expt.constF$rec
            constF[["biomass.oneplus.obs"]][sim,] <- expt.constF$biomass.oneplus.obs
            constF[["biomass.total.true"]][sim,] <- expt.constF$biomass.total.true
            constF[["no.fishing.tb"]] <- no.fishing
            }
            save(constF,file=paste("All",s,"constF",".RData",sep="_"))
            }
            
          if(HCR=="C1"){
            set.seed(123) # same seed
            for (sim in 1:nsims){
                rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
                expt.c1 <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "C1",equilib = equilib,steepness=steepness,obs.type = obs.type,R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s)
               
            C1[["biomass.oneplus.true"]][sim,] <- expt.c1$biomass.oneplus.true
            C1[["total.catch"]][sim,] <- expt.c1$total.catch
            C1[["fishing"]][sim,] <- expt.c1$fishing
            C1[["intended.f"]][sim,] <- expt.c1$intended.f    
            C1[["rec"]][sim,] <- expt.c1$rec
            C1[["biomass.oneplus.obs"]][sim,] <- expt.c1$biomass.oneplus.obs
            C1[["biomass.total.true"]][sim,] <- expt.c1$biomass.total.true
            C1[["no.fishing.tb"]] <- no.fishing
            }
            save(C1,file=paste("All",s,"C1",".RData",sep="_")) 
          }
          
          if(HCR=="C2"){
            set.seed(123) # same seed
            for (sim in 1:nsims){
              rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
              expt.c2 <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "C2",equilib = equilib,steepness=steepness,obs.type = obs.type,R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s)
              
              C2[["biomass.oneplus.true"]][sim,] <- expt.c2$biomass.oneplus.true
              C2[["total.catch"]][sim,] <- expt.c2$total.catch
              C2[["fishing"]][sim,] <- expt.c2$fishing
              C2[["intended.f"]][sim,] <- expt.c2$intended.f    
              C2[["rec"]][sim,] <- expt.c2$rec
              #C2[["depl"]][sim,] <- expt.c2$depl
              C2[["biomass.oneplus.obs"]][sim,] <- expt.c2$biomass.oneplus.obs
              C2[["biomass.total.true"]][sim,] <- expt.c2$biomass.total.true
              C2[["no.fishing.tb"]] <- no.fishing
            }
            save(C2,file=paste("All",s,"C2",".RData",sep="_")) 
          }
                  
          if(HCR=="C3"){
            set.seed(123) # same seed
            for (sim in 1:nsims){
              rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
              expt.c3 <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "C3",equilib = equilib,steepness=steepness,obs.type = obs.type,R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s)
              
              C3[["biomass.oneplus.true"]][sim,] <- expt.c3$biomass.oneplus.true
              C3[["total.catch"]][sim,] <- expt.c3$total.catch
              C3[["fishing"]][sim,] <- expt.c3$fishing
              C3[["intended.f"]][sim,] <- expt.c3$intended.f    
              C3[["rec"]][sim,] <- expt.c3$rec
              C3[["biomass.oneplus.obs"]][sim,] <- expt.c3$biomass.oneplus.obs
              C3[["biomass.total.true"]][sim,] <- expt.c3$biomass.total.true
              C3[["no.fishing.tb"]] <- no.fishing
            }
            save(C3,file=paste("All",s,"C3",".RData",sep="_")) 
          }
          
          
          if(HCR=="trend"){
            set.seed(123) # same seed
            for (sim in 1:nsims){
              rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
              expt.trend <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "trend",const.f.rate = 0.6,equilib = equilib,steepness=steepness,obs.type = obs.type,R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s)
              
              trend[["biomass.oneplus.true"]][sim,] <- expt.trend$biomass.oneplus.true
              trend[["total.catch"]][sim,] <- expt.trend$total.catch
              trend[["fishing"]][sim,] <- expt.trend$fishing
              trend[["intended.f"]][sim,] <- expt.trend$intended.f    # True total biomass with no fishing
              trend[["rec"]][sim,] <- expt.trend$rec
              trend[["biomass.oneplus.obs"]][sim,] <- expt.trend$biomass.oneplus.obs # Observed one-plus biomass
              trend[["biomass.total.true"]][sim,] <- expt.trend$biomass.total.true
              trend[["no.fishing.tb"]] <- no.fishing
            }
            save(trend,file=paste("All",s,"trend",".RData",sep="_")) 
          }
                  
                  if(HCR=="constF_HI"){
                    const.f.rate = equilib$Fmsy
                    set.seed(123) # same seed
                    for (sim in 1:nsims){
                      rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
                      expt.constF_HI <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,hcr.type = "constF", const.f.rate = const.f.rate, steepness = steepness,obs.type = obs.type,equilib=equilib,R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s)
                      
                      constF_HI[["biomass.oneplus.true"]][sim,] <- expt.constF_HI$biomass.oneplus.true
                      constF_HI[["total.catch"]][sim,] <- expt.constF_HI$total.catch
                      constF_HI[["fishing"]][sim,] <- expt.constF_HI$fishing
                      constF_HI[["intended.f"]][sim,] <- expt.constF_HI$intended.f   
                      constF_HI[["rec"]][sim,] <- expt.constF_HI$rec
                      constF_HI[["biomass.oneplus.obs"]][sim,] <- expt.constF_HI$biomass.oneplus.obs
                      constF_HI[["biomass.total.true"]][sim,] <- expt.constF_HI$biomass.total.true
                      constF_HI[["no.fishing.tb"]] <- no.fishing
                    }
                    save(constF_HI,file=paste("All",s,"constF_HI",".RData",sep="_"))
                  }
        } # end of scenarios loop
        
        #print(proc.time() - tm)
      }
  