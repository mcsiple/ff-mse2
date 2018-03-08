# This function estimates F_true, the true fishing rate based on the TAC. This is because the true F is the connection between the estimation and operating models.
# This "true F" is different than F from the HCR (the "target F"), because the TAC is determined by "managers" in the model with the assumption that observed biomass is correct.

# EXAMPLE W/ SEVEN AGES TOTAL
      # tac = 20000
      # M = 0.5
      # sel <-  c(0.004, 0.143, 0.994, 0.840, 0.191, 0.024, 0.000)
      # obs.at.age <- c(40000, 39869.064, 22180.144, 13353.527,  7767.168,  6295.520,  3154.487) # this is a vector of 6, has to be bc of errors
      # w.at.age <- c(0.0001281, 0.0002317, 0.0003285, 0.0003711, 0.0005371, 0.0004481)
      # calc.true.f(Btrue.fn = obs.at.age)

calc.true.f <- function(tac.fn = tac, M.fn = M, sel.fn = sel, Btrue.fn = NA, w.at.age = w.at.age){
        #' @description given the TAC, and the true biomass, what is the true F that will result in catches = TAC?
        #' @param tac.fn - vector of allowable catches at age (calculate from B_obs and F.intended)
        #' @param M.fn - natural mortality; drawn from OM 
        #' @param sel.fn - vector of selectivity at age
        #' @param Btrue.fn - vector of true biomass at age in yr y
        #' @param w.at.age - weight at age; drawn from OM
        fn.to.minimize <- function(param,tac = tac.fn, M = M.fn, sel = sel.fn,
                                   w.at.age.fn = w.at.age, Btrue = Btrue.fn){
          F.true <- param       # solving for F.true
          death_rate <- M + sel * F.true
          catch_age <- (Btrue * sel * F.true * (1-exp(-death_rate))) / death_rate # the biomass at age caught by the fishery
          pred.tac <- sum(catch_age)
          diff <- (sum(tac) - pred.tac)^2
          return(diff)
        }
  est <- optim(par=0.5, fn=fn.to.minimize, method="BFGS")
  return(est$par)
}

