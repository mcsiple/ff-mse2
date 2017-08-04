# This function estimates F_true, the true fishing rate based on the TAC. This is different than F from the HCR, because the TAC is determined with the assumption that the observed biomass is correct, i.e.,

# EXAMPLE - SEVEN AGES TOTAL
      # tac=20000
      # M = 0.5
      # sel <-  c(0.004, 0.143, 0.994, 0.840, 0.191, 0.024, 0.000)
      # obs.at.age.oneplus <- c(39869.064, 22180.144, 13353.527,  7767.168,  6295.520,  3154.487) # this is a vector of 6, has to be bc of errors
      # w.at.age <- c(0.0001281, 0.0002317, 0.0003285, 0.0003711, 0.0005371, 0.0004481)

calc.true.f <- function(tac.fn = tac, M.fn = M, sel.fn = sel, Btrue = NA, w.at.age =w.at.age){
        fn.to.minimize <- function(param,tac = tac.fn, M = M.fn, sel = sel.fn, Bobs_a.fn = Bobs_a,
                                   w.at.age.fn = w.at.age) {
          F.true <- param
          death_rate <- M + sel * F.true
          catch_age <- (Btrue * sel * F.true * (1-exp(-death_rate))) / death_rate
          #catch_age <- Bobs_a * (1-exp(M+sel*F.true)) / (M + F.true) # allowable catch at age
          pred.tac <- sum(catch_age)
          diff <- (sum(tac) - pred.tac)^2
          return(diff)
        }
  est <- optim(par=0.5, fn=fn.to.minimize, method="BFGS")
  return(est$par)
}

