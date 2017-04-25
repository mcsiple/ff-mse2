# Trend-based HCR
# This is a basic control rule that decides whether to fish or not based on the trend in observed biomass in the past three years. 
# Essentially, this is like a data-poor rule that people can follow when there is some survey data but no mechanistic model for estimating biomass
trend.rule <- function(B.2yr, B.prev, B.curr, F.const){
  #' @param B.2yr Observed biomass two years before the current year
  #' @param B.prev Observed biomass last year
  #' @param B.curr Observed biomass this year
  #' @param F.const Fishing rate when fishing is allowed
  #' @return The fishing rate F.const, given \code{B.2yr}, \code{B.prev}, and \code{B.curr} (trends in obs biomass in the past two years).
  r1 <- B.prev/B.2yr
  r2 <- B.curr/B.prev
  if(r1 > 1 && r2 > 1){F.const <- F.const}
  else F.const <- 0
  return(F.const)
}