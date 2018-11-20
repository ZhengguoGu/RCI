######### some functions used in the simulation #############

GRM_func <- function(abil, itempar){
  
  numeritor <- exp(sweep((abil - itempar[, -1]), 1, itempar[, 1], "*"))
  P_star <- numeritor/(1+numeritor) # this is the "true response"
  
  response <- rowSums(P_star >= runif(nrow(itempar), min=0, max=1))
  
  return(response)
}
