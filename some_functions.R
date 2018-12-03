######### some functions used in the simulation #############

# 1. GRM
GRM_func <- function(abil, itempar){
  
  numeritor <- exp(sweep((abil - itempar[, -1]), 1, itempar[, 1], "*"))
  P_star <- numeritor/(1+numeritor) # this is the "true response"
  
  response <- rowSums(P_star >= runif(nrow(itempar), min=0, max=1))
  
  return(response)
}


# 2. p-value
pvalue <- function(z_score, two_tail){
  
  if(missing(two_tail)){
    two_tail <- 2  #i.e., 2 tails
  }
  
  if(two_tail == 2){
    if(z_score < 0){
      p_value <- 2 * stats::pnorm(z_score)
    }else if (z_score > 0){
      p_value <- 2 * (1 - stats::pnorm(z_score))
    }else{
      p_value <- 0
    }
  }else if (two_tail == 1){
    if(z_score < 0){
      p_value <- stats::pnorm(z_score)
    }else if (z_score > 0){
      p_value <- 1 - stats::pnorm(z_score)
    }else{
      p_value <- 0
    }
  }else{
    print("two_tail allows for 2 values only: 1: one-tail test, and 2: two-tail test")
  }
    
  return(p_value)
}


# 3. carry-over effect

#--------------------------------------------------------------------------
# Mimicking carry-over effects
#
# Last update: 2018/10/24, percentage of persons showing carry-over effects
#--------------------------------------------------------------------------

carry_over <- function(pre, post, proc_N){
  
  #pre: pretest scores (a vector)
  #post: posttest scores (a vector)
  #proc_N: percentage of persons showing carry-over effects (a scalar)
  
  strong_post <- post
  weak_post <- post  
  
  ind1 <- (pre - post < -1) 
  strong_post[ind1] <- pre[ind1] + 1
  weak_post[ind1] <- post[ind1] - 1
  
  ind2 <- (pre - post > 1)
  strong_post[ind2] <- pre[ind2] - 1
  weak_post[ind2] <- post[ind2] + 1
  
  ind3 <- (abs(pre-post)==1 )
  strong_post[ind3] <- pre[ind3]
  weak_post[ind3] <- post[ind3]
  
  
  #taking into account the persentage of persons showing carry-over effects
  N <- length(pre)
  
  rand_index <- sample(N, floor(N * (1-proc_N)), replace = FALSE)  #records the index of persons showing NO carry-over effects! see below
  
  strong_post[rand_index] <- post[rand_index]  #those who dont show effects are kept unchanged. 
  weak_post[rand_index] <- post[rand_index]
  
  return(list(strong_post, weak_post))
}