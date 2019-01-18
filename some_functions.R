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


# 3 carry-over effect for profile change

carry_over <- function(pre, post, rand_index){
  
  #pre: pretest scores (a vector)
  #post: posttest scores (a vector)
  #rand_index: the index marking the persons showing no carry-over effects
  
  new_post <- post
  
  ind1 <- (pre - post < -1) 
  new_post[ind1] <- pre[ind1] + 1
  
  ind2 <- (pre - post > 1)
  new_post[ind2] <- pre[ind2] - 1
  
  ind3 <- (abs(pre-post)==1 )
  new_post[ind3] <- pre[ind3]
  
  new_post[rand_index] <- post[rand_index]  #those who dont show effects are kept unchanged. 
  
  return(new_post)
}


# 4 Bonferroni
bonfer <- function(p_val, alpha_M){
  
}