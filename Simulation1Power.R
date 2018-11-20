######################################################
#########                                   ##########
#########  Simulation 1: unidimensional     ##########
#########  Power                            ##########
#########                                   ##########
######################################################

library(psychometric)
library(doSNOW)
library(doRNG)

source(file = "some_functions.R")
######################################################
#########  Conditions
######################################################

test_length <- c(5, 15, 40)
item_character <- c("parallel", "non-parallel")
change_theta = c(0.5, 1)    

condition <- expand.grid(test_length, item_character, change_theta)
colnames(condition) <- c("test_length", "item_character", "magnitude_change")


Final_result <- list()
num_test <- 1
while(num_test <= dim(condition)[1]){
  
  if (condition[num_test, 2] == "parallel") {
    
    itempar <- matrix(NA,condition[num_test, 1],5)
    itempar[,1] <- runif(1,1.5,2.5)   # discrimination
    avg_beta <- runif(1, 0, 1.25)
    itempar[,2] <- avg_beta - .75
    itempar[,3] <- avg_beta - .25
    itempar[,4] <- avg_beta + .25
    itempar[,5] <- avg_beta + .75
    
  } else {
    
    itempar <- matrix(NA,condition[num_test, 1],5)
    itempar[,1] <- runif(condition[num_test, 1],1.5,2.5)  # discrimination
    avg_beta <- runif(condition[num_test, 1], 0, 1.25)
    itempar[,2] <- avg_beta - .75
    itempar[,3] <- avg_beta - .25
    itempar[,4] <- avg_beta + .25
    itempar[,5] <- avg_beta + .75
    
  }
  
  
  theta_pre <- seq(-3, 3, length.out = 1000)
  theta_post <- theta_pre + condition[num_test, 3] 
    
  cl <- makeCluster(2)
  registerDoSNOW(cl)
  sim_result <- foreach(i = 1:100) %dorng% {
    
    pretest <- t(sapply(theta_pre, FUN = GRM_func,  itempar = itempar))
    posttest <- t(sapply(theta_post, FUN = GRM_func,  itempar = itempar))  
    
    sum_pre <- rowSums(pretest)
    sum_post <- rowSums(posttest)
    
    
    r12 <- cor(sum_pre, sum_post)
    var_pre <- var(sum_pre)
    sd_pre <- sd(sum_pre)
    var_post <- var(sum_post)
    sd_post <- sd(sum_post)
    
    D_score <- sum_post - sum_pre
    sd_D <- sd(D_score)
    rDD <- psychometric::alpha(posttest - pretest)
    
    r11 <- psychometric::alpha(pretest)
    r22 <- psychometric::alpha(posttest)
    # 0. orginal equation for sigma_E_D_v 
    SE0 <- sqrt(2 * (1 - r12)) * sd_pre 
    # 1. alternative equation 1
    SE1 <- sqrt(var_pre * (1 - r11) + var_post * (1 - r22))
    # 2. alternative equation 2
    SE2 <- sd_D * sqrt(1 - (r11 * var_pre + r22 * var_post - 2 * r12 * sd_pre * sd_post) / (var_pre + var_post - 2 * r12 * sd_pre * sd_post))
    # 3. alternative equation 3
    SE3 <- sd_D * sqrt(1 - rDD)
    
    sig_eq0 <- (abs(D_score / SE0) > 1.645)  #thus, false positive
    sig_eq1 <- (abs(D_score / SE1) > 1.645)
    sig_eq2 <- (abs(D_score / SE2) > 1.645)
    sig_eq3 <- (abs(D_score / SE3) > 1.645)
    
    result <- cbind(sig_eq0, sig_eq1, sig_eq2, sig_eq3)
    return(result)
    
  }
  stopCluster(cl)
  
  Powers <- Reduce('+', sim_result) / 100  # parallel-generated 100 matrices, and we add these matrices together, and then compute the empirical Type 1 error rate 
  result <- cbind(theta_pre, Powers)
  colnames(result) <- c("theta", "eq0", "eq1", "eq2", "eq3")
  Final_result[[num_test]] <- result
  
  num_test = num_test + 1
}

