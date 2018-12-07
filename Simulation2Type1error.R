######################################################
#########                                   ##########
#########  Simulation 2: profile            ##########
#########  Type 1 error rate                ##########
#########                                   ##########
######################################################

library(psychometric)
library(doSNOW)
library(doRNG)
library(MASS)

source(file = "some_functions.R")

set.seed(1)

######################################################
#########  Conditions
######################################################
test_length <- c(5, 15, 40)  #all subtests consist of 5/15/40 items
item_character <- c("parallel", "non-parallel")
CO_effect <- c("non", "30%", "50%")  # carry-over effects

condition <- expand.grid(test_length, item_character, CO_effect)
colnames(condition) <- c("test_length", "item_character", "carry-over")

Final_result <- list()
num_test <- 1
while(num_test <= dim(condition)[1]){
  
  if (condition[num_test, 2] == "parallel") {
    
    #generate items for subtest 1
    itempar_sub1 <- matrix(NA,condition[num_test, 1],5)
    itempar_sub1[,1] <- runif(1,1.5,2.5)   # discrimination
    avg_beta <- runif(1, 0, 1.25)
    itempar_sub1[,2] <- avg_beta - .75
    itempar_sub1[,3] <- avg_beta - .25
    itempar_sub1[,4] <- avg_beta + .25
    itempar_sub1[,5] <- avg_beta + .75
    
    #generate items for subtest 2
    itempar_sub2 <- matrix(NA,condition[num_test, 1],5)
    itempar_sub2[,1] <- runif(1,1.5,2.5)   # discrimination
    avg_beta <- runif(1, 0, 1.25)
    itempar_sub2[,2] <- avg_beta - .75
    itempar_sub2[,3] <- avg_beta - .25
    itempar_sub2[,4] <- avg_beta + .25
    itempar_sub2[,5] <- avg_beta + .75
    
    #generate items for subtest 3
    itempar_sub3 <- matrix(NA,condition[num_test, 1],5)
    itempar_sub3[,1] <- runif(1,1.5,2.5)   # discrimination
    avg_beta <- runif(1, 0, 1.25)
    itempar_sub3[,2] <- avg_beta - .75
    itempar_sub3[,3] <- avg_beta - .25
    itempar_sub3[,4] <- avg_beta + .25
    itempar_sub3[,5] <- avg_beta + .75
    
    itempar <- list(itempar_sub1, itempar_sub2, itempar_sub3)
    
  } else {#non-parallel
    
    #generate items for subtest 1
    itempar_sub1 <- matrix(NA,condition[num_test, 1],5)
    itempar_sub1[,1] <- runif(condition[num_test, 1],1.5,2.5)  # discrimination
    avg_beta <- runif(condition[num_test, 1], 0, 1.25)
    itempar_sub1[,2] <- avg_beta - .75
    itempar_sub1[,3] <- avg_beta - .25
    itempar_sub1[,4] <- avg_beta + .25
    itempar_sub1[,5] <- avg_beta + .75
    
    #generate items for subtest 2
    itempar_sub2 <- matrix(NA,condition[num_test, 1],5)
    itempar_sub2[,1] <- runif(condition[num_test, 1],1.5,2.5)  # discrimination
    avg_beta <- runif(condition[num_test, 1], 0, 1.25)
    itempar_sub2[,2] <- avg_beta - .75
    itempar_sub2[,3] <- avg_beta - .25
    itempar_sub2[,4] <- avg_beta + .25
    itempar_sub2[,5] <- avg_beta + .75
    
    #generate items for subtest 3
    itempar_sub3 <- matrix(NA,condition[num_test, 1],5)
    itempar_sub3[,1] <- runif(condition[num_test, 1],1.5,2.5)  # discrimination
    avg_beta <- runif(condition[num_test, 1], 0, 1.25)
    itempar_sub3[,2] <- avg_beta - .75
    itempar_sub3[,3] <- avg_beta - .25
    itempar_sub3[,4] <- avg_beta + .25
    itempar_sub3[,5] <- avg_beta + .75
    
    itempar <- list(itempar_sub1, itempar_sub2, itempar_sub3)
    
  }
  
  
  theta <- MASS::mvrnorm(1000, mu = c(0, 0, 0), Sigma = matrix(c(1, .1, .1, .1, 1, .1, .1, .1, 1), 3, 3), empirical = FALSE)
  
  # some of the persons may show carry-over effects
  if(condition[num_test, 3] == "30%"){  #introducing carry-over effects, if any
    NoCarry_index <- sample(1000, floor(1000 * (1-0.3)), replace = FALSE)  #here we fix the persons who do NOT show carryover
  }else if (condition[num_test, 3] == "50%"){
    NoCarry_index <- sample(1000, floor(1000 * (1-0.5)), replace = FALSE)  #here we fix the persons who do NOT show carryover
  }
  
 
  cl <- makeCluster(2)
  registerDoSNOW(cl)
  sim_result <- foreach(i = 1:100) %dorng% {
    
    
    RCI_0 <- 0
    RCI_1 <- 0
    RCI_2 <- 0
    RCI_3 <- 0
    for(j in 1:3){ #jth subtest
      
      pretest <- t(sapply(theta[, j], FUN = GRM_func,  itempar = itempar[[j]]))
      posttest <- t(sapply(theta[, j], FUN = GRM_func,  itempar = itempar[[j]]))  #there is no change, and if there would be no carry-over effects
      
      if(condition[num_test, 3] == "30%"){  #introducing carry-over effects, if any
        posttest <- carry_over(pretest, posttest, NoCarry_index)
      }else if (condition[num_test, 3] == "50%"){
        posttest <- carry_over(pretest, posttest, NoCarry_index)
      }
      
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
      
      RCI_0 <- RCI_0 + (D_score/SE0)^2
      RCI_1 <- RCI_1 + (D_score/SE1)^2
      RCI_2 <- RCI_2 + (D_score/SE2)^2
      RCI_3 <- RCI_3 + (D_score/SE3)^2
    }
    
    
    
    sig_eq0 <- (RCI_0 > 5.991 | RCI_0 < 0.103)  #chi-square test at .1 level, two tails
    sig_eq1 <- (RCI_1 > 5.991 | RCI_1 < 0.103)
    sig_eq2 <- (RCI_2 > 5.991 | RCI_2 < 0.103)
    sig_eq3 <- (RCI_3 > 5.991 | RCI_3 < 0.103)
    
    result <- cbind(sig_eq0, sig_eq1, sig_eq2, sig_eq3)
    return(result)
    
  }
  stopCluster(cl)
  
  Type1error <- Reduce('+', sim_result) / 100  # parallel-generated 100 matrices, and we add these matrices together, and then compute the empirical Type 1 error rate 
  result <- cbind(theta, Type1error)
  colnames(result) <- c("theta", "eq0", "eq1", "eq2", "eq3")
  Final_result[[num_test]] <- result
  
  num_test = num_test + 1
}