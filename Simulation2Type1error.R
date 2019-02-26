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
library(stats)


source(file = "some_functions.R")

set.seed(3)

#some values that should not be changed

alpha_M <- .1/3 #experimental alpha = .1, M=3 (i.e., 3 subtests)
Q <- .15 #false discovery rate

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
  theta <- theta[order(stats::mahalanobis(theta, 0, cov=matrix(c(1, .1, .1, .1, 1, .1, .1, .1, 1), 3, 3), inverted = FALSE)),] #accending order in terms of mahalanobis distance

  
  # some of the persons may show carry-over effects
  if(condition[num_test, 3] == "30%"){  #introducing carry-over effects, if any
    NoCarry_index <- sample(1000, floor(1000 * (1-0.3)), replace = FALSE)  #here we fix the persons who do NOT show carryover
  }else if (condition[num_test, 3] == "50%"){
    NoCarry_index <- sample(1000, floor(1000 * (1-0.5)), replace = FALSE)  #here we fix the persons who do NOT show carryover
  }
  
 
  cl <- makeCluster(4)
  registerDoSNOW(cl)
  sim_result <- foreach(i = 1:100) %dorng% {
    
    RCI_0 <- 0
    RCI_1 <- 0
    RCI_2 <- 0
    RCI_3 <- 0
    
    p_vecregister <- list()
      
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
      
      p_vecregister[[j]] <- cbind(t(sapply(D_score/SE0, FUN = pvalue, two_tail = 2, alpha_M = alpha_M)),  #notice that t(sapply(..., FUN=pvalue),... ) generate two columns, one for p values, one for sig/nonsig
                                  t(sapply(D_score/SE1, FUN = pvalue, two_tail = 2, alpha_M = alpha_M)),
                                  t(sapply(D_score/SE2, FUN = pvalue, two_tail = 2, alpha_M = alpha_M)),
                                  t(sapply(D_score/SE3, FUN = pvalue, two_tail = 2, alpha_M = alpha_M))
                                  )  #p_vecregister[[j]] is a 1000x8 matrix, rows represent persons
      
       
    }
    
    sig_eq0 <- (RCI_0 > 4.6052)  #chi-square test at .1 level
    sig_eq1 <- (RCI_1 > 4.6052)
    sig_eq2 <- (RCI_2 > 4.6052)
    sig_eq3 <- (RCI_3 > 4.6052)
    
  
    posthoc_bonf <- cbind(p_vecregister[[1]][, 2], p_vecregister[[2]][, 2], p_vecregister[[3]][, 2],  #eq0: subtest 1, 2, 3
                          p_vecregister[[1]][, 4], p_vecregister[[2]][, 4], p_vecregister[[3]][, 4],  #eq1: subtest 1, 2, 3
                          p_vecregister[[1]][, 6], p_vecregister[[2]][, 6], p_vecregister[[3]][, 6],  #eq2: subtest 1, 2, 3
                          p_vecregister[[1]][, 8], p_vecregister[[2]][, 8], p_vecregister[[3]][, 8])  #eq3: subtest 1, 2, 3

    colnames(posthoc_bonf) <- c("bonf_eq0_sub1", "bonf_eq0_sub2", "bonf_eq0_sub3", 
                                 "bonf_eq1_sub1", "bonf_eq1_sub2", "bonf_eq1_sub3",
                                 "bonf_eq2_sub1", "bonf_eq2_sub2", "bonf_eq2_sub3",
                                 "bonf_eq3_sub1", "bonf_eq3_sub2", "bonf_eq3_sub3")
    
    posthoc_pvalues <- cbind(p_vecregister[[1]][, 1], p_vecregister[[2]][, 1], p_vecregister[[3]][, 1],  #eq0: subtest 1, 2, 3
                             p_vecregister[[1]][, 3], p_vecregister[[2]][, 3], p_vecregister[[3]][, 3],  #eq1: subtest 1, 2, 3
                             p_vecregister[[1]][, 5], p_vecregister[[2]][, 5], p_vecregister[[3]][, 5],  #eq2: subtest 1, 2, 3
                             p_vecregister[[1]][, 7], p_vecregister[[2]][, 7], p_vecregister[[3]][, 7])  #eq3: subtest 1, 2, 3
    
    BenHochresults <- cbind(apply(posthoc_pvalues[, 1:3], MARGIN = 2, FUN = Ben_Hoch, Q = Q),  #eq0: subtest 1, 2, 3
                            apply(posthoc_pvalues[, 4:6], MARGIN = 2, FUN = Ben_Hoch, Q = Q),  #eq1: subtest 1, 2, 3
                            apply(posthoc_pvalues[, 7:9], MARGIN = 2, FUN = Ben_Hoch, Q = Q),  #eq2: subtest 1, 2, 3
                            apply(posthoc_pvalues[, 10:12], MARGIN = 2, FUN = Ben_Hoch, Q = Q))  #eq3: subtest 1, 2, 3
                            
    colnames(BenHochresults) <- c("BenH_eq0_sub1", "BenH_eq0_sub2", "BenH_eq0_sub3", 
                                  "BenH_eq1_sub1", "BenH_eq1_sub2", "BenH_eq1_sub3",
                                  "BenH_eq2_sub1", "BenH_eq2_sub2", "BenH_eq2_sub3",
                                  "BenH_eq3_sub1", "BenH_eq3_sub2", "BenH_eq3_sub3")
    
    result <- cbind(sig_eq0, sig_eq1, sig_eq2, sig_eq3, posthoc_bonf, BenHochresults)
    return(result)
    
  }
  stopCluster(cl)
  
  Type1error <- Reduce('+', sim_result) / 100  # parallel-generated 100 matrices, and we add these matrices together, and then compute the empirical Type 1 error rate 
  result <- cbind(theta, Type1error)
  colnames(result) <- c("theta1", "theta2", "theta3", "omni_eq0", "omni_eq1", "omni_eq2", "omni_eq3",
                        "bonf_eq0_sub1", "bonf_eq0_sub2", "bonf_eq0_sub3", 
                        "bonf_eq1_sub1", "bonf_eq1_sub2", "bonf_eq1_sub3",
                        "bonf_eq2_sub1", "bonf_eq2_sub2", "bonf_eq2_sub3",
                        "bonf_eq3_sub1", "bonf_eq3_sub2", "bonf_eq3_sub3",
                        "BenH_eq0_sub1", "BenH_eq0_sub2", "BenH_eq0_sub3", 
                        "BenH_eq1_sub1", "BenH_eq1_sub2", "BenH_eq1_sub3",
                        "BenH_eq2_sub1", "BenH_eq2_sub2", "BenH_eq2_sub3",
                        "BenH_eq3_sub1", "BenH_eq3_sub2", "BenH_eq3_sub3")
  Final_result[[num_test]] <- result
  
  print(num_test)
  num_test = num_test + 1
}

