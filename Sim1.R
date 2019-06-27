######################################################
#########                                   ##########
#########  Simulation 1: unidimensional     ##########
#########                                   ##########
#########                                   ##########
######################################################

library(psychometric)
library(doSNOW)
library(doRNG)

source(file = "some_functions.R")

set.seed(2)

######################################################
#########  Conditions
######################################################

test_length <- c(5, 15, 40)
item_character <- c("parallel", "non-parallel")
CO_effect <- c("non", "30%", "50%")  # carry-over effects
perc_change = c(0.7, 0.5)    #70% showed change or 50% showed change

condition <- expand.grid(test_length, item_character, perc_change,  CO_effect)
colnames(condition) <- c("test_length", "item_character", "perc_change", "carry-over")

Final_result <- list()
Power_mean_median <- matrix(NA, nrow(condition), 4) 
Type1_mean_median <- matrix(NA, nrow(condition), 4)
num_test <- 1
while(num_test <= dim(condition)[1]){
  
  theta_pre <- rnorm(1000, mean = 0, sd = 1)
  #theta_pre <- sort(theta_pre)  #no need to sort it anymore
  theta_d <- replicate(1000, 0)
  index_change <- sample(1:1000, (condition[num_test, 3])*1000, replace = F)
  theta_d[index_change] <- 1  #these people show change
  theta_post <- theta_pre + theta_d 
  
  # some of the persons may show carry-over effects
  if(condition[num_test, 4] == "30%"){  #introducing carry-over effects, if any
    NoCarry_index <- sample(1000, floor(1000 * (1-0.3)), replace = FALSE)  #here we fix the persons who do not show carryover
  }else if (condition[num_test, 4] == "50%"){
    NoCarry_index <- sample(1000, floor(1000 * (1-0.5)), replace = FALSE)  #here we fix the persons who do not show carryover
  }
  
  
  ##############################################################
  
  n_rep <- 1 #note that we repeat 50 times, so that in each repitition, a new test (i.e., new set of item parameters) is generated.
  results <- matrix(0, 1000, 4)
  while(n_rep <= 50){
    
    if (condition[num_test, 2] == "parallel") {#identical item parameters
      
      itempar <- matrix(NA,condition[num_test, 1],5)
      itempar[,1] <- runif(1,1.5,2.5)   # discrimination
      avg_beta <- runif(1, 0, 1.25)
      itempar[,2] <- avg_beta - .75
      itempar[,3] <- avg_beta - .25
      itempar[,4] <- avg_beta + .25
      itempar[,5] <- avg_beta + .75
      
    } else { #non-identical item parameters
      
      itempar <- matrix(NA,condition[num_test, 1],5)
      itempar[,1] <- runif(condition[num_test, 1],1.5,2.5)  # discrimination
      avg_beta <- runif(condition[num_test, 1], 0, 1.25)
      itempar[,2] <- avg_beta - .75
      itempar[,3] <- avg_beta - .25
      itempar[,4] <- avg_beta + .25
      itempar[,5] <- avg_beta + .75
      
    }
    
    cl <- makeCluster(12)
    registerDoSNOW(cl)
    sim_result <- foreach(i = 1:100) %dorng% {
      
      pretest <- t(sapply(theta_pre, FUN = GRM_func,  itempar = itempar))
      posttest <- t(sapply(theta_post, FUN = GRM_func,  itempar = itempar))  
      
      if(condition[num_test, 4] != "non"){    #i.e., there is carryover effect
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
      # 1. using r11 instead of r12
      SE1 <- sqrt(2 * (1 - r11)) * sd_pre
      # 2. alternative equation 1
      SE2 <- sqrt(var_pre * (1 - r11) + var_post * (1 - r22))
      # 3. alternative equation 2 (equavalent to SE2)
      # SE2 <- sd_D * sqrt(1 - (r11 * var_pre + r22 * var_post - 2 * r12 * sd_pre * sd_post) / (var_pre + var_post - 2 * r12 * sd_pre * sd_post))
      # 3. alternative equation 3
      SE3 <- sd_D * sqrt(1 - rDD)
      
      sig_eq0 <- (abs(D_score / SE0) > 1.645)  
      sig_eq1 <- (abs(D_score / SE1) > 1.645)
      sig_eq2 <- (abs(D_score / SE2) > 1.645)
      sig_eq3 <- (abs(D_score / SE3) > 1.645)
      
      result <- cbind(sig_eq0, sig_eq1, sig_eq2, sig_eq3)
      return(result)
      
    }
    stopCluster(cl)
    
    results <- results + Reduce('+', sim_result) / 100  # parallel-generated 100 matrices, and we add these matrices together
    n_rep <- n_rep + 1
    
    print(n_rep)
  }
  
  results <- results/50 #average across 50 tests 

  result_power <- cbind(theta_pre[index_change],  results[index_change, ])  #these people show change --> power
  colnames(result_power) <- c("theta_pre", "eq0", "eq1", "eq2", "eq3")
  Power_mean_median[num_test, ] <- c(apply(result_power[, 2:5], 2, mean))  #note, I do not compute medians anymore, but the matrix is still Power_mean_median
  
  
  result_type1 <- cbind(theta_pre[setdiff(1:1000, index_change)], results[setdiff(1:1000, index_change), ]) #these people do not change --> type 1
  colnames(result_type1) <- c("theta_pre", "eq0", "eq1", "eq2", "eq3")
  Type1_mean_median[num_test, ] <- c(apply(result_type1[, 2:5], 2, mean))
  
  final_res <- list(result_power, result_type1) 
  
  Final_result[[num_test]] <- final_res
  
  print(num_test)
  num_test = num_test + 1
}

save(Final_result, file = "simulation1.RData")

################## summarizing results ###############################################

#1. Power: test length against carryover effect, when identical items, and 70% vs. 50% of people changes ####

condition[condition$item_character=="parallel" & condition$perc_change==0.7, ]
index_table_70 <- c(1, 13, 25, 2, 14, 26, 3, 15, 27)
condition[index_table_70, ]

condition[condition$item_character=="parallel" & condition$perc_change==0.5, ]
index_table_50 <- c(7, 19, 31, 8, 20, 32, 9, 21, 33)
condition[index_table_50, ]

PowerTable_parallel <- cbind(Power_mean_median[index_table_70, ], Power_mean_median[index_table_50, ])
write.csv(PowerTable_parallel, file = "PowerTable_parallel.csv")

#2. Power: test length against carryover effect, when non-identical items and 70% people change ####
condition[condition$item_character=="non-parallel" & condition$perc_change==0.7, ]
index_table_70 <- c(4, 16, 28, 5, 17, 29, 6, 18, 30)
condition[index_table_70, ]

condition[condition$item_character=="non-parallel" & condition$perc_change==0.5, ]
index_table_50 <- c(10, 22, 34, 11, 23, 35, 12, 24, 36)
condition[index_table_50, ]
PowerTable_nonparallel <- cbind(Power_mean_median[index_table_70, 1:4], Power_mean_median[index_table_50, 1:4])
write.csv(PowerTable_nonparallel, file = "PowerTable_nonparallel.csv")

#3. Type1: test length against carryover effect, when identical items, and 70% vs. 50% of people changes ####

condition[condition$item_character=="parallel" & condition$perc_change==0.7, ]
index_table_70 <- c(1, 13, 25, 2, 14, 26, 3, 15, 27)
condition[index_table_70, ]

condition[condition$item_character=="parallel" & condition$perc_change==0.5, ]
index_table_50 <- c(7, 19, 31, 8, 20, 32, 9, 21, 33)
condition[index_table_50, ]

Type1Table_parallel <- cbind(Type1_mean_median[index_table_70, 1:4], Type1_mean_median[index_table_50, 1:4])
write.csv(Type1Table_parallel, file = "Type1Table_parallel.csv")

#4. Type1: test length against carryover effect, when non-identical items and 70% people change ####
condition[condition$item_character=="non-parallel" & condition$perc_change==0.7, ]
index_table_70 <- c(4, 16, 28, 5, 17, 29, 6, 18, 30)
condition[index_table_70, ]

condition[condition$item_character=="non-parallel" & condition$perc_change==0.5, ]
index_table_50 <- c(10, 22, 34, 11, 23, 35, 12, 24, 36)
condition[index_table_50, ]
Type1Table_nonparallel <- cbind(Type1_mean_median[index_table_70, 1:4], Type1_mean_median[index_table_50, 1:4])
write.csv(Type1Table_nonparallel, file = "Type1Table_nonparallel.csv")


############################## END ############################################################


############################# below are not used anymore  ###################################
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyr)
pic_power <- function(data){
  data = data.frame(data)
  colnames(data)[2:5] <- c("Eq(5)", "Eq(11)", "Eq(14)", "Eq(15)")
  longdata<- gather(data, equation, result, "Eq(5)":"Eq(15)", factor_key=TRUE)
  p <- ggplot(longdata, aes(x = longdata$theta, y = longdata$result, colour = longdata$equation)) + 
    geom_line(aes(group = longdata$equation))  +
    scale_y_continuous(breaks = c(0, .5, 1), limits = c(0,1)) +
    xlim(-3, 3) + 
    theme(legend.position = "none", axis.text=element_text(size=12, color="black"), 
          axis.title.x = element_text(size=12), axis.title.y = element_text(size=12)) +
    labs(x=expression(theta[1]), y="Power") +
    facet_grid(equation ~ .)
  return(p)
} #!!! xlim and ylim are truncated, therefore when generating plots, we see warning messages.

pic_type1 <- function(data){
  data = data.frame(data)
  colnames(data)[2:5] <- c("Eq(5)", "Eq(11)", "Eq(14)", "Eq(15)")
  longdata<- gather(data, equation, result, "Eq(5)":"Eq(15)", factor_key=TRUE)
  p <- ggplot(longdata, aes(x = longdata$theta, y = longdata$result, colour = longdata$equation)) + 
    geom_line(aes(group = longdata$equation))  +
    scale_y_continuous(breaks = c(0, 0.1, .4), limits = c(0,.4)) +
    geom_hline(yintercept=0.1, linetype="dashed", color = "black") +
    xlim(-3, 3) + 
    theme(legend.position = "none", axis.text=element_text(size=12, color="black"), 
          axis.title.x = element_text(size=12), axis.title.y = element_text(size=12)) +
    labs(x=expression(theta[1]), y="Type-I Error Rate") +
    facet_grid(equation ~ .)
  return(p)
} #!!! xlim and ylim are truncated, therefore when generating plots, we see warning messages.


#############################################
############ Plots for Power ################
#############################################
### 1. test length against carryover effect, when identical items and 70% of people changes ####
condition[27, ]

p1 <- pic_power(data.frame(Final_result[[1]][[1]]))   #test length = 5, parallel, 70% people change, no carryover
p2 <- pic_power(data.frame(Final_result[[13]][[1]]))  #test length = 5, parallel, 70% people change, 30% carryover
p3 <- pic_power(data.frame(Final_result[[25]][[1]]))  #test length = 5, parallel, 70% people change, 50% carryover
p4 <- pic_power(data.frame(Final_result[[2]][[1]]))  #test length = 15, parallel, 70% people change, no carryover
p5 <- pic_power(data.frame(Final_result[[14]][[1]]))  #test length = 15, parallel, 70% people change, 30% carryover
p6 <- pic_power(data.frame(Final_result[[26]][[1]])) #test length = 15, parallel, 70% people change, 50% carryover
p7 <- pic_power(data.frame(Final_result[[3]][[1]]))  #test length = 40, parallel, 70% people change, no carryover
p8 <- pic_power(data.frame(Final_result[[15]][[1]]))  #test length = 40, parallel, 70% people change, 30% carryover
p9 <- pic_power(data.frame(Final_result[[27]][[1]])) #test length = 40, parallel, 70% people change, 50% carryover

grid.arrange(
  p1, p2, p3,
  p4, p5, p6,
  p7, p8, p9,
  nrow = 3
) #warnings are due to the truncation of xlim

### 2. test length against carryover effect, when non-identical items and 70% people change ####
condition[4, ]

p1 <- pic_power(data.frame(Final_result[[4]][[1]]))   #test length = 5, non-parallel, 70% people change, no carryover
p2 <- pic_power(data.frame(Final_result[[16]]))  #test length = 5, non-parallel, magnitude_change = .5, 30% carryover
p3 <- pic_power(data.frame(Final_result[[28]]))  #test length = 5, non-parallel, magnitude_change = .5, 50% carryover
p4 <- pic_power(data.frame(Final_result[[5]]))  #test length = 15, non-parallel, magnitude_change = .5, no carryover
p5 <- pic_power(data.frame(Final_result[[17]]))  #test length = 15, non-parallel, magnitude_change = .5, 30% carryover
p6 <- pic_power(data.frame(Final_result[[29]])) #test length = 15, non-parallel, magnitude_change = .5, 50% carryover
p7 <- pic_power(data.frame(Final_result[[6]]))  #test length = 40, non-parallel, magnitude_change = .5, no carryover
p8 <- pic_power(data.frame(Final_result[[18]]))  #test length = 40, non-parallel, magnitude_change = .5, 30% carryover
p9 <- pic_power(data.frame(Final_result[[30]])) #test length = 40, non-parallel, magnitude_change = .5, 50% carryover

grid.arrange(
  p1, p2, p3,
  p4, p5, p6,
  p7, p8, p9,
  nrow = 3
) #warnings are due to the truncation of xlim

### 3. test length against carryover effect, when identical items and change = 1 ####
condition[33, ]

p1 <- pic_function(data.frame(Final_result[[7]]))   #test length = 5, parallel, magnitude_change = 1, no carryover
p2 <- pic_function(data.frame(Final_result[[19]]))  #test length = 5, parallel, magnitude_change = 1, 30% carryover
p3 <- pic_function(data.frame(Final_result[[31]]))  #test length = 5, parallel, magnitude_change = 1, 50% carryover
p4 <- pic_function(data.frame(Final_result[[8]]))  #test length = 15, parallel, magnitude_change = 1, no carryover
p5 <- pic_function(data.frame(Final_result[[20]]))  #test length = 15, parallel, magnitude_change = 1, 30% carryover
p6 <- pic_function(data.frame(Final_result[[32]])) #test length = 15, parallel, magnitude_change = 1, 50% carryover
p7 <- pic_function(data.frame(Final_result[[9]]))  #test length = 40, parallel, magnitude_change = 1, no carryover
p8 <- pic_function(data.frame(Final_result[[21]]))  #test length = 40, parallel, magnitude_change = 1, 30% carryover
p9 <- pic_function(data.frame(Final_result[[33]])) #test length = 40, parallel, magnitude_change = 1, 50% carryover

grid.arrange(
  p1, p2, p3,
  p4, p5, p6,
  p7, p8, p9,
  nrow = 3
) #warnings are due to the truncation of xlim

### 4. test length against carryover effect, when non-identical items and change = 1 ####
condition[36, ]

p1 <- pic_function(data.frame(Final_result[[10]]))   #test length = 5, non-parallel, magnitude_change = 1, no carryover
p2 <- pic_function(data.frame(Final_result[[22]]))  #test length = 5, non-parallel, magnitude_change = 1, 30% carryover
p3 <- pic_function(data.frame(Final_result[[34]]))  #test length = 5, non-parallel, magnitude_change = 1, 50% carryover
p4 <- pic_function(data.frame(Final_result[[11]]))  #test length = 15, non-parallel, magnitude_change = 1, no carryover
p5 <- pic_function(data.frame(Final_result[[23]]))  #test length = 15, non-parallel, magnitude_change = 1, 30% carryover
p6 <- pic_function(data.frame(Final_result[[35]])) #test length = 15, non-parallel, magnitude_change = 1, 50% carryover
p7 <- pic_function(data.frame(Final_result[[12]]))  #test length = 40, non-parallel, magnitude_change = 1, no carryover
p8 <- pic_function(data.frame(Final_result[[24]]))  #test length = 40, non-parallel, magnitude_change = 1, 30% carryover
p9 <- pic_function(data.frame(Final_result[[36]])) #test length = 40, non-parallel, magnitude_change = 1, 50% carryover

grid.arrange(
  p1, p2, p3,
  p4, p5, p6,
  p7, p8, p9,
  nrow = 3
) #warnings are due to the truncation of xlim