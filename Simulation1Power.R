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

set.seed(2)

######################################################
#########  Conditions
######################################################

test_length <- c(5, 15, 40)
item_character <- c("parallel", "non-parallel")
CO_effect <- c("non", "30%", "50%")  # carry-over effects
change_theta = c(0.5, 1)    

condition <- expand.grid(test_length, item_character, change_theta,  CO_effect)
colnames(condition) <- c("test_length", "item_character", "magnitude_change", "carry-over")

Final_result <- list()
num_test <- 1
while(num_test <= dim(condition)[1]){
  
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
  
  
  theta_pre <- rnorm(1000, mean = 0, sd = 1)
  theta_pre <- sort(theta_pre)
  theta_post <- theta_pre + condition[num_test, 3] 
  
  # some of the persons may show carry-over effects
  if(condition[num_test, 4] == "30%"){  #introducing carry-over effects, if any
    NoCarry_index <- sample(1000, floor(1000 * (1-0.3)), replace = FALSE)  #here we fix the persons who do not show carryover
  }else if (condition[num_test, 4] == "50%"){
    NoCarry_index <- sample(1000, floor(1000 * (1-0.5)), replace = FALSE)  #here we fix the persons who do not show carryover
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
    # 1. alternative equation 1
    SE1 <- sqrt(var_pre * (1 - r11) + var_post * (1 - r22))
    # 2. alternative equation 2
    SE2 <- sd_D * sqrt(1 - (r11 * var_pre + r22 * var_post - 2 * r12 * sd_pre * sd_post) / (var_pre + var_post - 2 * r12 * sd_pre * sd_post))
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
  
  Powers <- Reduce('+', sim_result) / 100  # parallel-generated 100 matrices, and we add these matrices together
  result <- cbind(theta_pre, Powers)
  colnames(result) <- c("theta", "eq0", "eq1", "eq2", "eq3")
  Final_result[[num_test]] <- result
  
  print(num_test)
  num_test = num_test + 1
}

save(Final_result, file = "simulation1_power.RData")
################## summarizing results ############

library(ggplot2)
library(gridExtra)
library(grid)
library(tidyr)
pic_function <- function(data){
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



### 1. test length against carryover effect, when identical items and change = .5 ####
condition[27, ]

p1 <- pic_function(data.frame(Final_result[[1]]))   #test length = 5, parallel, magnitude_change = .5, no carryover
p2 <- pic_function(data.frame(Final_result[[13]]))  #test length = 5, parallel, magnitude_change = .5, 30% carryover
p3 <- pic_function(data.frame(Final_result[[25]]))  #test length = 5, parallel, magnitude_change = .5, 50% carryover
p4 <- pic_function(data.frame(Final_result[[2]]))  #test length = 15, parallel, magnitude_change = .5, no carryover
p5 <- pic_function(data.frame(Final_result[[14]]))  #test length = 15, parallel, magnitude_change = .5, 30% carryover
p6 <- pic_function(data.frame(Final_result[[26]])) #test length = 15, parallel, magnitude_change = .5, 50% carryover
p7 <- pic_function(data.frame(Final_result[[3]]))  #test length = 40, parallel, magnitude_change = .5, no carryover
p8 <- pic_function(data.frame(Final_result[[15]]))  #test length = 40, parallel, magnitude_change = .5, 30% carryover
p9 <- pic_function(data.frame(Final_result[[27]])) #test length = 40, parallel, magnitude_change = .5, 50% carryover

grid.arrange(
  p1, p2, p3,
  p4, p5, p6,
  p7, p8, p9,
  nrow = 3
) #warnings are due to the truncation of xlim

### 2. test length against carryover effect, when non-identical items and change = .5 ####
condition[30, ]

p1 <- pic_function(data.frame(Final_result[[4]]))   #test length = 5, non-parallel, magnitude_change = .5, no carryover
p2 <- pic_function(data.frame(Final_result[[16]]))  #test length = 5, non-parallel, magnitude_change = .5, 30% carryover
p3 <- pic_function(data.frame(Final_result[[28]]))  #test length = 5, non-parallel, magnitude_change = .5, 50% carryover
p4 <- pic_function(data.frame(Final_result[[5]]))  #test length = 15, non-parallel, magnitude_change = .5, no carryover
p5 <- pic_function(data.frame(Final_result[[17]]))  #test length = 15, non-parallel, magnitude_change = .5, 30% carryover
p6 <- pic_function(data.frame(Final_result[[29]])) #test length = 15, non-parallel, magnitude_change = .5, 50% carryover
p7 <- pic_function(data.frame(Final_result[[6]]))  #test length = 40, non-parallel, magnitude_change = .5, no carryover
p8 <- pic_function(data.frame(Final_result[[18]]))  #test length = 40, non-parallel, magnitude_change = .5, 30% carryover
p9 <- pic_function(data.frame(Final_result[[30]])) #test length = 40, non-parallel, magnitude_change = .5, 50% carryover

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
