######################################################
#########                                   ##########
#########  Simulation 1: unidimensional     ##########
#########  Type 1 error rate                ##########
#########                                   ##########
######################################################

library(psychometric)
library(doSNOW)
library(doRNG)

source(file = "some_functions.R")

set.seed(1)

######################################################
#########  Conditions
######################################################

test_length <- c(5, 15, 40)
item_character <- c("parallel", "non-parallel")
CO_effect <- c("non", "30%", "50%")  # carry-over effects
# change_theta = 0  # this is for calculating type 1 error rate

condition <- expand.grid(test_length, item_character, CO_effect)
colnames(condition) <- c("test_length", "item_character", "carry-over")

Final_result <- list()
num_test <- 1
while(num_test <= dim(condition)[1]){
  
  if (condition[num_test, 2] == "parallel") { #identical item parameters
    
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
  
  Semi_result_eq0 <- list()
  Semi_result_eq1 <- list()
  Semi_result_eq2 <- list()
  Semi_result_eq3 <- list()
  #Semi_result <- list()
  rep_n <- 1
  while(rep_n <= 50){
    theta <- rnorm(1000, mean = 0, sd = 1)
    theta <- sort(theta)
  
    theta_index <- 1:1000
    theta_RCI <- array(NA)  #night person indices in theta sample, whose RCI will be computed. 
    quantile_values <- quantile(theta, c(.1, .2, .3, .4, .5, .6, .7, .8, .9)) 
    theta_RCI[1] <- theta_index[abs(theta - quantile_values[1]) ==min(abs(theta - quantile_values[1]))][1]  #note, in case more than one theta satisfies, we choose the first one (i.e., [1]). 
    theta_RCI[2] <- theta_index[abs(theta - quantile_values[2]) ==min(abs(theta - quantile_values[2]))][1]
    theta_RCI[3] <- theta_index[abs(theta - quantile_values[3]) ==min(abs(theta - quantile_values[3]))][1]
    theta_RCI[4] <- theta_index[abs(theta - quantile_values[4]) ==min(abs(theta - quantile_values[4]))][1]
    theta_RCI[5] <- theta_index[abs(theta - quantile_values[5]) ==min(abs(theta - quantile_values[5]))][1]
    theta_RCI[6] <- theta_index[abs(theta - quantile_values[6]) ==min(abs(theta - quantile_values[6]))][1]
    theta_RCI[7] <- theta_index[abs(theta - quantile_values[7]) ==min(abs(theta - quantile_values[7]))][1]
    theta_RCI[8] <- theta_index[abs(theta - quantile_values[8]) ==min(abs(theta - quantile_values[8]))][1]
    theta_RCI[9] <- theta_index[abs(theta - quantile_values[9]) ==min(abs(theta - quantile_values[9]))][1]
  # some of the persons may show carry-over effects
  
    if(condition[num_test, 3] == "30%"){  #introducing carry-over effects, if any
      NoCarry_index <- sample(1000, floor(1000 * (1-0.3)), replace = FALSE)  #here we fix the persons who do not show carryover
    }else if (condition[num_test, 3] == "50%"){
      NoCarry_index <- sample(1000, floor(1000 * (1-0.5)), replace = FALSE)  #here we fix the persons who do not show carryover
    }
  
  
    cl <- makeCluster(12)  #for parallel computing: set the number of cores
    registerDoSNOW(cl)
    sim_result <- foreach(i = 1:100) %dorng% {
    
      pretest <- t(sapply(theta, FUN = GRM_func,  itempar = itempar))
      posttest <- t(sapply(theta, FUN = GRM_func,  itempar = itempar))  #there is no change, and if there would be no carry-over effects
    
      if(condition[num_test, 3] != "non"){  #i.e., there is carryover effect
        posttest_temp <- carry_over(pretest, posttest, NoCarry_index)
        posttest_temp[theta_RCI, ] <- posttest[theta_RCI, ]  # this is to make sure that carryover effect does not happen to the 9 persons (10th, 20th, ..., 90th percetiles)
        posttest <- posttest_temp
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
    
      sig_eq0 <- (abs(D_score / SE0) > 1.645)  #thus, false positive
      sig_eq1 <- (abs(D_score / SE1) > 1.645)
      sig_eq2 <- (abs(D_score / SE2) > 1.645)
      sig_eq3 <- (abs(D_score / SE3) > 1.645)
    
      result <- cbind(sig_eq0, sig_eq1, sig_eq2, sig_eq3)
      return(result)
    
    }
    stopCluster(cl)
  
    Type1error <- Reduce('+', sim_result) / 100  # parallel-generated 100 matrices, and we add these matrices together, and then compute the empirical Type 1 error rate 
    result <- cbind(Type1error)
    colnames(result) <- c("eq0", "eq1", "eq2", "eq3")
    
    #Semi_result[[rep_n]] <- result[theta_RCI, ] # we need the 9 persons. 
    Semi_result_eq0[[rep_n]] <- result[theta_RCI, 1] # we need the 9 persons. 
    Semi_result_eq1[[rep_n]] <- result[theta_RCI, 2] # we need the 9 persons. 
    Semi_result_eq2[[rep_n]] <- result[theta_RCI, 3] # we need the 9 persons. 
    Semi_result_eq3[[rep_n]] <- result[theta_RCI, 4] # we need the 9 persons. 
    print(paste("rep: ", rep_n))
    rep_n = rep_n + 1
    
  }
  
  eq0_result <- apply(matrix(unlist(Semi_result_eq0), nrow = 9 , byrow = F), MARGIN = 1, mean)
  eq1_result <- apply(matrix(unlist(Semi_result_eq1), nrow = 9 , byrow = F), MARGIN = 1, mean)
  eq2_result <- apply(matrix(unlist(Semi_result_eq2), nrow = 9 , byrow = F), MARGIN = 1, mean)
  eq3_result <- apply(matrix(unlist(Semi_result_eq3), nrow = 9 , byrow = F), MARGIN = 1, mean)
  
  Final_result[[num_test]] <- cbind(eq0_result, eq1_result, eq2_result, eq3_result)
  colnames(Final_result[[num_test]])  <- c("eq0", "eq1", "eq2", "eq3")
  print(num_test)
  num_test = num_test + 1
}

save(Final_result, file = "simulation1_type1error.RData")

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
    scale_y_continuous(breaks = c(0, 0.1, .4), limits = c(0,.4)) +
    geom_hline(yintercept=0.1, linetype="dashed", color = "black") +
    xlim(-3, 3) + 
    theme(legend.position = "none", axis.text=element_text(size=12, color="black"), 
          axis.title.x = element_text(size=12), axis.title.y = element_text(size=12)) +
    labs(x=expression(theta[1]), y="Type-I Error Rate") +
    facet_grid(equation ~ .)
  return(p)
} #!!! xlim and ylim are truncated, therefore when generating plots, we see warning messages.



### 1. test length against carryover effect, when identical items  ####

condition[15, ]
p1 <- pic_function(data.frame(Final_result[[1]]))  #test length = 5, no carryover
p2 <- pic_function(data.frame(Final_result[[7]]))  #test length = 5, 30% carryover
p3 <- pic_function(data.frame(Final_result[[13]])) #test length = 5, 50% carryover
p4 <- pic_function(data.frame(Final_result[[2]]))  #test length = 15, no carryover
p5 <- pic_function(data.frame(Final_result[[8]]))  #test length = 15, 30% carryover
p6 <- pic_function(data.frame(Final_result[[14]])) #test length = 15, 50% carryover
p7 <- pic_function(data.frame(Final_result[[3]]))  #test length = 40, no carryover
p8 <- pic_function(data.frame(Final_result[[9]]))  #test length = 40, 30% carryover
p9 <- pic_function(data.frame(Final_result[[15]])) #test length = 40, 50% carryover

grid.arrange(
  p1, p2, p3,
  p4, p5, p6,
  p7, p8, p9,
  nrow = 3
) #warnings are due to the truncation of xlim
  
###  2. test length against carryover effect, when non-identical items  ####
condition[18, ]
p1 <- pic_function(data.frame(Final_result[[4]]))  #test length = 5, no carryover
p2 <- pic_function(data.frame(Final_result[[10]]))  #test length = 5, 30% carryover
p3 <- pic_function(data.frame(Final_result[[16]])) #test length = 5, 50% carryover
p4 <- pic_function(data.frame(Final_result[[5]]))  #test length = 15, no carryover
p5 <- pic_function(data.frame(Final_result[[11]]))  #test length = 15, 30% carryover
p6 <- pic_function(data.frame(Final_result[[17]])) #test length = 15, 50% carryover
p7 <- pic_function(data.frame(Final_result[[6]]))  #test length = 40, no carryover
p8 <- pic_function(data.frame(Final_result[[12]]))  #test length = 40, 30% carryover
p9 <- pic_function(data.frame(Final_result[[18]])) #test length = 40, 50% carryover

grid.arrange(
  p1, p2, p3,
  p4, p5, p6,
  p7, p8, p9,
  nrow = 3
)



####### backup functions, not used #########
pic_function <- function(data){
  p <- ggplot(data, aes(data$theta)) + 
    geom_line(aes(y = data$eq0, colour = "Eq0")) + 
    geom_line(aes(y = data$eq1, colour = "Eq1")) +
    geom_line(aes(y = data$eq2, colour = "Eq2")) +
    geom_line(aes(y = data$eq3, colour = "Eq3")) +
    geom_hline(yintercept=0.1, linetype="dashed", color = "gray") +
    scale_y_continuous(breaks = c(0, 0.1, .25, .5, .75, 1), limits = c(0,1)) +
    xlim(-3, 3) +  #!!! by truncating the range of x, the plot will generate warnings 
    labs(x=expression(theta[pre]), y="Type-I error rate") +
    theme(legend.title = element_blank(), legend.position = c(.85, .8), legend.text=element_text(size=10), axis.text=element_text(size=10), 
          axis.title.x = element_text(size=10), axis.title.y = element_text(size=10))
  return(p)
}
