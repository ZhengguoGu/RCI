######################################################
#########                                   ##########
#########  Simulation 2: profile            ##########
#########                                   ##########
#########                                   ##########
######################################################

library(psychometric)
library(doSNOW)
library(doRNG)
library(MASS)
library(stats)


source(file = "some_functions.R")

set.seed(4)

#some values that should not be changed
alpha_M <- .1/3 #experimental alpha = .1, M=3 (i.e., 3 subtests)
Q <- .15 #false discovery rate

######################################################
#########  Conditions
######################################################
test_length <- c(5, 15, 40)  #all subtests consist of 5/15/40 items
item_character <- c("parallel", "non-parallel")
#change_theta <- list(c(.5, .5, .5), c(1, 1, 1), c(0.1,0.5, 0.9))
perc_change <- c(0.7, 0.5) # percentage of persons change
CO_effect <- c("non", "30%", "50%")  # carry-over effects
condition <- expand.grid(test_length, item_character, perc_change, CO_effect)
colnames(condition) <- c("test_length", "item_character", "perc_change","carry-over")

cov_mat <- matrix(c(1, .5, .5, .5, 1, .5, .5, .5, 1), 3, 3) # covariance matrix for theta_pretest

OMNI_power <- matrix(NA, nrow(condition), 4) 
OMNI_type1 <- matrix(NA, nrow(condition), 4)
POSThoc_power <- matrix(NA, nrow(condition), 24) 
POSThoc_type1 <- matrix(NA, nrow(condition), 24)
num_test <- 1
while(num_test <= dim(condition)[1]){
  
  theta_pre <- MASS::mvrnorm(1000, mu = c(0, 0, 0), Sigma = cov_mat, empirical = FALSE)
  index_change <- sample(1:1000, (condition[num_test, 3])*1000, replace = F) #these people change
  theta_post <- theta_pre
  theta_post[index_change, ] <- sweep(theta_pre[index_change, ], MARGIN = 2, c(0.1,0.5, 0.9), "+")  # theta_d = c(0.1,0.5, 0.9)
  
  # some of the persons may show carry-over effects
  if(condition[num_test, 4] == "30%"){  #introducing carry-over effects, if any
    NoCarry_index <- sample(1000, floor(1000 * (1-0.3)), replace = FALSE)  #here we fix the persons who do NOT show carryover
  }else if (condition[num_test, 4] == "50%"){
    NoCarry_index <- sample(1000, floor(1000 * (1-0.5)), replace = FALSE)  #here we fix the persons who do NOT show carryover
  }
  
  
  ##############################################################
  
  n_rep <- 1 #note that we repeat 50 times, so that in each repitition, a new test (i.e., new set of item parameters) is generated.
  
  Results <- 0 # note, this will be broadcasted and becomes a matrix. 
  while(n_rep <= 50){
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
    
    cl <- makeCluster(12)
    registerDoSNOW(cl)
    sim_result <- foreach(i = 1:100) %dorng% {
      
      RCI_0 <- 0
      RCI_1 <- 0
      RCI_2 <- 0
      RCI_3 <- 0
      
      p_vecregister <- list()
      
      for(j in 1:3){ #jth subtest
        
        pretest <- t(sapply(theta_pre[, j], FUN = GRM_func,  itempar = itempar[[j]]))
        posttest <- t(sapply(theta_post[, j], FUN = GRM_func,  itempar = itempar[[j]]))  # if there would be no carry-over effects
        
        if(condition[num_test, 4] == "30%"){  #introducing carry-over effects, if any
          posttest <- carry_over(pretest, posttest, NoCarry_index)
        }else if (condition[num_test, 4] == "50%"){
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
      
      sig_eq0 <- (RCI_0 > 6.251)  #chi-square test at .1 level with df = 3
      sig_eq1 <- (RCI_1 > 6.251)
      sig_eq2 <- (RCI_2 > 6.251)
      sig_eq3 <- (RCI_3 > 6.251)
      
      
      posthoc_bonf <- cbind(p_vecregister[[1]][, 2], p_vecregister[[2]][, 2], p_vecregister[[3]][, 2],  #eq0: subtest 1, 2, 3
                            p_vecregister[[1]][, 4], p_vecregister[[2]][, 4], p_vecregister[[3]][, 4],  #eq1: subtest 1, 2, 3
                            p_vecregister[[1]][, 6], p_vecregister[[2]][, 6], p_vecregister[[3]][, 6],  #eq2: subtest 1, 2, 3
                            p_vecregister[[1]][, 8], p_vecregister[[2]][, 8], p_vecregister[[3]][, 8])  #eq3: subtest 1, 2, 3
      
      #colnames(posthoc_bonf) <- c("bonf_eq0_sub1", "bonf_eq0_sub2", "bonf_eq0_sub3", 
      #"bonf_eq1_sub1", "bonf_eq1_sub2", "bonf_eq1_sub3",
      #"bonf_eq2_sub1", "bonf_eq2_sub2", "bonf_eq2_sub3",
      #"bonf_eq3_sub1", "bonf_eq3_sub2", "bonf_eq3_sub3")
      
      posthoc_pvalues <- cbind(p_vecregister[[1]][, 1], p_vecregister[[2]][, 1], p_vecregister[[3]][, 1],  #eq0: subtest 1, 2, 3
                               p_vecregister[[1]][, 3], p_vecregister[[2]][, 3], p_vecregister[[3]][, 3],  #eq1: subtest 1, 2, 3
                               p_vecregister[[1]][, 5], p_vecregister[[2]][, 5], p_vecregister[[3]][, 5],  #eq2: subtest 1, 2, 3
                               p_vecregister[[1]][, 7], p_vecregister[[2]][, 7], p_vecregister[[3]][, 7])  #eq3: subtest 1, 2, 3
      
      BenHochresults <- cbind(apply(posthoc_pvalues[, 1:3], MARGIN = 2, FUN = Ben_Hoch, Q = Q),  #eq0: subtest 1, 2, 3
                              apply(posthoc_pvalues[, 4:6], MARGIN = 2, FUN = Ben_Hoch, Q = Q),  #eq1: subtest 1, 2, 3
                              apply(posthoc_pvalues[, 7:9], MARGIN = 2, FUN = Ben_Hoch, Q = Q),  #eq2: subtest 1, 2, 3
                              apply(posthoc_pvalues[, 10:12], MARGIN = 2, FUN = Ben_Hoch, Q = Q))  #eq3: subtest 1, 2, 3
      
      #olnames(BenHochresults) <- c("BenH_eq0_sub1", "BenH_eq0_sub2", "BenH_eq0_sub3", 
      # "BenH_eq1_sub1", "BenH_eq1_sub2", "BenH_eq1_sub3",
      #"BenH_eq2_sub1", "BenH_eq2_sub2", "BenH_eq2_sub3",
      # "BenH_eq3_sub1", "BenH_eq3_sub2", "BenH_eq3_sub3")
      
      result <- cbind(sig_eq0, sig_eq1, sig_eq2, sig_eq3, posthoc_bonf, BenHochresults)
      return(result)
      
    }
    stopCluster(cl)
    
    Results <- Results + Reduce('+', sim_result) / 100  # parallel-generated 100 matrices, and we add these matrices together, and then compute the empirical Type 1 error rate and power
  
    n_rep <- n_rep + 1
    print(n_rep)  
  }
  
  Results <- Results/50 #average across 50 tests 
  
  OMNI_power[num_test, ] <- colMeans(Results[index_change, 1:4])
  OMNI_type1[num_test, ] <- colMeans(Results[setdiff(1:1000, index_change), 1:4])
  
  POSThoc_power[num_test, ] <- colMeans(Results[index_change, 5:28])
  colnames(POSThoc_power) <- c("bonf_eq0_sub1", "bonf_eq0_sub2", "bonf_eq0_sub3", 
                               "bonf_eq1_sub1", "bonf_eq1_sub2", "bonf_eq1_sub3",
                               "bonf_eq2_sub1", "bonf_eq2_sub2", "bonf_eq2_sub3",
                               "bonf_eq3_sub1", "bonf_eq3_sub2", "bonf_eq3_sub3", 
                               "BenH_eq0_sub1", "BenH_eq0_sub2", "BenH_eq0_sub3", 
                               "BenH_eq1_sub1", "BenH_eq1_sub2", "BenH_eq1_sub3",
                               "BenH_eq2_sub1", "BenH_eq2_sub2", "BenH_eq2_sub3",
                               "BenH_eq3_sub1", "BenH_eq3_sub2", "BenH_eq3_sub3")
  POSThoc_type1[num_test, ] <- colMeans(Results[setdiff(1:1000, index_change), 5:28])
  colnames(POSThoc_type1) <- c("bonf_eq0_sub1", "bonf_eq0_sub2", "bonf_eq0_sub3", 
                               "bonf_eq1_sub1", "bonf_eq1_sub2", "bonf_eq1_sub3",
                               "bonf_eq2_sub1", "bonf_eq2_sub2", "bonf_eq2_sub3",
                               "bonf_eq3_sub1", "bonf_eq3_sub2", "bonf_eq3_sub3", 
                               "BenH_eq0_sub1", "BenH_eq0_sub2", "BenH_eq0_sub3", 
                               "BenH_eq1_sub1", "BenH_eq1_sub2", "BenH_eq1_sub3",
                               "BenH_eq2_sub1", "BenH_eq2_sub2", "BenH_eq2_sub3",
                               "BenH_eq3_sub1", "BenH_eq3_sub2", "BenH_eq3_sub3")
  
  print(num_test)
  num_test = num_test + 1
}

save(OMNI_power, OMNI_type1, POSThoc_power, POSThoc_type1,  file = "Sim2results.RData")

########################### summarizing results ################################

#1. Power: test length against carryover effect, when identical items, and 70% vs. 50% of people changes ####

condition[condition$item_character=="parallel" & condition$perc_change==0.7, ]
index_table_70 <- c(1, 13, 25, 2, 14, 26, 3, 15, 27)
condition[index_table_70, ]

condition[condition$item_character=="parallel" & condition$perc_change==0.5, ]
index_table_50 <- c(7, 19, 31, 8, 20, 32, 9, 21, 33)
condition[index_table_50, ]

OMNI_powerTable_parallel <- cbind(OMNI_power[index_table_70, ], OMNI_power[index_table_50, ])
write.csv(OMNI_powerTable_parallel, file = "OMNI_powerTable_parallel.csv")

#2. Power: test length against carryover effect, when non-identical items and 70% people change ####
condition[condition$item_character=="non-parallel" & condition$perc_change==0.7, ]
index_table_70 <- c(4, 16, 28, 5, 17, 29, 6, 18, 30)
condition[index_table_70, ]

condition[condition$item_character=="non-parallel" & condition$perc_change==0.5, ]
index_table_50 <- c(10, 22, 34, 11, 23, 35, 12, 24, 36)
condition[index_table_50, ]
OMNI_powerTable_nonparallel <- cbind(OMNI_power[index_table_70, 1:4], OMNI_power[index_table_50, 1:4])
write.csv(OMNI_powerTable_nonparallel, file = "OMNI_powerTable_nonparallel.csv")

#3. Type1: test length against carryover effect, when identical items, and 70% vs. 50% of people changes ####

condition[condition$item_character=="parallel" & condition$perc_change==0.7, ]
index_table_70 <- c(1, 13, 25, 2, 14, 26, 3, 15, 27)
condition[index_table_70, ]

condition[condition$item_character=="parallel" & condition$perc_change==0.5, ]
index_table_50 <- c(7, 19, 31, 8, 20, 32, 9, 21, 33)
condition[index_table_50, ]

OMNI_type1Table_parallel <- cbind(OMNI_type1[index_table_70, 1:4], OMNI_type1[index_table_50, 1:4])
write.csv(OMNI_type1Table_parallel, file = "OMNI_type1Table_parallel.csv")

#4. Type1: test length against carryover effect, when non-identical items and 70% people change ####
condition[condition$item_character=="non-parallel" & condition$perc_change==0.7, ]
index_table_70 <- c(4, 16, 28, 5, 17, 29, 6, 18, 30)
condition[index_table_70, ]

condition[condition$item_character=="non-parallel" & condition$perc_change==0.5, ]
index_table_50 <- c(10, 22, 34, 11, 23, 35, 12, 24, 36)
condition[index_table_50, ]

OMNI_type1Table_nonparallel <- cbind(OMNI_type1[index_table_70, 1:4], OMNI_type1[index_table_50, 1:4])
write.csv(OMNI_type1Table_nonparallel, file = "OMNI_type1Table_nonparallel.csv")

#########

#5. power: posthoc test, when item length = 40, non parallel items
condition[condition$test_length == 40 & condition$item_character == "non-parallel" & condition$perc_change == .7,]
condition[condition$test_length == 40 & condition$item_character == "non-parallel" & condition$perc_change == .5,]
POSTHOC_table <- rbind(POSThoc_power[condition$test_length == 40 & condition$item_character == "non-parallel" & condition$perc_change == .7, ], 
      POSThoc_power[condition$test_length == 40 & condition$item_character == "non-parallel" & condition$perc_change == .5, ])
write.csv(POSTHOC_table, file = "POSTHOC_table.csv")

condition[condition$test_length == 40 & condition$item_character == "non-parallel" & condition$perc_change == .7,]
condition[condition$test_length == 40 & condition$item_character == "non-parallel" & condition$perc_change == .5,]
type1POSTHOC_table <- rbind(POSThoc_type1[condition$test_length == 40 & condition$item_character == "non-parallel" & condition$perc_change == .7, ], 
                            POSThoc_type1[condition$test_length == 40 & condition$item_character == "non-parallel" & condition$perc_change == .5, ])
write.csv(type1POSTHOC_table, file = "type1POSTHOC_table.csv")

################## END ##########################################################

#######################################################################################################
##############################################################
################## summarizing results #######################
##############################################################
#######################################################################################################
library(ggplot2)
library(reshape2)
library(gridExtra)
library(latex2exp)

categorize_mahalanobis <- function(mydata, percentiles, labels_perc){
  
  mydata <- data.frame(mydata)
  perc_index <- quantile(mydata$mahalanobis, probs = percentiles)
  mydata$category <- cut(mydata$mahalanobis, 
                         breaks=c(-Inf, perc_index, Inf), 
                         labels=labels_perc, 
                         right = TRUE)
  #new_results <- mydata
  new_results <- aggregate(mydata[, 1:(ncol(mydata)-1)], by = list(mydata$category), FUN = median)
  return(new_results)
}

cate_final <- list()
for(i in 1:162){
  cate_final[[i]] <- categorize_mahalanobis(Final_result[[i]], percentiles = c(.3333, .6666), 
                                            labels_perc = c("close", "medium", "far"))
}  #note, 1000 persons per cell. The 1000 persons are grouped into 3 categories. We then compute, for each cell, the median of each category.
# extra comment! cate_final is very important, it is used for generating boxplots (for omnibus) and for generating tables (for posthoc)

box_omnibus <- function(row_index, cate_final, x_title){
  #row_index: this refers to the rows in the matrix "condition". Each row is a cell. Thus, we pick the cells that we want to plot
  #           for example, row_index <- condition$test_length == 5, all cells with test lengths = 5.
  #cate_final:The data matrix obtained by using the function categorize_mahalanobis()
  omni_result <- matrix(NA, 1, 5)
  colnames(omni_result) <- c("Group.1", "omni_eq0", "omni_eq1", "omni_eq2", "omni_eq3") #the names are in line with the names generated by categorize_mahalanobis() 
  for(i in 1:length(row_index)){
    if(row_index[i]){
      omni_result <- rbind(omni_result, cate_final[[i]][, c(1, 6, 7, 8, 9)])
    }
  }
  
  omni_result <- omni_result[-1, ] #the first row is always NA
  colnames(omni_result)[1:5] <- c("category", "Omnibus+Eq(5)", "Omnibus+Eq(11)", "Omnibus+Eq(14)", "Omnibus+Eq(15)")
  
  dat_temp <- melt(omni_result,id.vars="category", measure.vars=c("Omnibus+Eq(5)", "Omnibus+Eq(11)", "Omnibus+Eq(14)", "Omnibus+Eq(15)"))
  dat_temp$category <- factor(dat_temp$category, levels = c("close", "medium", "far"))  #in this way we make sure the correct order of labels
  p <- ggplot(dat_temp, aes(x = category, y = value, fill = category)) +
    geom_boxplot()+
    theme_bw() +
    scale_x_discrete(name = x_title) +
    theme(text = element_text(size = 14,  color="black"),
          axis.text.x=element_blank(), 
          axis.ticks.x = element_blank(),
          axis.text.y=element_text(size = 12, color = "black"), 
          legend.title=element_text(size=12), 
          legend.text=element_text(size=12)) +
    scale_y_continuous(name = "Power", limits = c(0, 1)) +
    facet_grid(. ~ variable) +
    guides(fill=guide_legend(title= expression(paste("Closeness to ", mu [theta][1]))))
  
  return(p)
  
}


# boxplots: test length  #####
row_index <- condition$test_length == 5
box_omnibus(row_index, cate_final, x_title = "Test Length: 5 Items")
row_index <- condition$test_length == 15
box_omnibus(row_index, cate_final, x_title = "Test Length: 15 Items")
row_index <- condition$test_length == 40
box_omnibus(row_index, cate_final, x_title = "Test Length: 40 Items")

# boxplots: parallel/non-parallel
row_index <- condition$item_character == "parallel"
box_omnibus(row_index, cate_final, x_title = "Item Characteristics: Identical Item Parameters")
row_index <- condition$item_character == "non-parallel"
box_omnibus(row_index, cate_final, x_title = "Item Characteristics: Non-Identical Item Parameters")

# boxplots: change: small/large/mixed
row_index <- array()
for(i in 1:nrow(condition)){
  row_index[i] <- ifelse(condition$change_in_theta[[i]][3] == 0.5, TRUE, FALSE)  #small change
}
box_omnibus(row_index, cate_final, x_title = "Small Change")

row_index <- array()
for(i in 1:nrow(condition)){
  row_index[i] <- ifelse(condition$change_in_theta[[i]][3] == 1, TRUE, FALSE)  #large change
}
box_omnibus(row_index, cate_final, x_title = "Large Change")

row_index <- array()
for(i in 1:nrow(condition)){
  row_index[i] <- ifelse(condition$change_in_theta[[i]][3] == 0.9, TRUE, FALSE)  #mixed change
}
box_omnibus(row_index, cate_final, x_title = "Mixed Change")

# boxplots: carry-over effect
row_index <- condition$`carry-over` == 'non'
box_omnibus(row_index, cate_final, x_title = "No Carry-Over Effects")
row_index <- condition$`carry-over` == "30%"
box_omnibus(row_index, cate_final, x_title = "30% of Persons Showing Carry-Over Effects")
row_index <- condition$`carry-over` == "50%"
box_omnibus(row_index, cate_final, x_title = "50% of Persons Showing Carry-Over Effects")

# boxplots: effects correlations among dimensions
row_index <- condition$eff_size_cor_sub_attr == "small"
box_omnibus(row_index, cate_final, x_title = "Small Effect Size of Correlations Among Dimensions")
row_index <- condition$eff_size_cor_sub_attr == "medium"
box_omnibus(row_index, cate_final, x_title = "Medium Effect Size of Correlations Among Dimensions")
row_index <- condition$eff_size_cor_sub_attr == "large"
box_omnibus(row_index, cate_final, x_title = "Large Effect Size of Correlations Among Dimensions")


# table: 40 items, non-identical parameters, mixed change large effect size
condition_number <- 1:nrow(condition)
mixed_change <- array()
for(i in 1:nrow(condition)){
  mixed_change[i] <- ifelse(condition$change_in_theta[[i]][3] == 0.9, TRUE, FALSE)  #mixed change
}
row_index <- (condition$test_length==40) & (condition$item_character=="non-parallel") & (condition$eff_size_cor_sub_attr=="large") & (mixed_change == TRUE) & (condition$`carry-over`=='non')
condition_number[row_index]
part1 <- cate_final[[condition_number[row_index]]]
row_index <- (condition$test_length==40) & (condition$item_character=="non-parallel") & (condition$eff_size_cor_sub_attr=="large") & (mixed_change == TRUE) & (condition$`carry-over`=="30%")
condition_number[row_index]
part2 <- cate_final[[condition_number[row_index]]]
row_index <- (condition$test_length==40) & (condition$item_character=="non-parallel") & (condition$eff_size_cor_sub_attr=="large") & (mixed_change == TRUE) & (condition$`carry-over`=="50%")
condition_number[row_index]
part3 <- cate_final[[condition_number[row_index]]]
part_together <- rbind(part1, part2, part3)
save(part_together, file = "post_hoc_power.RData")
write.csv(part_together, file = "post_hoc_power.csv")


########### old tables for the posthot tests ######################################
sd_median <- function(cate_data, row_index){
  #cate_data: the data matrix where persons are categorized into a few groups. 
  #row_index: this refers to the rows in the matrix "condition". Each row is a cell. Thus, we pick the cells that we want to plot
  #           for example, row_index <- condition$test_length == 5, all cells with test lengths = 5.
  omni_result <- matrix(NA, 1, 25)
  colnames(omni_result) <- c("Group.1", "bonf_eq0_sub1", "bonf_eq0_sub2", "bonf_eq0_sub3", 
                             "bonf_eq1_sub1", "bonf_eq1_sub2", "bonf_eq1_sub3",
                             "bonf_eq2_sub1", "bonf_eq2_sub2", "bonf_eq2_sub3",
                             "bonf_eq3_sub1", "bonf_eq3_sub2", "bonf_eq3_sub3",
                             "BenH_eq0_sub1", "BenH_eq0_sub2", "BenH_eq0_sub3",
                             "BenH_eq1_sub1", "BenH_eq1_sub2", "BenH_eq1_sub3",
                             "BenH_eq2_sub1", "BenH_eq2_sub2", "BenH_eq2_sub3",
                             "BenH_eq3_sub1", "BenH_eq3_sub2", "BenH_eq3_sub3") #the names are in line with the names generated by categorize_mahalanobis() 
  for(i in 1:length(row_index)){
    if(row_index[i]){
      omni_result <- rbind(omni_result, cate_data[[i]][, c(1, 10:33)])
    }
  }
  omni_result <- data.frame(omni_result[-1, ]) #the first row is always NA
  omni_result$Group.1 <- factor(omni_result$Group.1, levels = c("close", "medium", "far"))  #in this way we make sure the correct order of labels
  out_median <- aggregate(omni_result[, 2:ncol(omni_result)], by = list(omni_result$Group.1), FUN = median)
  out_se <- aggregate(omni_result[, 2:ncol(omni_result)], by = list(omni_result$Group.1), FUN = sd)
  out_results <- list(median = out_median, se = out_se)
  
  return(out_results)
}


# table for test length
row_index <- condition$test_length == 5
temp_mat1 <- sd_median(cate_final, row_index)
row_index <- condition$test_length == 15
temp_mat2 <- sd_median(cate_final, row_index)
row_index <- condition$test_length == 40
temp_mat3 <- sd_median(cate_final, row_index)
testlength_mat <- rbind(temp_mat1$median, temp_mat2$median, temp_mat3$median)
test_cond <- c("test length 5", "test length 5", "test length 5",
               "test length 15", "test length 15", "test length 15",
               "test length 40", "test length 40", "test length 40")
testlength_mat <- cbind(test_cond, testlength_mat)

testlength_mat_se <- rbind(temp_mat1$se, temp_mat2$se, temp_mat3$se)
testlength_mat_se <- cbind(test_cond, testlength_mat_se)
save(testlength_mat, testlength_mat_se, file = 'table_testlength.RData')
write.csv(testlength_mat, file = 'table_testlength_median.csv')
write.csv(testlength_mat_se, file = 'table_testlength_se.csv')

# table for item characters
row_index <- condition$item_character == "parallel"
temp_mat1 <- sd_median(cate_final, row_index)
row_index <- condition$item_character == "non-parallel"
temp_mat2 <- sd_median(cate_final, row_index)
itemCharacter_mat <- rbind(temp_mat1$median, temp_mat2$median)
item_cond <- c("parallel", "parallel", "parallel",
               "non-parallel", "non-parallel", "non-parallel")
itemCharacter_mat <- cbind(item_cond, itemCharacter_mat)
itemCharacter_mat_se <- rbind(temp_mat1$se, temp_mat2$se)
itemCharacter_mat_se <- cbind(item_cond, itemCharacter_mat_se)
save(itemCharacter_mat, itemCharacter_mat_se, file = "item_characteristic.RData")
write.csv(itemCharacter_mat, file = "table_itemcharacter_median.csv")
write.csv(itemCharacter_mat_se, file = "table_itemcharacter_se.csv")

# table for change in theta
row_index <- array()
for(i in 1:nrow(condition)){
  row_index[i] <- ifelse(condition$change_in_theta[[i]][3] == 0.5, TRUE, FALSE)  #small change
}
temp_mat1 <- sd_median(cate_final, row_index)
row_index <- array()
for(i in 1:nrow(condition)){
  row_index[i] <- ifelse(condition$change_in_theta[[i]][3] == 1, TRUE, FALSE)  #large change
}
temp_mat2 <- sd_median(cate_final, row_index)
row_index <- array()
for(i in 1:nrow(condition)){
  row_index[i] <- ifelse(condition$change_in_theta[[i]][3] == 0.9, TRUE, FALSE)  #mixed change
}
temp_mat3 <- sd_median(cate_final, row_index)
change_theta_mat <- rbind(temp_mat1$median, temp_mat2$median, temp_mat3$median)
test_cond <- c("Small change", "Small change", "Small change",
               "Large change", "Large change", "Large change",
               "Mixed change", "Mixed change", "Mixed change")
change_theta_mat <- cbind(test_cond, change_theta_mat)
change_theta_mat_se <- rbind(temp_mat1$se, temp_mat2$se, temp_mat3$se)
change_theta_mat_se <- cbind(test_cond, change_theta_mat_se)
save(change_theta_mat, change_theta_mat_se, file = "table_changeintheta.RData")
write.csv(change_theta_mat, file = "table_changeintheta_median.csv")
write.csv(change_theta_mat_se, file = "table_changeintheta_se.csv")

# table for carry-over
row_index <- condition$`carry-over` == 'non'
temp_mat1 <- sd_median(cate_final, row_index)
row_index <- condition$`carry-over` == "30%"
temp_mat2 <- sd_median(cate_final, row_index)
row_index <- condition$`carry-over` == "50%"
temp_mat3 <- sd_median(cate_final, row_index)
carry_mat <- rbind(temp_mat1$median, temp_mat2$median, temp_mat3$median)
test_cond <- c("No carry", "No carry", "No carry",
               "30%", "30%", "30%",
               "50%", "50%", "50%")
carry_mat <- cbind(test_cond, carry_mat)
carry_mat_se <- rbind(temp_mat1$se, temp_mat2$se, temp_mat3$se)
carry_mat_se <- cbind(test_cond, carry_mat_se)

save(carry_mat, carry_mat_se, file = "table_carryover.RData")
write.csv(carry_mat, file = "table_carryover_median.csv")
write.csv(carry_mat_se, file = "table_carryover_se.csv")

# table for corre effect
row_index <- condition$eff_size_cor_sub_attr == "small"
temp_mat1 <- sd_median(cate_final, row_index)
row_index <- condition$eff_size_cor_sub_attr == "medium"
temp_mat2 <- sd_median(cate_final, row_index)
row_index <- condition$eff_size_cor_sub_attr == "large"
temp_mat3 <- sd_median(cate_final, row_index)
cor_effect_mat <- rbind(temp_mat1$median, temp_mat2$median, temp_mat3$median)
test_cond <- c("Small", "Small", "Small",
               "Medium", "Medium", "Medium",
               "Large", "Large", "Large")
cor_effect_mat <- cbind(test_cond, cor_effect_mat)
cor_effect_mat_se <- rbind(temp_mat1$se, temp_mat2$se, temp_mat3$se)
cor_effect_mat_se <- cbind(test_cond, cor_effect_mat_se)

save(cor_effect_mat, cor_effect_mat_se, file = "table_cor_effect.RData")
write.csv(cor_effect_mat, file = "table_coreffect_median.csv")
write.csv(cor_effect_mat_se, file = "table_coreffect_se.csv")
