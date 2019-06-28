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
cor_subdim <- c("small", "medium", "large")  #effect size of correlation in Sigma_1 (i.e., cov matrix pretest)
cov_mat_sig1 <- list(matrix(c(1, .1, .1, .1, 1, .1, .1, .1, 1), 3, 3), 
                   matrix(c(1, .3, .3, .3, 1, .3, .3, .3, 1), 3, 3), 
                   matrix(c(1, .5, .5, .5, 1, .5, .5, .5, 1), 3, 3))
CO_effect <- c("non", "30%", "50%")  # carry-over effects

condition <- expand.grid(test_length, item_character, CO_effect, cor_subdim)
colnames(condition) <- c("test_length", "item_character", "carry-over", "eff_size_cor_sub_attr")

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
  
  if(condition[num_test, 4] == "small"){
    cov_mat <- cov_mat_sig1[[1]]
  }else if(condition[num_test, 4] == "medium"){
    cov_mat <- cov_mat_sig1[[2]]
  }else{
    cov_mat <- cov_mat_sig1[[3]]
  }
  
  
  theta <- MASS::mvrnorm(1000, mu = c(0, 0, 0), Sigma = cov_mat, empirical = FALSE)
  m_distance <- stats::mahalanobis(theta, 0, cov = cov_mat, inverted = FALSE)
  theta <- theta[order(m_distance),] #accending order in terms of mahalanobis distance
  m_distance_ordered <- m_distance[order(m_distance)]
  
  # some of the persons may show carry-over effects
  if(condition[num_test, 3] == "30%"){  #introducing carry-over effects, if any
    NoCarry_index <- sample(1000, floor(1000 * (1-0.3)), replace = FALSE)  #here we fix the persons who do NOT show carryover
  }else if (condition[num_test, 3] == "50%"){
    NoCarry_index <- sample(1000, floor(1000 * (1-0.5)), replace = FALSE)  #here we fix the persons who do NOT show carryover
  }
  
 
  cl <- makeCluster(6)
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
    #comment: how the data are organized in posthoc_bonf is a bit weird, but this is for later (see posthoc_pvalues and BenHochresults)
    
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
  result <- cbind(theta, m_distance_ordered,  Type1error)
  colnames(result) <- c("theta1", "theta2", "theta3", "mahalanobis", "omni_eq0", "omni_eq1", "omni_eq2", "omni_eq3",
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

save(Final_result, file = "results_simulation2type1error.RData")


###################################################################################
########### boxplots for omnibus tests ############################################
###################################################################################

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
for(i in 1:54){
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
          geom_hline(yintercept=0.1, linetype="dashed", color = "black") +
          scale_y_continuous(name = "Type-I Error Rate", limits = c(0, 1), breaks = c(0, 0.1, .25, .5, .75, 1)) +
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





# table: 40 items, non-identical parameters, large effect size
condition_number <- 1:54
row_index <- (condition$test_length==40) & (condition$item_character=="non-parallel") & (condition$eff_size_cor_sub_attr=="large") & (condition$`carry-over`=='non')
condition_number[row_index]
part1 <- cate_final[[condition_number[row_index]]]
row_index <- (condition$test_length==40) & (condition$item_character=="non-parallel") & (condition$eff_size_cor_sub_attr=="large") & (condition$`carry-over`=="30%")
condition_number[row_index]
part2 <- cate_final[[condition_number[row_index]]]
row_index <- (condition$test_length==40) & (condition$item_character=="non-parallel") & (condition$eff_size_cor_sub_attr=="large") & (condition$`carry-over`=="50%")
condition_number[row_index]
part3 <- cate_final[[condition_number[row_index]]]
part_together <- rbind(part1, part2, part3)
save(part_together, file = "post_hoc_type.RData")
write.csv(part_together, file = "post_hoc_type.csv")

################################# below are old script for creating more tables #################################################

########### tables for the posthot tests ######################################
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
save(testlength_mat, testlength_mat_se, file = 'type1_table_testlength.RData')
write.csv(testlength_mat, file = 'type1_table_testlength_median.csv')
write.csv(testlength_mat_se, file = 'type1_table_testlength_se.csv')

# table for item characteristics
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
save(itemCharacter_mat, itemCharacter_mat_se, file = "type1_item_characteristic.RData")
write.csv(itemCharacter_mat, file = "type1_table_itemcharacter_median.csv")
write.csv(itemCharacter_mat_se, file = "type1_table_itemcharacter_se.csv")

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

save(carry_mat, carry_mat_se, file = "type1_table_carryover.RData")
write.csv(carry_mat, file = "type1_table_carryover_median.csv")
write.csv(carry_mat_se, file = "type1_table_carryover_se.csv")

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

save(cor_effect_mat, cor_effect_mat_se, file = "type1_table_cor_effect.RData")
write.csv(cor_effect_mat, file = "type1_table_coreffect_median.csv")
write.csv(cor_effect_mat_se, file = "type1_table_coreffect_se.csv")
