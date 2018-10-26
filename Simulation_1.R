###############################################################
#########                Simulation 1                 #########
###############################################################

set.seed(110)

#------------------------------------------------------------------------------ 
# parameteres 
#------------------------------------------------------------------------------

###################################################################################
test_length <- c(5, 15, 40)
parallel_item <- c(1, 0) # 1== yes, parallel, 0 == no
theta_change <- c(0, .5, 1) # 0 == no change, .5 == small change, 1 == large change

num_condition <- length(test_length) * length(parallel_item) * length(theta_change)
conditions <- list()

p <- 1
for(i in 1:length(test_length)){
  for (j in 1:length(parallel_item)){
    for (k in 1:length(theta_change)){
          
          conditions[[p]] <- c(test_length[i], parallel_item[j], theta_change[k])
          p <- p + 1
          
    }
  }
}


df <- data.frame(matrix(unlist(conditions), nrow=num_condition, byrow = T))  # now we have in total 18 cells
colnames(df) <- c('test length', 'parallel item', 'theta change')

###############################################################################

num_test <- 1
while (num_test <= nrow(df)){
  
  num_items <- df[num_test, 1]
  parallel_items <- df[num_test, 2]
  change_of_theta <- df[num_test, 3]
  
  # Step 1: simulating item parameters
  
  if (parallel_items == 1) {
    
    itempar <- matrix(NA,num_items,5)
    itempar[,1] <- runif(1,1.5,2.5)   # discrimination
    avg_beta <- runif(1, 0, 1.25)
    itempar[,2] <- avg_beta - .75
    itempar[,3] <- avg_beta - .25
    itempar[,4] <- avg_beta + .25
    itempar[,5] <- avg_beta + .75
    
  } else {
    
    itempar <- matrix(NA,num_items,5)
    itempar[,1] <- runif(num_items,1.5,2.5)  # discrimination
    avg_beta <- runif(num_items, 0, 1.25)
    itempar[,2] <- avg_beta - .75
    itempar[,3] <- avg_beta - .25
    itempar[,4] <- avg_beta + .25
    itempar[,5] <- avg_beta + .75
    
  }
  
  
}
