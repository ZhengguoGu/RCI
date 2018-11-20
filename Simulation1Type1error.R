######################################################
#########                                   ##########
#########  Simulation 1: unidimensional     ##########
#########  Type 1 error rate                ##########
#########                                   ##########
######################################################

######################################################
#########  Conditions
######################################################

test_length <- c(5, 15, 40)
item_character <- c("parallel", "non-parallel")
# change_theta = 0  # this is for calculating type 1 error rate

condition <- expand.grid(test_length, item_character)
colnames(condition) <- c("test_length", "item_character")

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
  
  theta <- runif(1000, -3, 3)
  
  
  
}