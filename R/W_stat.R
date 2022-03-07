##########
# W_stat #
##########

# function for the Wald test

W_stat <- function(y, X, Vi){
  
  t(y) %*% Vi %*% X %*% 
    chol2inv(chol(t(X) %*% Vi %*% X)) %*% t(X) %*% Vi %*% y
  
}