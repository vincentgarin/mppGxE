##############
# pair_index #
##############

pair_index <- function(x, type = 'whole_mat'){
  
  n_el <- length(x)
  
  x_1 <- c()
  x_2 <- c()
  
  if(type == 'whole_mat'){
    
    for(i in 1:n_el){
      
      x_1 <- c(x_1, rep(i, n_el))
      x_2 <- c(x_2, 1:n_el)
      
    }
    
  } else if (type == 'upper_trg_diag'){
    
    for(i in 1:n_el){
      
      x_1 <- c(x_1, rep(i, (n_el + 1) - i))
      x_2 <- c(x_2, i:n_el)
      
    }
    
  } else if (type == 'upper_trg'){
    
    for(i in 1:n_el){
      
      x_1 <- c(x_1, rep(i, (n_el - i)))
      x_2 <- c(x_2, (i + 1):n_el)
      
    }
    
  }
  
  
  p_i <- cbind(x[x_1], x[x_2])
  
  cov_id <- paste0(p_i[, 1], p_i[, 2])
  
  return(list(p_i = p_i, cov_id = cov_id))
  
}