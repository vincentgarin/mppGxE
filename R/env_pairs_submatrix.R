#######################
# env_pairs_submatrix #
#######################

# function to create the sub-matrix elements used to construct the
# unstructured matrix

env_pairs_submatrix <- function(nEnv){
  
  env_id <- paste0('E', 1:nEnv)
  env_gr <- expand.grid(env_id, env_id)
  env_pairs <- paste0(env_gr[, 2], env_gr[, 1])
  
  e_p_m <- matrix(env_pairs, nEnv)
  e_p_m[lower.tri(e_p_m)] <-  t(e_p_m)[lower.tri(e_p_m)]
  all_e_p <- c(e_p_m)
  
  un_e_p <- e_p_m[lower.tri(e_p_m)]
  
  n_un_p <- length(un_e_p)
  m_list <- vector(mode = 'list', length = n_un_p)
  
  for(i in 1:n_un_p){
    
    m_list[[i]] <- matrix((all_e_p %in% un_e_p[i]*1), nEnv)
    
  }
  
  names(m_list) <- un_e_p
  
  return(m_list)
  
}