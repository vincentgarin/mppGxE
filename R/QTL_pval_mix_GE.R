###################
# QTL_pval_mix_GE #
###################

# Function to collect the p-values from MPP GxE analyses

QTL_pval_mix_GE <- function(model, nEnv, Q.eff, QTL.el, x, par.clu, ref.name,
                            par.names, fct) {
  
  if(fct == "SIM"){
    start.ind <- 2; end.ind <- 1
  } else if(fct == "CIM") {
    start.ind <- 3; end.ind <- 2}
  
  sign <- sign(rev(model$coefficients$fixed[1:QTL.el]))
  pval <- wald(model)[start.ind:(QTL.el + end.ind), 4]
  pval <- pval * sign
  pval[pval == 0] <- 1
  names(pval) <- ref.name
  
  if(Q.eff == "par"){
    
    Env_name <- rep(paste0("_E", 1:nEnv), each = length(par.names))
    par.name <- paste0(rep(par.names, nEnv), Env_name)
    
    pval <- pval[par.name]
    
  } else if (Q.eff == "anc") {
    
    # project into parents
    
    ref.all <- paste0("A.allele", par.clu[x, ])
    Env_name <- rep(paste0("_E", 1:nEnv), each = length(par.names))
    ref.all <- paste0(ref.all, Env_name)
    
    pval <- pval[ref.all]
    
    
  } else if (Q.eff == "biall"){
    
    pval <- pval[seq(1, (2*nEnv), 2)]
    
  }
  
  return(pval)
  
  
}