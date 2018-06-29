###############
# QTL_pval_GE #
###############

# function to get the QTL decomposed genetic effect for MPP GxE linear models.

QTL_pval_GE <- function(mppData, nEnv, model, Q.eff, x, par.clu) {
  
  coeffs <- coef(model)
  index <- which(substr(names(coeffs), 1, 3) == "QTL")
  coeffs <- coeffs[index]
  
  var.comp <- sqrt(diag(vcov(model)))
  index <- which(substr(names(var.comp), 1, 3) == "QTL")
  var.comp <- var.comp[index]
  
  var.comp.full <- rep(NA, length(coeffs))
  var.comp.full[match(names(var.comp), names(coeffs))] <- var.comp
  
  pval <- 2 * pt(q = abs(coeffs/var.comp.full),
                 df = df.residual(model), lower.tail = FALSE)
  pval <- pval * sign(coeffs)
  
  names(pval) <- substr(names(pval), 4, nchar(names(pval)))
  
  if (Q.eff == "cr") {
    
    Env_name <- rep(paste0("_E", 1:nEnv), each = mppData$n.cr)
    Cr.name <- paste0(rep(paste0("Cr", unique(mppData$cross.ind)), nEnv),
                      Env_name)
    
    pval <- pval[Cr.name]
    
    
  } else if (Q.eff == "par") {
    
    Env_name <- rep(paste0("_E", 1:nEnv), each = mppData$n.par)
    par.name <- paste0(rep(mppData$parents, nEnv), Env_name)
    
    pval <- pval[par.name]
    
  } else if (Q.eff == "anc") {
    
    ref.all <- paste0("A.allele", par.clu[x, ])
    Env_name <- rep(paste0("_E", 1:nEnv), each = mppData$n.par)
    ref.all <- paste0(ref.all, Env_name)
    
    pval <- pval[ref.all]
    
  }
  
  return(pval)
  
}