###########
# MM_comp #
###########

# Function to calculate the different mixed models without the QTL effect

# the initial mixed model is:

# y = (Env + cr) + g*e + err

MM_comp <- function(mppData, nEnv, y, VCOV, maxIter = 100, msMaxIter = 100){
  
  nGeno <- dim(mppData$pheno)[1]
  env <- rep(paste0('E', 1:nEnv), each = nGeno)
  cross <- rep(mppData$cross.ind, nEnv)
  geno <- rep(rownames(mppData$pheno), nEnv)
  cross_env <- paste0(cross, '_', env)
  
  d <- data.frame(trait = y, env = env, cross_env = cross_env, geno = geno)
  d[, 2:4] <- lapply(d[, 2:4], as.factor)
  
  ### CS model
  if (VCOV == 'CS'){
    
    m <- tryCatch(lme(trait ~ cross_env, random = ~ 1 | geno,
                  control = list(opt = "optim", maxIter = maxIter,
                             msMaxIter = msMaxIter),
                  data = d, na.action = na.omit), error = function(x) NULL)
    
  ### CSE model
  } else if (VCOV == 'CSE'){
    
    m <- tryCatch(gls(trait ~ cross_env,
                  weights = varIdent(form = ~ 1 | cross_env),
                  control = list(opt = "optim", maxIter = maxIter,
                            msMaxIter = msMaxIter),
                  data = d, na.action = na.omit), error = function(x) NULL)
    
  ### CS + CSE model  
  } else if (VCOV == 'CS_CSE'){
    
    m <- tryCatch(lme(trait ~ cross_env, random = ~ 1 | geno,
                      weights = varIdent(form = ~ 1 | cross_env),
                      control = list(opt = "optim", maxIter = maxIter,
                                     msMaxIter = msMaxIter),
                      data = d, na.action = na.omit), error = function(x) NULL)
    
  ### Unstructured model  
  } else if ((VCOV == 'UN') | (VCOV == 'UN_K')){
    
    
    m <- tryCatch(lme(trait ~ cross_env,
                      random = list(geno = pdSymm(form = ~ -1 + env)),
                      control = list(opt = "optim", maxIter = maxIter,
                                     msMaxIter = msMaxIter),
                      data = d, na.action = na.omit), error = function(x) NULL)
    
  }
  
  return(list(model = m, data = d))
  
}