###################
# QTLModelBeta_GE #
###################

# function to compute the genetic effects (Beta) of a list of QTL in a MPP GxE
# model.

QTLModelBeta_GE <- function(mppData, trait, nEnv, Q.list, VCOV, names.QTL,
                            workspace){

  nQTL <- length(Q.list)

  # form the dataset

  nGeno <- length(mppData$geno.id)

  env_ind <- rep(paste0("Env", 1:nEnv), each = nGeno)
  env_ind <- factor(env_ind, levels = paste0("Env", 1:nEnv))

  cr_ind <- rep(mppData$cross.ind, times = nEnv)
  cr_ind <- factor(cr_ind, levels = unique(mppData$cross.ind))

  geno_id <- rep(mppData$geno.id, times = nEnv)
  geno_id <- factor(x = geno_id, levels = mppData$geno.id)

  dataset <- data.frame(QTL = do.call(cbind, Q.list), trait,
                        env = env_ind, genotype = geno_id, cross = cr_ind)

  dataset$cross_env <- factor(paste0(as.character(dataset$cross),
                                     as.character(dataset$env)))

  colnames(dataset)[1:length(names.QTL)] <- names.QTL

  # formula

  f <- "trait ~ -1 + env:cross + grp(QTLs)"

  if (VCOV == 'CSRT'){

    formula.random <- "~ NULL"
    formula.rcov <- "~ at(cross_env):units"

  } else if (VCOV == "CS_CSRT"){

    formula.random <- "~ genotype"
    formula.rcov <- "~ at(cross_env):units"

  }

  model <- asreml::asreml(fixed = as.formula(f),
                  random = as.formula(formula.random),
                  rcov = as.formula(formula.rcov),
                  group = list(QTLs = 1:length(names.QTL)),
                  data = dataset, trace = FALSE,
                  na.method.Y = "include",
                  na.method.X = "include", keep.order = TRUE,
                  workspace = workspace)

  return(model)

}
