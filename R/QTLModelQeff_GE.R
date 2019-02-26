###################
# QTLModelQeff_GE #
###################

QTLModelQeff_GE <- function(mppData, trait, nEnv, Q.list, VCOV, names.QTL,
                            workspace){

  nQTL <- length(Q.list)

  if(VCOV == "ID"){ # linear model

    # form the formula

    Qenv <- do.call(cbind, Q.list)
    colnames(Qenv) <- names.QTL

    cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
    CrMatEnv <- diag(nEnv) %x% cross.mat

    model <- lm(trait ~ - 1 + CrMatEnv + Qenv)

  } else { # mixed models

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

    f <- paste("trait ~ -1 + env:cross +", paste(names.QTL, collapse = "+"))

    if (VCOV == 'CSRT'){

      formula.rcov <- "~ at(cross_env):units"

      model <- asreml::asreml(fixed = as.formula(f),
                      rcov = as.formula(formula.rcov),
                      data = dataset, trace = FALSE,
                      na.method.Y = "include",
                      na.method.X = "omit", keep.order = TRUE,
                      workspace = workspace)

    } else if (VCOV == "CS_CSRT"){

      formula.random <- "~ genotype"
      formula.rcov <- "~ at(cross_env):units"

      model <- asreml::asreml(fixed = as.formula(f),
                      random = as.formula(formula.random),
                      rcov = as.formula(formula.rcov),
                      data = dataset, trace = FALSE,
                      na.method.Y = "include",
                      na.method.X = "omit", keep.order = TRUE,
                      workspace = workspace)

    }


  }

  return(model)

}
