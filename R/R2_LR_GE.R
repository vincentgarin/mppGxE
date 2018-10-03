############
# R2_LR_GE #
############

# function to compute likelihood R2 for GxE list of QTLs

R2_LR_GE <- function(mppData, trait, nEnv, QTL, VCOV, workspace){

  nGeno <- length(mppData$geno.id)

  env_ind <- rep(paste0("Env", 1:nEnv), each = nGeno)
  env_ind <- factor(env_ind, levels = paste0("Env", 1:nEnv))

  cr_ind <- rep(mppData$cross.ind, times = nEnv)
  cr_ind <- factor(cr_ind, levels = unique(mppData$cross.ind))

  geno_id <- rep(mppData$geno.id, times = nEnv)
  geno_id <- factor(x = geno_id, levels = mppData$geno.id)

  dataset <- data.frame(QTL = QTL, trait, env = env_ind, genotype = geno_id,
                        cross = cr_ind)

  dataset$cross_env <- factor(paste0(as.character(dataset$cross),
                                     as.character(dataset$env)))

  # n <- length(trait)

  f.full <- "trait ~ -1 + env:cross + grp(QTL)"
  f.red <- "trait ~ -1 + env:cross"

  nQTLel <- dim(QTL)[2]

  if (VCOV == 'CSRT'){

    formula.rcov <- "~ at(cross_env):units"

    model.full <- tryCatch(asreml::asreml(fixed = as.formula(f.full),
                                  rcov = as.formula(formula.rcov),
                                  data = dataset, trace = FALSE,
                                  group = list(QTL = 1:nQTLel),
                                  na.method.Y = "include",
                                  na.method.X = "omit", keep.order = TRUE,
                                  workspace = workspace),
                           error = function(e) NULL)

    model.red <- tryCatch(asreml::asreml(fixed = as.formula(f.red),
                                 rcov = as.formula(formula.rcov),
                                 data = dataset, trace = FALSE,
                                 na.method.Y = "include",
                                 na.method.X = "omit", keep.order = TRUE,
                                 workspace = workspace),
                          error = function(e) NULL)

  } else if (VCOV == "CS_CSRT"){

    formula.random <- "~ genotype"
    formula.rcov <- "~ at(cross_env):units"

    model.full <- tryCatch(asreml::asreml(fixed = as.formula(f.full),
                                  random = as.formula(formula.random),
                                  rcov = as.formula(formula.rcov),
                                  data = dataset, trace = FALSE,
                                  group = list(QTL = 1:nQTLel),
                                  na.method.Y = "include",
                                  na.method.X = "omit", keep.order = TRUE,
                                  workspace = workspace),
                           error = function(e) NULL)

    model.red <- tryCatch(asreml::asreml(fixed = as.formula(f.red),
                                 random = as.formula(formula.random),
                                 rcov = as.formula(formula.rcov),
                                 data = dataset, trace = FALSE,
                                 na.method.Y = "include",
                                 na.method.X = "omit", keep.order = TRUE,
                                 workspace = workspace),
                          error = function(e) NULL)

  }


  if(is.null(model.full) || is.null(model.red)){ R2 <- R2.adj <-  NA

  } else {

    n <- dim(dataset[complete.cases(dataset),])[1]

    R2 <- 1-exp(-(2/n)*(model.full$loglik - model.red$loglik))

    # adjust R2

    z <- dim(QTL)[2]
    N <- model.full$nedf
    R2.adj <- R2 - ((z/N)*(1-R2))

    R2 <- 100 * R2; R2.adj <- 100 * R2.adj


  }

  return(c(R2, R2.adj))

}
