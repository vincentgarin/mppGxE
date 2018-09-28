########################
# QTLModelQeff_GE_1QTL #
########################

# Function to compute the effect of a single QTL position given a list of QTL
# taken as cofactors. Q_sel indicates for which QTL of Q.list the effects should
# be calculated.

QTLModelQeff_GE_1QTL <- function(mppData, trait, nEnv, Q.list, VCOV, names.QTL,
                                 workspace, Q_sel){

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

    # form the QTL groups

    n.QTL.el <- unlist(lapply(Q.list, function(x) dim(x)[2]))
    QTL.seq <- lapply(n.QTL.el, function(x) seq(1:x))
    cum.sum <- cumsum(n.QTL.el)
    add.el <- c(0, cum.sum[-length(cum.sum)])
    QTL.seq <- mapply(function(x, y) x + y, x = QTL.seq, y = add.el,
                      SIMPLIFY = FALSE)

    cof_id <- 1:nQTL

    cof_id <- cof_id[-Q_sel]

    seq_cof <- unlist(QTL.seq[cof_id])

    seq_QTL <- unlist(QTL.seq[Q_sel])

    # formula

    f <- paste("trait ~ -1 + env:cross + grp(cof) +",
               paste(names.QTL[seq_QTL], collapse = "+"))

    if (VCOV == 'CSRT'){

      formula.rcov <- "~ at(cross_env):units"

      model <- tryCatch(asreml(fixed = as.formula(f),
                               rcov = as.formula(formula.rcov),
                               data = dataset, trace = FALSE,
                               group = list(cof = seq_cof),
                               na.method.Y = "include",
                               na.method.X = "omit", keep.order = TRUE,
                               workspace = workspace),
                        error = function(e) NULL)

    } else if (VCOV == "CS_CSRT"){

      formula.random <- "~ genotype"
      formula.rcov <- "~ at(cross_env):units"

      model <- asreml(fixed = as.formula(f),
                      random = as.formula(formula.random),
                      rcov = as.formula(formula.rcov),
                      group = list(cof = seq_cof),
                      data = dataset, trace = FALSE,
                      na.method.Y = "include",
                      na.method.X = "omit", keep.order = TRUE,
                      workspace = workspace)

    }


  }

  return(list(model = model, names.QTL = names.QTL[seq_QTL]))

}
