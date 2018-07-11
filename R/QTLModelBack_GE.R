###################
# QTLModelBack_GE #
###################

# single model for MPP GxE QTL backward elimination

QTLModelBack_GE <- function(x, mppData, trait, nEnv, Q.list, VCOV){


  if(VCOV == "ID"){ # linear model

    cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
    CrMatEnv <- diag(nEnv) %x% cross.mat

    an.table <- anova(lm(as.formula(x), data = Q.list))
    res <- an.table[(dim(an.table)[1] - 1), 5]

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

    # form the QTL groups

    n.QTL.el <- unlist(lapply(Q.list, function(x) dim(x)[2]))
    QTL.seq <- lapply(n.QTL.el, function(x) seq(1:x))
    cum.sum <- cumsum(n.QTL.el)
    add.el <- c(0, cum.sum[-length(cum.sum)])
    QTL.seq <- mapply(function(x, y) x + y, x = QTL.seq, y = add.el,
                      SIMPLIFY = FALSE)

    # random and rcov formulas

    if (VCOV == "CS"){

      formula.random <- "~ genotype"
      formula.rcov <- "~ units"

    } else if (VCOV == "DG"){

      formula.rcov <- "~ at(env):units"


    } else if (VCOV == "UCH"){

      formula.random <- "~ genotype"
      formula.rcov <- "~ at(env):units"


    } else if (VCOV == "UN"){

      formula.rcov <- "~ us(env):genotype"

    }

    # compute the model

    if(VCOV %in% c("DG", "UN")){

      model <- asreml(fixed = as.formula(x),
                      rcov = as.formula(formula.rcov),
                      group = QTL.seq,
                      data = dataset, trace = FALSE,
                      na.method.Y = "include",
                      na.method.X = "omit", keep.order = TRUE)

    } else {

      model <- asreml(fixed = as.formula(x),
                      random = as.formula(formula.random),
                      rcov = as.formula(formula.rcov),
                      group = QTL.seq,
                      data = dataset, trace = FALSE,
                      na.method.Y = "include",
                      na.method.X = "omit", keep.order = TRUE)

    }

    w.table <- wald(model)
    res <- w.table[(dim(w.table)[1] - 1), 4]

  }

  return(res)

}
