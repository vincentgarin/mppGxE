#####################
# QTLModelQeff_oneS #
#####################

# Compute a single multi-QTL model to estimate QTL genetic effects in a one
# stage MPP GxE analysis

QTLModelQeff_oneS2 <- function(plot_data, mppData, trait, Q.list,
                              VCOV, exp_des_form, names.QTL, workspace){


  if(VCOV == "ID"){ # linear model

    # Not available

    # cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
    # CrMatEnv <- diag(nEnv) %x% cross.mat
    #
    # an.table <- anova(lm(as.formula(x), data = Q.list))
    # res <- an.table[(dim(an.table)[1] - 1), 5]

  } else { # mixed models

    # form the dataset

    t_sel <- plot_data[, trait]

    dataset <- data.frame(QTL = do.call(cbind, Q.list), trait = t_sel,
                          plot_data)

    colnames(dataset)[1:length(names.QTL)] <- names.QTL

    dataset$cross_env <- factor(paste0(as.character(dataset$cross),
                                       as.character(dataset$env)))

    # formula

    f <- "trait ~ -1 + env:cross + grp(QTLs)"


    # random and rcov formulas

    if (VCOV == 'CSRT'){

      formula.random <- paste0('~ ', exp_des_form)
      formula.rcov <- "~ at(cross_env):units"

    } else if (VCOV == "CS_CSRT"){

      formula.random <- paste0('~ genotype + ', exp_des_form)
      formula.rcov <- "~ at(cross_env):units"

    }

    # compute the model

    model <- tryCatch(asreml::asreml(fixed = as.formula(f),
                             random = as.formula(formula.random),
                             rcov = as.formula(formula.rcov), data = dataset,
                             group = list(QTLs = 1:length(names.QTL)),
                             trace = FALSE, na.method.Y = "include",
                             na.method.X = "omit",
                             keep.order = TRUE, workspace = workspace),
                      error = function(e) NULL)

  }

  return(model)

}
