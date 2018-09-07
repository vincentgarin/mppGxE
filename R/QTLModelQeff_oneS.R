#####################
# QTLModelQeff_oneS #
#####################

# Compute a single multi-QTL model to estimate QTL genetic effects in a one
# stage MPP GxE analysis

QTLModelQeff_oneS <- function(plot_data, mppData, trait, Q.list,
                              VCOV, names.QTL){


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

    # formula

    f <- paste("trait ~ -1 + env:cross +", paste(names.QTL, collapse = "+"))


    # random and rcov formulas

    if(VCOV == "CS"){

      formula.random <- "~ genotype + env:Rep + env:Rep:Block"
      formula.rcov <- "~ units"

    } else if (VCOV == "DG"){

      formula.random <- "~ env:Rep + env:Rep:Block"
      formula.rcov <- "~ at(env):units"


    } else if (VCOV == "UCH"){

      formula.random <- "~ genotype + env:Rep + env:Rep:Block"
      formula.rcov <- "~ at(env):units"


    } else if (VCOV == "UN"){

      formula.random <- "~ env:Rep + env:Rep:Block"
      formula.rcov <- "~ us(env):genotype"

      # make sure that each genotype appear in each environment

    }

    # compute the model

    model <- tryCatch(asreml(fixed = as.formula(f),
                             random = as.formula(formula.random),
                             rcov = as.formula(formula.rcov), data = dataset,
                             trace = FALSE, na.method.Y = "include",
                             na.method.X = "omit",
                             keep.order = TRUE),
                      error = function(e) NULL)


  }

  return(model)

}
