#####################
# QTLModelQeff_oneS #
#####################

# Compute a single multi-QTL model to estimate QTL genetic effects in a one
# stage MPP GxE analysis

QTLModelQeff_oneS <- function(plot_data, mppData, trait, Q.list,
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

    dataset$genotype[dataset$check != 'genotype'] <- NA

    if(VCOV %in% c('CSRT', 'CS_CSRT')){

      dataset <- dataset[order(dataset$cross), ]

    } else { # AR1xAR1 VCOVs

      dataset <- dataset[order(dataset$env, dataset$col, dataset$row), ]

      dataset$env <- factor(as.character(dataset$env))
      dataset$col <- factor(as.character(dataset$col))
      dataset$row <- factor(as.character(dataset$row))

    }


    # formula

    f <- paste("trait ~ -1 + check + env:cross +", paste(names.QTL, collapse = "+"))


    # random and rcov formulas

    formulas <- mod_formulas_oneS(VCOV = VCOV, exp_des_form = exp_des_form)

    # compute the model

    model <- tryCatch(asreml(fixed = as.formula(f),
                             random = as.formula(formulas[1]),
                             rcov = as.formula(formulas[2]), data = dataset,
                             trace = FALSE, na.method.Y = "include",
                             na.method.X = "include",
                             keep.order = TRUE, workspace = workspace),
                      error = function(e) NULL)

  }

  return(model)

}
