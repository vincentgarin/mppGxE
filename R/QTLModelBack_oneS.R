#####################
# QTLModelBack_oneS #
#####################

# single model for MPP GxE QTL backward elimination

QTLModelBack_oneS <- function(x, plot_data, mppData, trait, nEnv, Q.list, VCOV,
                              exp_des_form, workspace){


  if(VCOV == "ID"){ # linear model

    cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
    CrMatEnv <- diag(nEnv) %x% cross.mat

    an.table <- anova(lm(as.formula(x), data = Q.list))
    res <- an.table[(dim(an.table)[1] - 1), 5]

  } else { # mixed models

    # form the dataset

    t_sel <- plot_data[, trait]

    dataset <- data.frame(QTL = do.call(cbind, Q.list), trait = t_sel,
                          plot_data)

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

    # form the QTL groups

    n.QTL.el <- unlist(lapply(Q.list, function(x) dim(x)[2]))
    QTL.seq <- lapply(n.QTL.el, function(x) seq(1:x))
    cum.sum <- cumsum(n.QTL.el)
    add.el <- c(0, cum.sum[-length(cum.sum)])
    QTL.seq <- mapply(function(x, y) x + y, x = QTL.seq, y = add.el,
                      SIMPLIFY = FALSE)

    # random and rcov formulas

    formulas <- mod_formulas_oneS(VCOV = VCOV, exp_des_form = exp_des_form)

    # compute the model

    model <- tryCatch(asreml::asreml(fixed = as.formula(x),
                             random = as.formula(formulas[1]),
                             rcov = as.formula(formulas[2]), data = dataset,
                             group = QTL.seq,
                             trace = FALSE, na.method.Y = "include",
                             na.method.X = "include",
                             keep.order = TRUE, workspace = workspace),
                      error = function(e) NULL)

   if(!is.null(model)){

     w.table <- asreml::wald(model)
     res <- w.table[(dim(w.table)[1] - 1), 4]

   } else {

     message("Problem mixed model computation.")
     res <-  1 # if the model can not be computed we remove the position


     }

  }

  return(res)

}
