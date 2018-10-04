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

    dataset$genotype[dataset$check == 'check'] <- NA

    dataset <- dataset[order(dataset$cross), ]

    # form the QTL groups

    n.QTL.el <- unlist(lapply(Q.list, function(x) dim(x)[2]))
    QTL.seq <- lapply(n.QTL.el, function(x) seq(1:x))
    cum.sum <- cumsum(n.QTL.el)
    add.el <- c(0, cum.sum[-length(cum.sum)])
    QTL.seq <- mapply(function(x, y) x + y, x = QTL.seq, y = add.el,
                      SIMPLIFY = FALSE)

    # random and rcov formulas

    if (VCOV == 'CSRT'){

      formula.random <- paste0('~ ', exp_des_form)
      formula.rcov <- "~ at(cross_env):units"

    } else if (VCOV == "CS_CSRT"){

      formula.random <- paste0('~ genotype + ', exp_des_form)
      formula.rcov <- "~ at(cross_env):units"

    }

    # compute the model

    model <- tryCatch(asreml(fixed = as.formula(x),
                             random = as.formula(formula.random),
                             rcov = as.formula(formula.rcov), data = dataset,
                             group = QTL.seq,
                             trace = FALSE, na.method.Y = "include",
                             na.method.X = "include",
                             keep.order = TRUE, workspace = workspace),
                      error = function(e) NULL)

   if(!is.null(model)){

     w.table <- wald(model)
     res <- w.table[(dim(w.table)[1] - 1), 4]

   } else {

     message("Problem mixed model computation.")
     res <-  1 # if the model can not be computed we remove the position


     }

  }

  return(res)

}
