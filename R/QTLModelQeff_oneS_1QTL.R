##########################
# QTLModelQeff_oneS_1QTL #
##########################

# Compute a single multi-QTL model to estimate QTL genetic effects in a one
# stage MPP GxE analysis

QTLModelQeff_oneS_1QTL <- function(plot_data, mppData, trait, Q.list,
                              VCOV, exp_des_form, names.QTL, workspace,
                              Q_sel){


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

    dataset <- dataset[order(dataset$cross), ]

    #########################

    nQTL <- length(Q.list)

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

    f <- paste("trait ~ -1 + check + env:cross + grp(cof) +",
               paste(names.QTL[seq_QTL], collapse = "+"))

    ####################################

    # random and rcov formulas

    if (VCOV == 'CSRT'){

      formula.random <- paste0('~ ', exp_des_form)
      formula.rcov <- "~ at(cross_env):units"

    } else if (VCOV == "CS_CSRT"){

      formula.random <- paste0('~ genotype + ', exp_des_form)
      formula.rcov <- "~ at(cross_env):units"

    }

    # compute the model

    model <- tryCatch(asreml(fixed = as.formula(f),
                             random = as.formula(formula.random),
                             rcov = as.formula(formula.rcov), data = dataset,
                             group = list(cof = seq_cof),
                             trace = FALSE, na.method.Y = "include",
                             na.method.X = "include",
                             keep.order = TRUE, workspace = workspace),
                      error = function(e) NULL)

  }

  return(list(model = model, names.QTL = names.QTL[seq_QTL]))

}
