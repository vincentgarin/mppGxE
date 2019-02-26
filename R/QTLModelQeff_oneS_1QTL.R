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

    if(VCOV %in% c('CSRT', 'CS_CSRT')){

      dataset <- dataset[order(dataset$cross), ]

    } else { # AR1xAR1 VCOVs

      dataset <- dataset[order(dataset$env, dataset$col, dataset$row), ]

      dataset$env <- factor(as.character(dataset$env))
      dataset$col <- factor(as.character(dataset$col))
      dataset$row <- factor(as.character(dataset$row))

    }

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

    f <- paste("trait ~ -1 + env:check + env:cross + grp(cof) +",
               paste(names.QTL[seq_QTL], collapse = "+"))

    ####################################

    # random and rcov formulas

    # random and rcov formulas

    formulas <- mod_formulas_oneS(VCOV = VCOV, exp_des_form = exp_des_form)

    # compute the model

    model <- tryCatch(asreml::asreml(fixed = as.formula(f),
                             random = as.formula(formulas[1]),
                             rcov = as.formula(formulas[2]), data = dataset,
                             group = list(cof = seq_cof),
                             trace = FALSE, na.method.Y = "include",
                             na.method.X = "include",
                             keep.order = TRUE, workspace = workspace),
                      error = function(e) NULL)

  }

  return(list(model = model, names.QTL = names.QTL[seq_QTL]))

}
