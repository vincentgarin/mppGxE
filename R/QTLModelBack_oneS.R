#####################
# QTLModelBack_oneS #
#####################

# single model for MPP GxE QTL backward elimination

QTLModelBack_oneS <- function(x, plot_data, mppData, trait, nEnv, Q.list, VCOV){


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

    # form the QTL groups

    n.QTL.el <- unlist(lapply(Q.list, function(x) dim(x)[2]))
    QTL.seq <- lapply(n.QTL.el, function(x) seq(1:x))
    cum.sum <- cumsum(n.QTL.el)
    add.el <- c(0, cum.sum[-length(cum.sum)])
    QTL.seq <- mapply(function(x, y) x + y, x = QTL.seq, y = add.el,
                      SIMPLIFY = FALSE)

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

    model <- tryCatch(asreml(fixed = as.formula(x),
                             random = as.formula(formula.random),
                             rcov = as.formula(formula.rcov), data = dataset,
                             group = QTL.seq,
                             trace = FALSE, na.method.Y = "include",
                             na.method.X = "omit",
                             keep.order = TRUE),
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
