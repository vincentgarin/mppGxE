##################
# QTLModelSIM_GE #
##################

# function to compute a single position MPP GxE QTL model

QTLModelSIM_GE <- function(x, mppData, nEnv, TraitEnv, Q.eff, VCOV,
                           plot.gen.eff){

  # 1. formation of the QTL incidence matrix
  ###########################################

  QTL <- inc_mat_QTL(x = x, mppData = mppData, Q.eff = Q.eff, order.MAF = TRUE)

  # QTL main effect

  # QTLmain <- matrix((rep(1, nEnv)), nrow = nEnv) %x% QTL

  # QTL x Env effect

  ref.name <- colnames(QTL)
  Env_name <- rep(paste0("_E", 1:nEnv), each = length(ref.name))
  ref.name <- paste0(c(ref.name, ref.name), Env_name)

  QTLenv <- diag(nEnv) %x% QTL
  colnames(QTLenv) <- ref.name

  QTL.el <- dim(QTLenv)[2] # number of QTL elements

  if(VCOV != "ID"){ # mixed model prepare dataset

    nGeno <- length(mppData$geno.id)

    env_ind <- rep(paste0("Env", 1:nEnv), each = nGeno)
    env_ind <- factor(env_ind, levels = paste0("Env", 1:nEnv))

    cr_ind <- rep(mppData$cross.ind, times = nEnv)
    cr_ind <- factor(cr_ind, levels = unique(mppData$cross.ind))

    geno_id <- rep(mppData$geno.id, times = nEnv)
    geno_id <- factor(x = geno_id, levels = mppData$geno.id)

    dataset <- data.frame(trait = TraitEnv, env = env_ind,
                          genotype = geno_id, cross = cr_ind, QTLenv)

    colnames(dataset)[5:(QTL.el + 4)] <- paste0("Q", 1:QTL.el)

    # model fixed term (QTL) formula

    formula.QTL <- paste("+", paste0("Q", 1:QTL.el), collapse = " ")
    formula.fix <- paste("trait ~ -1 + env:cross", formula.QTL)

  }


  # 2. model computation
  ######################

  ### 2.1 Identity

  if(VCOV == "ID"){

    cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
    CrMatEnv <- diag(nEnv) %x% cross.mat

    model <- tryCatch(expr = lm(TraitEnv ~ -1 + CrMatEnv + QTLenv),
                      error = function(e) NULL)

    if (is.null(model)){

      if(plot.gen.eff) {

        if(Q.eff == "cr"){ results <- c(0, rep(1, nEnv * mppData$n.cr))

        } else if (Q.eff == "biall") { results <- c(0, rep(1, nEnv))

        } else { results <- c(0, rep(1, nEnv * mppData$n.par)) }

      } else { results <- 0 }

    } else {

      if(!("QTLenv" %in% rownames(anova(model)))){ # QTL effect could not be
        # estimated probably due to singularities.

        if(plot.gen.eff) {

          if(Q.eff == "cr"){ results <- c(0, rep(1, nEnv * mppData$n.cr))

          } else if (Q.eff == "biall") { results <- c(0, rep(1, nEnv))

          } else { results <- c(0, rep(1, nEnv * mppData$n.par)) }

        } else { results <- 0 }

      } else {

        if(plot.gen.eff){

          gen.eff <- QTL_pval_GE(mppData = mppData, nEnv = nEnv, model = model,
                                 Q.eff = Q.eff, x = x)

          results <- c(-log10(anova(model)$Pr[2]), gen.eff)

        } else { results <- -log10(anova(model)$Pr[2]) }

      }

    }


  } else if (VCOV == "CS"){

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

  if(VCOV != "ID"){ # compute the mixed model

    if(VCOV %in% c("DG", "UN")){

      model <- tryCatch(expr = asreml(fixed = as.formula(formula.fix),
                                      rcov = as.formula(formula.rcov),
                                      data = dataset, trace = FALSE,
                                      na.method.Y = "include",
                                      na.method.X = "omit", keep.order = TRUE),
                        error = function(e) NULL)

    } else {

      model <- tryCatch(expr = asreml(fixed = as.formula(formula.fix),
                                      random = as.formula(formula.random),
                                      rcov = as.formula(formula.rcov),
                                      data = dataset, trace = FALSE,
                                      na.method.Y = "include",
                                      na.method.X = "omit", keep.order = TRUE),
                        error = function(e) NULL)

    }

    # 3. organise the results for the mixed models similar for all models
    #####################################################################

    if (is.null(model)){

      if(plot.gen.eff) {

        if(Q.eff == "cr"){ results <- c(0, rep(1, nEnv * mppData$n.cr))

        } else if (Q.eff == "biall") { results <- c(0, rep(1, nEnv))

        } else { results <- c(0, rep(1, nEnv * mppData$n.par)) }

      } else { results <- 0 }

    } else {

      W.stat <- sum(wald(model)[2:(QTL.el+1), 3])

      if(W.stat == 0){

        if(plot.gen.eff) {

          if(Q.eff == "cr"){ results <- c(0, rep(1, nEnv * mppData$n.cr))

          } else if (Q.eff == "biall") { results <- c(0, rep(1, nEnv))

          } else { results <- c(0, rep(1, nEnv * mppData$n.par)) }

        } else { results <- 0 }

      } else {

        df <- sum(wald(model)[2:(QTL.el+1), 1])

        pval <- pchisq(W.stat, df, lower.tail = FALSE)

        results <- -log10(pval)

        if(plot.gen.eff){

          gen.eff  <- QTL_pval_mix_GE(model = model, nEnv = nEnv, Q.eff = Q.eff,
                                      QTL.el = QTL.el, x = x, ref.name = ref.name,
                                      par.names = mppData$parents, fct = "SIM")

          results  <- c(results, gen.eff)

        }

      }

    }

  }

  return(results)

}
