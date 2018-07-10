##################
# QTLModelCIM_GE #
##################

# function to compute a single position MPP GxE QTL CIM model

QTLModelCIM_GE <- function(x, mppData, nEnv, TraitEnv, Q.eff, VCOV, cof.list,
                           cof.part, plot.gen.eff){

  # 1. formation of the QTL incidence matrix
  ###########################################

  ### 2.2 QTL position

  QTL <- inc_mat_QTL(x = x, mppData = mppData, Q.eff = Q.eff, order.MAF = TRUE)

  ref.name <- colnames(QTL)
  Env_name <- rep(paste0("_E", 1:nEnv), each = length(ref.name))
  ref.name <- paste0(c(ref.name, ref.name), Env_name)

  QTLenv <- diag(nEnv) %x% QTL
  colnames(QTLenv) <- ref.name

  QTL.el <- dim(QTLenv)[2] # number of QTL elements

  # ref.name <- colnames(QTLenv)

  ### 2.1 cofactors

  cof.mat <- do.call(cbind, cof.list[which(cof.part[x, ])])

  # test if no cofactors

  if(is.null(cof.mat)){ cof.mat <- rep(0, length(TraitEnv)); cof.el <- 1

  } else {

    cof.mat <- diag(nEnv) %x% cof.mat

    cof.el <- dim(cof.mat)[2]

  }

  # prepare mixed model dataset

  if(VCOV != "ID"){

    nGeno <- length(mppData$geno.id)

    env_ind <- rep(paste0("Env", 1:nEnv), each = nGeno)
    env_ind <- factor(env_ind, levels = paste0("Env", 1:nEnv))

    cr_ind <- rep(mppData$cross.ind, times = nEnv)
    cr_ind <- factor(cr_ind, levels = unique(mppData$cross.ind))

    geno_id <- rep(mppData$geno.id, times = nEnv)
    geno_id <- factor(x = geno_id, levels = mppData$geno.id)

    dataset <- data.frame(cof.mat = cof.mat, QTL = QTLenv, trait = TraitEnv,
                          env = env_ind, genotype = geno_id, cross = cr_ind)

    colnames(dataset) <- c(paste0("cof", 1:cof.el), paste0("Q",1:QTL.el),
                           "trait", "env", "genotype", "cross")

    # model fixed term (QTL) formula

    formula.QTL <- paste("+", paste0("Q", 1:QTL.el), collapse = " ")
    formula.fix <- paste("trait ~ -1 + env:cross + grp(cof)", formula.QTL)

  }


  # 2. model computation
  ######################

  ### 2.1 homogeneous residual variance error

  if(VCOV == "ID"){

    cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
    CrMatEnv <- diag(nEnv) %x% cross.mat

    model <- tryCatch(expr = lm(TraitEnv ~ - 1 + CrMatEnv + cof.mat + QTLenv),
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

        results <- -log10(anova(model)$Pr[which(rownames(anova(model))=="QTLenv")])

        if(plot.gen.eff){

          gen.eff <- QTL_pval_GE(mppData = mppData, nEnv = nEnv, model = model,
                                 Q.eff = Q.eff, x = x)


          results <- c(results, gen.eff)

        }

      }

    }

    ### 2.2 HRT REML or cross-specific variance residual terms

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
                                      group = list(cof=1:cof.el),
                                      data = dataset, trace = FALSE,
                                      na.method.Y = "include",
                                      na.method.X = "omit", keep.order = TRUE),
                        error = function(e) NULL)

    } else {

      model <- tryCatch(expr = asreml(fixed = as.formula(formula.fix),
                                      random = as.formula(formula.random),
                                      rcov = as.formula(formula.rcov),
                                      group = list(cof=1:cof.el),
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

      W.stat <- sum(wald(model)[3:(QTL.el+2), 3])

      if(W.stat == 0){

        if(plot.gen.eff) {

          if(Q.eff == "cr"){ results <- c(0, rep(1, nEnv * mppData$n.cr))

          } else if (Q.eff == "biall") { results <- c(0, rep(1, nEnv))

          } else { results <- c(0, rep(1, nEnv * mppData$n.par)) }

        } else { results <- 0 }

      } else {

        df <- sum(wald(model)[3:(QTL.el+2), 1])

        pval <- pchisq(W.stat, df, lower.tail = FALSE)

        results <- -log10(pval)

        if(plot.gen.eff){

          gen.eff  <- QTL_pval_mix_GE(model = model, nEnv = nEnv, Q.eff = Q.eff,
                                      QTL.el = QTL.el, x = x,
                                      ref.name = ref.name,
                                      par.names = mppData$parents, fct = "CIM")

          results  <- c(results, gen.eff)

        }

      }

    }

  }

  return(results)

}
