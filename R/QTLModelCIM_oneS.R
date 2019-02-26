####################
# QTLModelCIM_oneS #
####################

# function to compute single position one stage QTL model with cofactors (CIM)

QTLModelCIM_oneS <- function(x, plot_data, mppData, trait, nEnv, EnvNames,
                             Q.eff, VCOV, exp_des_form, cof.list, cof.part,
                             plot.gen.eff, workspace){

  # process the QTL incidence matrix

  QTL <- inc_mat_QTL(x = x, mppData = mppData, Q.eff = Q.eff,
                     order.MAF = TRUE)

  rownames(QTL) <- mppData$geno.id

  ref.name <- colnames(QTL)
  Env_name <- rep(paste0("_E", 1:nEnv), each = length(ref.name))
  ref.name <- paste0(c(ref.name, ref.name), Env_name)

  QTLenv <- diag(nEnv) %x% QTL
  colnames(QTLenv) <- ref.name

  QTL.el <- dim(QTLenv)[2]

  # process the cofactors elements

  cof.mat <- do.call(cbind, cof.list[which(cof.part[x, ])])

  # test if no cofactors

  if(is.null(cof.mat)){ cof.mat <- rep(0, length(mppData$geno.id)); cof.el <- 1

  } else {

    cof.mat <- diag(nEnv) %x% cof.mat

    cof.el <- dim(cof.mat)[2]

  }

  # Form a unique dataset with cofactors and QTL to distribute it on the data

  QTLdat <- data.frame(genotype = rep(mppData$geno.id, nEnv), cof.mat, QTLenv,
                       stringsAsFactors = FALSE)

  # form the dataset

  ind_row <- split(1:dim(QTLdat)[1], factor(sort(rank(1:dim(QTLdat)[1])%%nEnv)))

  dataset <- c()

  for(i in 1:nEnv){

    data_i <- plot_data[plot_data$env == EnvNames[i], ]
    Q_data_i <- QTLdat[ind_row[[i]], ]
    data_i <- merge(data_i, Q_data_i, by = c("genotype"), all.x = TRUE)

    dataset <- rbind(dataset, data_i)

  }



  # add colnames cofactors and QTL

  n_cof <- dim(plot_data)[2]
  colnames(dataset) <- c(colnames(dataset)[1:n_cof], paste0("cof", 1:cof.el),
                         paste0("Q",1:QTL.el))

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

  # determine mixed model formula

  formula.QTL <- paste("+", paste0("Q", 1:QTL.el), collapse = " ")
  formula.fix <- paste(trait, "~ -1 + env:check + env:cross + grp(cof)", formula.QTL)

  # random and rcov formulas

  formulas <- mod_formulas_oneS(VCOV = VCOV, exp_des_form = exp_des_form)

  # compute the mixed model

  model <- tryCatch(asreml::asreml(fixed = as.formula(formula.fix),
                           random = as.formula(formulas[1]),
                           rcov = as.formula(formulas[2]), data = dataset,
                           group = list(cof = (n_cof + 1):(n_cof + cof.el)),
                           trace = FALSE, na.method.Y = "include",
                           na.method.X = "include",
                           keep.order = TRUE,
                           workspace = workspace),
                    error = function(e) NULL)

  # Get the results

  if (is.null(model)){

    if(plot.gen.eff) {

      if(Q.eff == "cr"){ results <- c(0, rep(1, nEnv * mppData$n.cr))

      } else if (Q.eff == "biall") { results <- c(0, rep(1, nEnv))

      } else { results <- c(0, rep(1, nEnv * mppData$n.par)) }

    } else { results <- 0 }

  } else {

    W.stat <- sum(asreml::wald(model)[4:(QTL.el+3), 3])

    if(W.stat == 0){

      if(plot.gen.eff) {

        if(Q.eff == "cr"){ results <- c(0, rep(1, nEnv * mppData$n.cr))

        } else if (Q.eff == "biall") { results <- c(0, rep(1, nEnv))

        } else { results <- c(0, rep(1, nEnv * mppData$n.par)) }

      } else { results <- 0 }

    } else {

      df <- sum(asreml::wald(model)[4:(QTL.el+3), 1])

      pval <- pchisq(W.stat, df, lower.tail = FALSE)

      results <- -log10(pval)

      if(plot.gen.eff){

        gen.eff  <- QTL_pval_mix_GE(mppData = mppData, model = model,
                                    nEnv = nEnv, Q.eff = Q.eff,
                                    QTL.el = QTL.el, x = x, ref.name = ref.name,
                                    par.names = mppData$parents, fct = "CIM",
                                    mod = 'M4')

        results  <- c(results, gen.eff)

      }

    }

  }

  return(results)

}
