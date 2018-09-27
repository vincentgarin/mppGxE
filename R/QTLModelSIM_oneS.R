####################
# QTLModelSIM_oneS #
####################

# function to compute single position one stage QTL model

QTLModelSIM_oneS <- function(x, plot_data, mppData, trait, nEnv, EnvNames,
                             Q.eff, VCOV, exp_des_form, plot.gen.eff,
                             workspace){

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

  QTLenv <- data.frame(genotype = mppData$geno.id, QTLenv,
                       stringsAsFactors = FALSE)

  # form the dataset

  ind_row <- split(1:dim(QTLenv)[1], factor(sort(rank(1:dim(QTLenv)[1])%%nEnv)))

  dataset <- c()

  for(i in 1:nEnv){

    data_i <- plot_data[plot_data$env == EnvNames[i], ]
    Q_data_i <- QTLenv[ind_row[[i]], ]
    data_i <- merge(data_i, Q_data_i, by = c("genotype"))

    dataset <- rbind(dataset, data_i)

  }

  dataset$cross_env <- factor(paste0(as.character(dataset$cross),
                                     as.character(dataset$env)))

  # determine mixed model formula

  n_cof <- dim(plot_data)[2]
  colnames(dataset)[(n_cof + 1):(QTL.el + n_cof)] <- paste0("Q", 1:QTL.el)

  formula.QTL <- paste("+", paste0("Q", 1:QTL.el), collapse = " ")
  formula.fix <- paste(trait, " ~ -1 + env:cross", formula.QTL)

  if (VCOV == 'CSRT'){

    formula.random <- paste0('~ ', exp_des_form)
    formula.rcov <- "~ at(cross_env):units"

  } else if (VCOV == "CS_CSRT"){

    formula.random <- paste0('~ genotype + ', exp_des_form)
    formula.rcov <- "~ at(cross_env):units"

  }

  # compute the mixed model

  model <- tryCatch(asreml(fixed = as.formula(formula.fix),
                           random = as.formula(formula.random),
                           rcov = as.formula(formula.rcov), data = dataset,
                           trace = FALSE, na.method.Y = "include",
                           na.method.X = "omit",
                           keep.order = TRUE, workspace = workspace),
                    error = function(e) NULL)

  # Get the results

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

  return(results)

}
