##############
# R2_LR_oneS #
##############

R2_LR_oneS <- function(plot_data, mppData, trait, nEnv, QTL, VCOV,
                       exp_des_form, workspace){

  t_sel <- plot_data[, trait]

  dataset <- data.frame(QTL = QTL, trait = t_sel, plot_data)

  dataset$cross_env <- factor(paste0(as.character(dataset$cross),
                                     as.character(dataset$env)))


  f.full <- "trait ~ -1 + env:cross + grp(QTL)"
  f.red <- "trait ~ -1 + env:cross"

  nQTLel <- dim(QTL)[2]

  if (VCOV == 'CSRT'){

    formula.random <- paste0('~ ', exp_des_form)
    formula.rcov <- "~ at(cross_env):units"

  } else if (VCOV == "CS_CSRT"){

    formula.random <- paste0('~ genotype + ', exp_des_form)
    formula.rcov <- "~ at(cross_env):units"

  }


  model.full <- tryCatch(asreml::asreml(fixed = as.formula(f.full),
                                        random = as.formula(formula.random),
                                        rcov = as.formula(formula.rcov),
                                        data = dataset, trace = FALSE,
                                        group = list(QTL = 1:nQTLel),
                                        na.method.Y = "include",
                                        na.method.X = "omit", keep.order = TRUE,
                                        workspace = workspace),
                         error = function(e) NULL)

  model.red <- tryCatch(asreml::asreml(fixed = as.formula(f.red),
                                       random = as.formula(formula.random),
                                       rcov = as.formula(formula.rcov),
                                       data = dataset, trace = FALSE,
                                       na.method.Y = "include",
                                       na.method.X = "omit", keep.order = TRUE,
                                       workspace = workspace),
                        error = function(e) NULL)


  if(is.null(model.full) || is.null(model.red)){ R2 <- R2.adj <-  NA

  } else {

    n <- dim(dataset[complete.cases(dataset),])[1]

    R2 <- 1-exp(-(2/n)*(model.full$loglik - model.red$loglik))

    # adjust R2

    z <- dim(QTL)[2]
    N <- model.full$nedf
    R2.adj <- R2 - ((z/N)*(1-R2))

    R2 <- 100 * R2; R2.adj <- 100 * R2.adj


  }

  return(c(R2, R2.adj))

}
