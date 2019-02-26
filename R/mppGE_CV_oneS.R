#################
# mppGE_CV_oneS #
#################

#' MPP GxE cross-validation one stage analysis
#'
#' MPP GxE cross-validation one stage analysis
#'
#' @param pop.name \code{Character} name of the studied population.
#' Default = "MPP".
#'
#' @param trait.name \code{Character} name of the studied trait.
#' Default = "trait1".
#'
#' @param plot_data \code{data.frame} containing the plot data with the following
#' columns: the trait(s), 'genotype' (genotype indicator), 'cross'
#' (cross indicator), 'env' (environment indicator), and all other experimental
#' design covariates (e.g. replicate, blocks, etc.). The column names of the
#' data.frame must be identical to the one specified ('genotype', 'cross',
#' 'env'). The names of the experimental design covariates must be the same as
#' the one used in 'exp_des_form'.
#'
#' @param mppData An object of class \code{mppData}.
#'
#' @param trait \code{Character} expression for the trait matching the trait
#' column in 'plot_data' argument.
#'
#' @param cv.ref \code{Numerical} or \code{character} indicator to specify which
#' trait of the \code{mppData} object should be used to check the prediction
#' in the CV process. Possibility to specify a vector to predict different
#' within environment genotypic vector values.
#'
#' @param Rep \code{Numeric} value representing the number of repetition of the
#' k-fold procedure. Default = 5.
#'
#' @param k \code{Numeric} value representing the number of folds for the within
#' cross partition of the population. Default = 3.
#'
#' @param EnvNames \code{character} expression indicating the environment or trait
#' labels. By default it is labeled Env_1, 2, 3, etc.
#'
#' @param Q.eff \code{Character} expression indicating the assumption concerning
#' the QTL effect: 1) "cr" for cross-specific effects; 2) "par" parental
#' effects; 3) "anc" for an ancestral effects; 4) "biall" for a bi-allelic
#' effects. Default = "cr".
#'
#' @param VCOV VCOV \code{Character} expression defining the type of variance
#' covariance structure used. "CSRT" for within environment
#' cross-specific residual term, "CS_CSRT" for compound symmetry with within
#' environment cross-specific residual term. Default = "CS_CSRT".
#'
#' @param exp_des_form \code{Character} expression for the random experimental
#' design effects in asreml-R format. For example,
#' 'env:replicate + env:replicate:block'. The column variables names used in
#' 'exp_des_form' should strictly match the names used in 'plot_data'.
#'
#' @param thre.cof \code{Numeric} value representing the -log10(p-value)
#' threshold above which a position can be peaked as a cofactor. Default = 4.
#'
#' @param win.cof \code{Numeric} value in centi-Morgan representing the minimum
#' distance between two selected cofactors. Default = 50 cM.
#'
#' @param N.cim \code{Numeric} value specifying the number of time the CIM
#' analysis is repeated. Default = 1.
#'
#' @param window \code{Numeric} distance (cM) on the left and the right of a
#' cofactor position where it is not included in the model. Default = 20.
#'
#' @param thre.QTL \code{Numeric} value representing the -log10(p-value)
#' threshold above which a position can be selected as QTL. Default = 4.
#'
#' @param win.QTL \code{Numeric} value in centi-Morgan representing the minimum
#' distance between two selected QTLs. Default = 20.
#'
#' @param alpha \code{Numeric} value indicating the level of significance for
#' the backward elimination. Default = 0.01.
#'
#' @param verbose \code{Logical} value indicating if the steps of mpp_proc should
#' be printed. It will not affect the printing of the other functions called by
#' \code{mpp_proc()}, especially the printing of \code{asreml()}. Default = TRUE.
#'
#' @param output.loc Path where a folder will be created to save the results.
#' Default = NULL.
#'
#' @param workspace size of workspace for the REML routines measured in double
#' precision words (groups of 8 bytes). The default is workspace = 8e6.
#'
#'
#' @return Return:
#'
#' \code{List} containing the following results items:
#'
#' \item{p_vs}{\code{Matrix} with : 1) the number of detected QTL;
#' 2) the proportion of predicted genetic variance
#' in the VS (p.vs) at the population level (average of within cross prediction)
#' per environment.}
#'
#' \item{QTL}{\code{Data.frame} containing: 1) the list of QTL position detected
#' at least one time during the entire CV process; 2) the number of times
#' the position has been detected}
#'
#' The results elements return as R object are also saved as text
#' files at the specified output location (\code{output.loc}).
#'
#'
#' @author Vincent Garin
#'
#'
#' @examples
#'
#' # Come later
#'
#'
#' @export
#'

# setwd("F:/Data_mppR/EUNAM_Flint")
#
# pheno <- read.csv("./data/pheno/pheno_red.csv", row.names = 1,
#                   stringsAsFactors = FALSE)
# lines_used <- read.csv("./data/pheno/List_lines_Flint_Lehermeier.csv")
# pheno$Genotype <- as.factor(as.character(pheno$Genotype))
# pheno$Fam <- as.factor(as.character(pheno$Fam))
# pheno <- pheno[order(pheno$Fam), ]
# par_names <- c("D152", "EC49A", "EP44", "EZ5", "F03802", "F2", "F283",
#                "F64", "UH006", "UH009", "DK105")
#
# colnames(pheno)[1] <- "genotype"
# colnames(pheno)[2] <- "cross"
# colnames(pheno)[6] <- "env"
#
# plot_data <- pheno[pheno$env %in% c("KWS", "CIAM"), ]
# load('./data/mpp_data/mppDataGE_toy.RData')
#
# load('./data/mpp_data/mppDataGE_toy.RData')
#
# pop.name = "MPP"
# trait.name = "trait1"
# trait = c("PH")
# Rep = 1
# k = 2
# EnvNames = c("KWS", "CIAM")
# Q.eff = "par"
# VCOV = "DG"
# plot.gen.eff = TRUE
# thre.cof = 4
# win.cof = 50
# N.cim = 1
# window = 20
# thre.QTL = 4
# win.QTL = 20
# alpha = 0.01
# text.size = 18
# verbose = TRUE
# output.loc = "F:/mppGxE/results/toy_example"


mppGE_CV_oneS <- function(pop.name = "MPP", trait.name = "trait1", plot_data,
                          mppData, trait, cv.ref, Rep = 5, k = 3, EnvNames = NULL,
                          Q.eff = "cr", VCOV = "CS_CSRT", exp_des_form,
                          thre.cof = 4, win.cof = 50, N.cim = 1, window = 20,
                          thre.QTL = 4, win.QTL = 20, alpha = 0.01, verbose = TRUE,
                          output.loc = NULL, workspace = 8e6) {


  # 1. Check the validity of the parameters that have been introduced
  ###################################################################

  # CheckModelGE(mppData = mppData, trait = trait, Q.eff = Q.eff, VCOV = VCOV,
  #              plot.gen.eff = plot.gen.eff,
  #              parallel = parallel, cluster = cluster, EnvNames = EnvNames,
  #              fct = "SIM")

  if(is.null(EnvNames)){

    EnvNames <- unique(plot_data$env)

  }


  # 2. Create a directory to store the results
  ############################################

  # create a directory to store the results of the QTL analysis

  folder.loc <- file.path(output.loc, paste("CV_GE_1S", pop.name, trait.name,
                                            Q.eff, VCOV, sep = "_"))

  dir.create(folder.loc)

  # 3. Create space to store the results
  ######################################

  # global results

  nEnv <- length(EnvNames)

  N.QTL <- rep(0, (k*Rep))

  p.ts <- matrix(0, nrow = (k*Rep), ncol = nEnv)
  p.ts.cr <- vector(mode = 'list', length = length(cv.ref))

  for(i in 1:nEnv){

    p.ts.cr[[i]] <- matrix(0, mppData$n.cr, (k*Rep))
    rownames(p.ts.cr[[i]]) <- unique(mppData$cross.ind)
    colnames(p.ts.cr[[i]]) <- paste0('rep', 1:(k*Rep))

  }

  p.vs <- matrix(0, nrow = (k*Rep), ncol = nEnv)
  p.vs.cr <- vector(mode = 'list', length = length(cv.ref))

  for(i in 1:nEnv){

    p.vs.cr[[i]] <- matrix(0, mppData$n.cr, (k*Rep))
    rownames(p.vs.cr[[i]]) <- unique(mppData$cross.ind)
    colnames(p.vs.cr[[i]]) <- paste0('rep', 1:(k*Rep))

  }


  ind.res <- 1 # index to feed the results later

  # individual QTL position results

  QTL.positions <- rep(0, dim(mppData$map)[1])
  # profiles <- c() could be done later

  # keep the marker and in between position full list

  mk.list <- mppData$map[, 1]

  # 4. start to loop from 1 to r replicates
  ########################################

  for (i in 1:Rep) {

    if(verbose){

      cat(paste("CV repetition", i))
      cat("\n")

    }


    ### 4.1 generate a CV partition

    folds <- CV_partition(cross.ind = mppData$cross.ind, k = k)

    ### 4.2 iterate through the CV partition

    for (j in 1:k) {

      if(verbose){

        cat(paste("fold", j))
        cat("\n")

      }

      # training set

      mppData.ts <- subset(x = mppData, gen.list = folds[[j]]$train.set)

      # validation set

      mppData.vs <- subset(x = mppData, gen.list = folds[[j]]$val.set)

      # QTL detection on training set (ts)

      CV_ij <- one_stage_proc(pop.name = paste0("run", ind.res),
                              mppData = mppData.ts, plot_data = plot_data,
                              trait = trait, EnvNames = EnvNames, Q.eff = Q.eff,
                              VCOV = VCOV, exp_des_form = exp_des_form,
                              plot.gen.eff = FALSE,
                              thre.cof = thre.cof, win.cof = win.cof,
                              N.cim = N.cim, window = window, thre.QTL = thre.QTL,
                              win.QTL = win.QTL, alpha = alpha, verbose = FALSE,
                              output.loc = tempdir(), workspace = workspace)

      if (!is.null(CV_ij$QTL)) {

        QTL <- CV_ij$QTL

        ### 4.3 compute the CV statistics (N.QTL, p.ts. etc.)

        # a) N.QTL

        N.QTL[ind.res] <- dim(QTL)[1]

        # b) QTL positions

        QTL.names <- QTL[, 1]

        QTL.positions <- QTL.positions + (is.element(mk.list, QTL.names) * 1)

        # store the profiles (could be done later)

        # compute predicted R squared
        ##############################

        # ts ----

        R2_ts <- QTL_pred_R2_GE_oneS(plot_data = plot_data,
                                     mppData.ts = mppData.ts,
                                     mppData.vs = mppData.ts, trait = trait,
                                     cv.ref = cv.ref, nEnv = nEnv,
                                     Q.eff = Q.eff, VCOV = VCOV,
                                     exp_des_form = exp_des_form, QTL = QTL,
                                     workspace = workspace)


        p.ts[ind.res, ] <- R2_ts$R2_av

        for(w in 1:nEnv){

          p.ts.cr[[w]][, ind.res] <- R2_ts$R2_cr[[w]][unique(mppData$cross.ind)]

        }

        # vs ----

        R2.vs <- QTL_pred_R2_GE_oneS(plot_data = plot_data,
                                     mppData.ts = mppData.ts,
                                     mppData.vs = mppData.vs, trait = trait,
                                     cv.ref = cv.ref, nEnv = nEnv,
                                     Q.eff = Q.eff, VCOV = VCOV,
                                     exp_des_form = exp_des_form, QTL = QTL,
                                     workspace = workspace)


        p.vs[ind.res, ] <- R2.vs$R2_av

        for(w in 1:nEnv){

          p.vs.cr[[w]][, ind.res] <- R2.vs$R2_cr[[w]][unique(mppData$cross.ind)]

        }


      }


      ind.res <- ind.res + 1

    }  # end jth fold loop


  }  # end ith repetition loop

  CV_res <- data.frame(N.QTL, p.ts, p.vs)
  CV_res <- round(CV_res, 2)

  colnames(CV_res) <- c("nQTL", paste0('pts_', EnvNames), paste0('pvs_', EnvNames))

  QTL <- data.frame(mk.names = mk.list, nDet = QTL.positions)

  # save the results

  write.table(x = CV_res, file = file.path(folder.loc, "CV_res.txt"), quote = FALSE,
              sep = "\t", row.names = FALSE)

  save(p.ts.cr, file = file.path(folder.loc, 'pts_cr.RData'))
  save(p.vs.cr, file = file.path(folder.loc, 'pvs_cr.RData'))

  write.table(x = QTL, file = file.path(folder.loc, "QTL.txt"), quote = FALSE,
              sep = "\t", row.names = FALSE)


  results <- list(p.ts = p.ts, p.vs = p.vs, p.vs.cr = p.vs.cr,
                  p.ts.cr = p.ts.cr,QTL = QTL)

  return(results)

}
