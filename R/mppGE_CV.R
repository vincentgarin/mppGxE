############
# mppGE_CV #
############

#' MPP GxE cross-validation
#'
#' MPP GxE cross-validation
#'
#' @param pop.name \code{Character} name of the studied population.
#' Default = "MPP".
#'
#' @param trait.name \code{Character} name of the studied trait.
#' Default = "trait1".
#'
#' @param mppData An object of class \code{mppData}.
#'
#' @param trait \code{Character vector} specifying which traits should be used.
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
#' covariance structure used. "ID" for identity, "CS" for compound symmetry,
#' "DG" for heterogeneous environmental (residual) variance,
#' "UCH" for uniform covariance with heterogeneous environmental variance,
#' and "UN" for unstructured. Default = "ID".
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
#' @param parallel \code{Logical} value specifying if the function should be
#' executed in parallel on multiple cores. To run function in parallel user must
#' provide cluster in the \code{cluster} argument. \strong{Parallelization is
#' only available for 'ID' model}. Default = FALSE.
#'
#' @param cluster Cluster object obtained with the function \code{makeCluster()}
#' from the parallel package. Default = NULL.
#'
#' @param verbose \code{Logical} value indicating if the steps of mpp_proc should
#' be printed. It will not affect the printing of the other functions called by
#' \code{mpp_proc()}, especially the printing of \code{asreml()}. Default = TRUE.
#'
#' @param output.loc Path where a folder will be created to save the results.
#' Default = NULL.
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
#' @seealso
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
# load('./data/mpp_data/mppDataGE_toy.RData')
#
# pop.name = "MPP"
# trait.name = "trait1"
# trait = c("PH_KWS", "PH_CIAM")
# Rep = 1
# k = 3
# EnvNames = c("KWS", "CIAM")
# Q.eff = "par"
# VCOV = "ID"
# plot.gen.eff = TRUE
# thre.cof = 4
# win.cof = 50
# N.cim = 1
# window = 20
# thre.QTL = 4
# win.QTL = 20
# alpha = 0.01
# text.size = 18
# parallel = FALSE
# cluster = NULL
# verbose = TRUE
# output.loc = "F:/mppGxE/results/toy_example"


mppGE_CV <- function(pop.name = "MPP", trait.name = "trait1", mppData, trait,
                     Rep = 5, k = 3, EnvNames = NULL,  Q.eff = "cr", VCOV = "ID",
                     thre.cof = 4, win.cof = 50, N.cim = 1, window = 20,
                     thre.QTL = 4, win.QTL = 20, alpha = 0.01,
                     parallel = FALSE, cluster = NULL, verbose = TRUE,
                     output.loc = NULL) {


  # 1. Check the validity of the parameters that have been introduced
  ###################################################################

  # CheckModelGE(mppData = mppData, trait = trait, Q.eff = Q.eff, VCOV = VCOV,
  #              plot.gen.eff = plot.gen.eff,
  #              parallel = parallel, cluster = cluster, EnvNames = EnvNames,
  #              fct = "SIM")

  if(is.null(EnvNames)){

    EnvNames <- paste0("Env_", 1:dim(trait)[2])

  }


  # 2. Create a directory to store the results
  ############################################

  # create a directory to store the results of the QTL analysis

  folder.loc <- file.path(output.loc, paste("CV_GE", pop.name, trait.name, Q.eff,
                                            VCOV, sep = "_"))

  dir.create(folder.loc)

  # 3. Create space to store the results
  ######################################

  # global results

  nEnv <- length(trait)

  N.QTL <- rep(0, (k*Rep))
  p.vs <- matrix(0, nrow = (k*Rep), ncol = nEnv)

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

      CV_ij <- mppGE_proc(pop.name = paste0("run", ind.res), mppData = mppData.ts,
                          trait = trait, EnvNames = EnvNames, Q.eff = Q.eff,
                          VCOV = VCOV, plot.gen.eff = FALSE,
                          thre.cof = thre.cof, win.cof = win.cof,
                          N.cim = N.cim, window = window, thre.QTL = thre.QTL,
                          win.QTL = win.QTL, alpha = alpha, parallel = parallel,
                          cluster = cluster, verbose = FALSE,
                        output.loc = tempdir())

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

          R2.vs <- QTL_pred_R2_GE(mppData.ts = mppData.ts,
                                  mppData.vs = mppData.vs, trait = trait,
                                  Q.eff = Q.eff, VCOV = VCOV, QTL = QTL)


          # global results
          ################

          # non adjusted

          p.vs[ind.res, ] <- R2.vs


        }


      ind.res <- ind.res + 1

    }  # end jth fold loop


  }  # end ith repetition loop

  p_vs <- cbind(N.QTL, p.vs)

  colnames(p_vs) <- c("nQTL", EnvNames)
  rownames(p_vs) <- paste0("CV_run_", 1:(Rep * k))

  QTL <- data.frame(mk.names = mk.list, nDet = QTL.positions)

  # save the results

  write.table(x = p_vs, file = file.path(folder.loc, "p_vs.txt"), quote = FALSE,
              sep = "\t", row.names = FALSE)

  write.table(x = QTL, file = file.path(folder.loc, "QTL.txt"), quote = FALSE,
              sep = "\t", row.names = FALSE)


  results <- list(p_vs = p_vs, QTL = QTL)

  return(results)

}