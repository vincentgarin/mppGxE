##############
# mppGE_proc #
##############

#' MPP GxE QTL analysis
#'
#' The procedure is the following:
#'
#' \enumerate{
#'
#' \item{Simple interval mapping (SIM) to select cofactor
#' (\code{\link{mppGE_SIM}}).}
#'
#' \item{Composite interval mapping (CIM) with selected cofactors
#' (\code{\link{mppGE_CIM}}).}
#'
#' \item{Backward elimination on the list of QTLs
#' (\code{\link{mppGE_back_elim}}).}
#'
#' }
#'
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
#' @param EnvNames \code{character} expression indicating the environment or trait
#' labels. By default it is labeled Env_1, 2, 3, etc.
#'
#' @param Q.eff \code{Character} expression indicating the assumption concerning
#' the QTL effect: 1) "cr" for cross-specific effects; 2) "par" parental
#' effects; 3) "anc" for an ancestral effects; 4) "biall" for a bi-allelic
#' effects. Default = "cr".
#'
#' @param VCOV VCOV \code{Character} expression defining the type of variance
#' covariance structure used. "ID" for identity, "CSRT" for within environment
#' cross-specific residual term, "CS_CSRT" for compound symmetry with within
#' environment cross-specific residual term. Default = "CS_CSRT".
#'
#' @param plot.gen.eff \code{Logical} value. If \code{plot.gen.eff = TRUE},
#' the function will save the decomposed genetic effects per cross/parent.
#' These results can be ploted with the function \code{\link{plot_genEffects_GE}}
#' to visualize a genome-wide decomposition of the genetic effects.
#' Default value = FALSE.
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
#' @param text.size \code{Numeric} value specifying the size of graph axis text
#' elements. Default = 18.
#'
#' @param parallel \code{Logical} value specifying if the function should be
#' executed in parallel on multiple cores. To run function in parallel user must
#' provide cluster in the \code{cluster} argument. \strong{Parallelization is
#' only available for 'ID' model}. Default = FALSE.
#'
#' @param cluster Cluster object obtained with the function \code{makeCluster()}
#' from the parallel package. Default = NULL.
#'
#' @param workspace size of workspace for the REML routines measured in double
#' precision words (groups of 8 bytes). The default is workspace = 8e6.
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
#' List containing the following items:
#'
#' \item{n.QTL}{Number of detected QTLs.}
#'
#' \item{cofactors}{\code{Data.frame} with cofactors positions.}
#'
#' \item{QTL}{\code{Data.frame} with QTL positions.}
#'
#' \item{R2}{\code{List} containing R squared statistics of the QTL effects.
#' For details see \code{\link{QTL_R2}} output section.}
#'
#' Some output files are also saved at the specified location
#' (\code{output.loc}):
#'
#' \enumerate{
#'
#' \item{The SIM and CIM results in a text file (SIM.txt, CIM.txt).}
#'
#' \item{The list of cofactors (cofactors.txt).}
#'
#' \item{The list of QTL (QTL.txt).}
#'
#' \item{The QTL R squared statistics (QTL_R2.txt) (for details see
#' \code{\link{QTL_R2}}).}
#'
#' \item{The plot of the CIM profile (QTL_profile.pdf) with dotted vertical
#' lines representing the cofactors positions. If \code{plot.gen.eff = TRUE},
#' plot of the genetic effects per cross or parents obtained with
#' \code{\link{plot_genEffects_GE}}  (gen_eff.pdf) with dashed
#' lines representing the QTL positions.}
#'
#' }
#'
#'
#' @author Vincent Garin
#'
#' @seealso
#'
#' \code{\link{mppGE_CIM}},
#' \code{\link{mppGE_SIM}},
#' \code{\link{plot_genEffects_GE}},
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
# EnvNames = c("KWS", "CIAM")
# Q.eff = "cr"
# VCOV = "ID"
# plot.gen.eff = TRUE
# thre.cof = 4
# win.cof = 50
# N.cim = 1
# window = 20
# thre.QTL = 4
# win.QTL = 20
# text.size = 18
# parallel = FALSE
# cluster = NULL
# verbose = TRUE
# output.loc = "F:/mppGxE/results/toy_example"


mppGE_proc <- function(pop.name = "MPP", trait.name = "trait1", mppData, trait,
                       EnvNames = NULL,  Q.eff = "cr", VCOV = "CS_CSRT",
                       plot.gen.eff = FALSE, thre.cof = 4,
                       win.cof = 50, N.cim = 1, window = 20, thre.QTL = 4,
                       win.QTL = 20, alpha = 0.01, text.size = 18, parallel = FALSE,
                       cluster = NULL, workspace = 8e6, verbose = TRUE,
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

  # if(is.null(win.cof)){
  #
  #   chr_fact <- factor(mppData$map[, 2], levels = unique(mppData$map[, 2]))
  #
  #   max_chr_len <- tapply(X = mppData$map[, 4], INDEX = chr_fact,
  #                         FUN = function(x) max(x))
  #   max_chr_len <- max(unlist(max_chr_len))
  #
  #   win.cof <- max_chr_len + 1
  #
  # }

  # 2. Create a directory to store the results
  ############################################

  # create a directory to store the results of the QTL analysis

  folder.loc <- file.path(output.loc, paste("QTLGE", pop.name, trait.name, Q.eff,
                                            VCOV, sep = "_"))

  dir.create(folder.loc)


  # 3. Cofactors selection - SIM
  ##############################

  if(verbose){

    cat("\n")
    cat("Cofactors selection - SIM")
    cat("\n")
    cat("\n")

  }

  SIM <- mppGE_SIM(trait = trait, mppData = mppData, Q.eff = Q.eff,
                   VCOV = VCOV, plot.gen.eff = plot.gen.eff,
                   parallel = parallel, cluster = cluster,
                   workspace = workspace)

  # save SIM results in output location

  write.table(SIM, file = file.path(folder.loc, "SIM.txt"), quote = FALSE,
              sep = "\t", row.names = FALSE)

  # cofactors selection

  cofactors <- QTL_select(Qprof = SIM, threshold = thre.cof, window = win.cof)


  if (is.null(cofactors)) { # test if cofactors have been selected

    message("No QTL/cofactor position detected based on the SIM profile.")

    return(NULL)

  }

  # 4. Multi-QTL model search - CIM
  #################################

  if(verbose){

    cat("\n")
    cat("Multi-QTL model search - CIM")
    cat("\n")
    cat("\n")

  }

  CIM <- mppGE_CIM(trait = trait, mppData = mppData, Q.eff = Q.eff,
                   VCOV = VCOV, cofactors = cofactors,
                   window = window, plot.gen.eff = plot.gen.eff,
                   parallel = parallel, cluster = cluster,
                   workspace = workspace)


  if (N.cim > 1) {

    for (i in 1:(N.cim - 1)) {

      # take the cofactors of the previous analysis

      cofactors <- QTL_select(Qprof = CIM, threshold = thre.cof,
                              window = win.cof)

      if (is.null(cofactors)) { # test if cofactors have been selected

        mes.text <- paste("No cofactor position detected in CIM profile nb", i)
        message(mes.text)

        return(NULL)

      } else {

        if(verbose){

          cat("\n")
          cat(paste("CIM scan", (i+1)))
          cat("\n")
          cat("\n")

        }

        CIM <- mppGE_CIM(trait = trait, mppData = mppData, Q.eff = Q.eff,
                         VCOV = VCOV, cofactors = cofactors,
                         window = window, plot.gen.eff = plot.gen.eff,
                         parallel = parallel, cluster = cluster,
                         workspace = workspace)

      }

    }

  }

  # save the list of cofactors

  write.table(cofactors[, 1:5], file = file.path(folder.loc, "cofactors.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)

  # save CIM results

  write.table(CIM, file = file.path(folder.loc, "CIM.txt"), quote = FALSE,
              sep = "\t", row.names = FALSE)

  # select QTL candidates

  QTL <- QTL_select(Qprof = CIM, threshold = thre.QTL, window = win.QTL)

  if (is.null(QTL)) { # test if QTL have been selected

    message("No QTL position detected based on the (last) CIM profile.")
    return(NULL)

  }

  # 5. Backward elimination
  #########################

  Q_back <- mppGE_back_elim(mppData = mppData, trait = trait, Q.eff = Q.eff,
                      VCOV = VCOV, QTL = QTL, alpha = alpha,
                      workspace = workspace)

  # save the list of QTLs

  if(!is.null(Q_back)){

    QTL <- QTL[QTL[, 1] %in% Q_back[, 1], ]

    write.table(QTL[, 1:5], file = file.path(folder.loc, "QTL.txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)

  }

  # 6. QTL R2
  ###########

  R2 <- mppGE_QTL_R2(mppData = mppData, trait = trait, Q.eff = Q.eff, VCOV = VCOV,
                  QTL = QTL[, 1], workspace = workspace)

  # save R2 results

  QTL.R2 <- data.frame(QTL[, 1:5], round(R2[[3]], 2), round(R2[[4]], 2),
                       round(R2[[5]], 2), round(R2[[6]], 2),
                       stringsAsFactors = FALSE)

  colnames(QTL.R2)[6:9] <- c("R2.diff", "adj.R2.diff", "R2.sg", "adj.R2.sg")

  write.table(QTL.R2, file = file.path(folder.loc, "QTL_R2.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)


  # 7. Results processing
  #######################

  if(verbose){

    cat("\n")
    cat("Results processing")
    cat("\n")
    cat("\n")

  }


  ### 9.2: Plots


  main.cim <- paste("CIM", pop.name, trait.name, Q.eff, VCOV)
  main.Qeff <- paste("QTL gen. effects", pop.name, trait.name, Q.eff, VCOV)

  if(Q.eff == "biall"){t_plot <- "h"} else {t_plot <- "l"}

  pdf(file.path(folder.loc, "QTL_profile.pdf"), height = 10, width = 16)

  print(plot(x = CIM, QTL = cofactors, type = t_plot, main = main.cim,
             threshold = thre.QTL, text.size = text.size))

  dev.off()

  if (plot.gen.eff) {

    pdf(file.path(folder.loc, "gen_eff.pdf"), height = 10, width = 16)

    print(plot_genEffects_GE(mppData = mppData, nEnv = length(trait),
                             EnvNames = EnvNames, Qprof = CIM, Q.eff = Q.eff,
                             QTL = QTL, main = main.Qeff, text.size = text.size))

    dev.off()

  }


  ### 9.4: Return R object


  results <- list(n.QTL = dim(QTL)[1], cofactors = cofactors[, 1:5],
                  QTL = QTL[, 1:5], R2 = R2)

  return(results)

}
