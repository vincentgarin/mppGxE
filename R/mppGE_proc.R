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
#' }
#'
#'
#' @param pop.name \code{Character} name of the studied population.
#' Default = "MPP".
#'
#' @param trait.name \code{Character} name of the studied trait.
#' Default = "trait1".
#'
#' @param trait \code{Data.frame} with phenotypic trait value. Each column
#' correspond to a different trait or environment.
#'
#' @param EnvNames \code{character} expression indicating the environment or trait
#' labels. By default it is labeled Env_1, 2, 3, etc.
#'
#' @param mppData An object of class \code{mppData}.
#'
#' @param Q.eff \code{Character} expression indicating the assumption concerning
#' the QTL effect: 1) "cr" for cross-specific effects; 2) "par" parental
#' effects; 3) "anc" for an ancestral effects; 4) "biall" for a bi-allelic
#' effects. Default = "cr".
#'
#' @param par.clu Required argument for the ancesral model \code{(Q.eff = "anc")}.
#' \code{Interger matrix} representing the results of a parents genotypes
#' clustering. The columns represent the parental lines and the rows
#' the different markers or in between positions. \strong{The columns names must
#' be the same as the parents list of the mppData object. The rownames must be
#' the same as the map marker list of the mppData object.} At a particular
#' position, parents with the same value are assumed to inherit from the same
#' ancestor. Default = NULL.
#'
#' @param VCOV \code{Character} expression defining the type of variance
#' covariance structure used: 1) "h.err" for an homogeneous variance residual term
#' (HRT) linear model; 2) "Env.err" for environmental (trait) specific variance.
#' Default = "h.err".
#'
#' @param plot.gen.eff \code{Logical} value. If \code{plot.gen.eff = TRUE},
#' the function will save the decomposed genetic effects per cross/parent.
#' These results can be ploted with the function \code{\link{plot_genEffects_GE}}
#' to visualize a genome-wide decomposition of the genetic effects.
#' \strong{This functionality is ony available for the cross-specific,
#' parental and ancestral models.}
#' Default value = FALSE.
#'
#' @param thre.cof \code{Numeric} value representing the -log10(p-value)
#' threshold above which a position can be peaked as a cofactor. Default = 3.
#'
#' @param win.cof \code{Numeric} value in centi-Morgan representing the minimum
#' distance between two selected cofactors. By default, the function select
#' maximum 1 cofactor per chromosome.
#'
#' @param N.cim \code{Numeric} value specifying the number of time the CIM
#' analysis is repeated. Default = 1.
#'
#' @param window \code{Numeric} distance (cM) on the left and the right of a
#' cofactor position where it is not included in the model. Default = 20.
#'
#' @param thre.QTL \code{Numeric} value representing the -log10(p-value)
#' threshold above which a position can be selected as QTL. Default = 3.
#'
#' @param win.QTL \code{Numeric} value in centi-Morgan representing the minimum
#' distance between two selected QTLs. Default = 20.
#'
#' @param text.size \code{Numeric} value specifying the size of graph axis text
#' elements. Default = 18.
#'
#' @param parallel \code{Logical} value specifying if the function should be
#' executed in parallel on multiple cores. To run function in parallel user must
#' provide cluster in the \code{cluster} argument. \strong{Parallelization is
#' only available for HRT (linear) models \code{VCOV = "h.err"}}.
#' Default = FALSE.
#'
#' @param cluster Cluster object obtained with the function \code{makeCluster()}
#' from the parallel package. Default = NULL.
#'
#' @param verbose \code{Logical} value indicating if the steps of mpp_proc should
#' be printed. It will not affect the printing of the other functions called by
#' \code{mpp_proc()}, especially the printing of \code{asreml()}. Default = TRUE.
#'
#' @param output.loc Path where a folder will be created to save the results.
#' By default the function uses the current working directory.
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

####### example

# argument

# setwd("F:/Data_mppR/EUNAM_Flint")
#
# # pheno
#
# pheno_KWS <- read.csv("./data/pheno/Adj_means_KWS.csv", row.names = 1)
# pheno_CIAM <- read.csv("./data/pheno/Adj_means_CIAM.csv", row.names = 1)
#
# # gather the trait in the same matrix. The number of column correspond to
# # the number of environments.
#
# trait <- cbind(pheno_KWS[, 1], pheno_CIAM[, 1])
# rownames(trait) <- rownames(pheno_CIAM)
# colnames(trait) <- c("Env_1", "Env_2")
#
# # mppData
#
# data_bi <- readRDS(file = "./data/mpp_data/mppData_bi.rds")
# data <- readRDS(file = "./data/mpp_data/mppData.rds")
#
# # par.clu
#
# par.clu <- read.table("./data/clustering/par_clu.txt", h = TRUE,
#                       check.names = FALSE)
# par.clu <- as.matrix(par.clu)
#
# # focus on the first chromosome
#
# mk.sel <- data$map[data$map[, 2]==1, 1]
#
# data <- mppData_subset(mppData = data, mk.list = mk.sel)
# data_bi <- mppData_subset(mppData = data_bi, mk.list = mk.sel)
# par.clu <- par.clu[mk.sel, ]
#
# # other arguments
#
# pop.name = "MPP"
# trait.name = "trait1"
# trait = trait
# EnvNames = c("KWS", "CIAM")
# mppData = data
# Q.eff = "cr"
# par.clu = par.clu
# VCOV = "h.err"
# plot.gen.eff = TRUE
# thre.cof = 3
# win.cof = 20
# N.cim = 1
# window = 20
# thre.QTL = 3
# win.QTL = 20
# text.size = 18
# parallel = FALSE
# cluster = NULL
# verbose = TRUE
# output.loc = "F:/mppGE/results"
#
# source('F:/mppGE/functions/CheckModelGE.R')

mppGE_proc <- function(pop.name = "MPP", trait.name = "trait1", trait,
                       EnvNames = NULL, mppData, Q.eff = "cr", par.clu = NULL,
                       VCOV = "h.err", plot.gen.eff = FALSE, thre.cof = 3,
                       win.cof = NULL, N.cim = 1, window = 20, thre.QTL = 3,
                       win.QTL = 20, text.size = 18, parallel = FALSE,
                       cluster = NULL, verbose = TRUE, output.loc = getwd()) {


  # 1. Check the validity of the parameters that have been introduced
  ###################################################################

  # check.mpp.proc(mppData = mppData, Q.eff = Q.eff, VCOV = VCOV,
  #                par.clu = par.clu, plot.gen.eff = plot.gen.eff,
  #                ref.par = ref.par, parallel = parallel, cluster = cluster,
  #                output.loc = output.loc)

  CheckModelGE(mppData = mppData, trait = trait, Q.eff = Q.eff, VCOV = VCOV,
               par.clu = par.clu, plot.gen.eff = plot.gen.eff,
               parallel = parallel, cluster = cluster, EnvNames = EnvNames,
               fct = "SIM")

  if(is.null(EnvNames)){

    EnvNames <- paste0("Env_", 1:dim(trait)[2])

  }

  if(is.null(win.cof)){

    chr_fact <- factor(mppData$map[, 2], levels = unique(mppData$map[, 2]))

    max_chr_len <- tapply(X = mppData$map[, 4], INDEX = chr_fact,
                          FUN = function(x) max(x))
    max_chr_len <- max(unlist(max_chr_len))

    win.cof <- max_chr_len + 1

  }

  # 2. Create a directory to store the results
  ############################################

  # create a directory to store the results of the QTL analysis

  end.char <- substr(output.loc, nchar(output.loc), nchar(output.loc))

  if(end.char == "/"){

    folder.loc <- paste0(output.loc, paste("QTLGE", pop.name, trait.name, Q.eff,
                                           VCOV, sep = "_"))

  } else {

    folder.loc <- paste0(output.loc, "/", paste("QTLGE", pop.name, trait.name,
                                                Q.eff, VCOV, sep = "_"))

  }

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
                   par.clu = par.clu, VCOV = VCOV, plot.gen.eff = plot.gen.eff,
                   parallel = parallel, cluster = cluster)

  # save SIM results in output location

  write.table(SIM, file = paste0(folder.loc, "/", "SIM.txt"), quote = FALSE,
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
                   par.clu = par.clu, VCOV = VCOV, cofactors = cofactors,
                   window = window, plot.gen.eff = plot.gen.eff,
                   parallel = parallel, cluster = cluster)


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
                       par.clu = par.clu, VCOV = VCOV, cofactors = cofactors,
                       window = window, plot.gen.eff = plot.gen.eff,
                       parallel = parallel, cluster = cluster)

      }

    }

  }

  # save the list of cofactors

  write.table(cofactors[, 1:5], file = paste0(folder.loc, "/", "cofactors.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)

  # save CIM results

  write.table(CIM, file = paste0(folder.loc, "/", "CIM.txt"), quote = FALSE,
              sep = "\t", row.names = FALSE)

  # select QTL candidates

  QTL <- QTL_select(Qprof = CIM, threshold = thre.QTL, window = win.QTL)

  if (is.null(QTL)) { # test if QTL have been selected

    message("No QTL position detected based on the (last) CIM profile.")
    return(NULL)

  }


  # 9. Results processing
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



    pdf(paste0(folder.loc, "/", "QTL_profile.pdf"), height = 10, width = 16)

    print(plot_QTLprof(Qprof = CIM, QTL = cofactors, type = t_plot, main = main.cim,
                       threshold = thre.QTL, text.size = text.size))

    dev.off()

    if (plot.gen.eff) {

      pdf(paste0(folder.loc, "/", "gen_eff.pdf"), height = 10, width = 16)

      print(plot_genEffects_GE(mppData = mppData, nEnv = dim(trait)[2],
                               EnvNames = EnvNames, Qprof = CIM, Q.eff = Q.eff,
                               QTL = QTL, main = main.Qeff, text.size = text.size))

      dev.off()

    }


  ### 9.4: Return R object


  results <- list(n.QTL = dim(QTL)[1], cofactors = cofactors[, 1:5],
                  QTL = QTL[, 1:5])

  return(results)

  }
