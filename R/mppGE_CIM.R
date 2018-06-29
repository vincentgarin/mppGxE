###########
# mpp_CIM #
###########

#' MPP GxE composite interval maping
#'
#' Compute QTL models along the genome using cofactors representing other
#' genetic positions for control.
#'
#' For more details about the different models, see documentation of the
#' function \code{\link{mpp_SIM}}. The function returns a -log10(p-value) QTL
#' profile.
#'
#' @param trait \code{Data.frame} with phenotypic trait value. Each column
#' correspond to a different trait or environment.
#'
#' @param mppData An object of class \code{mppData}.
#'
#' @param Q.eff \code{Character} expression indicating the assumption concerning
#' the QTL effects: 1) "cr" for cross-specific; 2) "par" for parental; 3) "anc"
#' for ancestral; 4) "biall" for a bi-allelic. Default = "cr".
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
#' @param cofactors Object of class \code{QTLlist} representing a list of
#' selected position obtained with the function QTL_select() or
#' vector of \code{character} marker or in between marker positions names.
#' Default = NULL.
#'
#' @param window \code{Numeric} distance (cM) on the left and the right of a
#' cofactor position where it is not included in the model. Default = 20.
#'
#' @param plot.gen.eff \code{Logical} value. If \code{plot.gen.eff = TRUE},
#' the function will save the decomposed genetic effects per cross/parent.
#' These results can be ploted with the function \code{\link{plot_genEffects_GE}}
#' to visualize a genome-wide decomposition of the genetic effects.
#' Default value = FALSE.
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
#'
#' @return Return:
#'
#' \item{CIM }{\code{Data.frame} of class \code{QTLprof}. with five columns :
#' 1) QTL marker or in between position names; 2) chromosomes;
#' 3) interger position indicators on the chromosome;
#' 4) positions in centi-Morgan; and 5) -log10(p-val). And if
#' \code{plot.gen.eff = TRUE}, p-values of the cross or parental QTL effects.}
#'
#' @author Vincent Garin
#'
#'
#' @seealso \code{\link{mppGE_SIM}}, \code{\link{plot_genEffects_GE}}
#'
#'
#' @examples
#'
#' # Come later
#'
#' @import ggplot2
#' @import mppR
#'
#' @export
#'

# library(mppR)
#
# source('G:/mppGE/functions/QTLModelSIM_GE.R')
# source('G:/mppGE/functions/QTL_pval_GE.R')
# source('G:/mppGE/functions/CheckModelGE.R')
# source('G:/mppGE/functions/check.inf.R')
# source('G:/mppGE/functions/mppGE_SIM.R')
#
# # Load data
#
# setwd("G:/Data_mppR/EUNAM_Flint")
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
# mppData <- data_bi
# Q.eff = "biall"
# par.clu = par.clu
# VCOV = "h.err"
# plot.gen.eff = TRUE
# parallel = FALSE
# cluster = NULL
#
# # parental
#
# SIM <- mppGE_SIM(trait = trait, mppData = data, Q.eff = "par",
#                  par.clu = par.clu, plot.gen.eff = TRUE)
#
# plot_QTLprof(Qprof = SIM, main = "parental")
#
# cofactors <- QTL_select(Qprof = SIM, threshold = 3, window = 50)
#
# window = 50

mppGE_CIM <- function(trait, mppData, Q.eff = "cr", par.clu = NULL,
                      VCOV = "h.err", cofactors = NULL,  window = 20,
                      plot.gen.eff = FALSE, parallel = FALSE, cluster = NULL)
{

  # 1. Check data format and arguments
  ####################################

  CheckModelGE(mppData = mppData, Q.eff = Q.eff, VCOV = VCOV,
               par.clu = par.clu, plot.gen.eff = plot.gen.eff,
               parallel = parallel, cluster = cluster,
               cofactors = cofactors, fct = "CIM")

  # 2. Form required elements for the analysis
  ############################################

  ### 2.1 cross matrix (cross intercept)

  nEnv <- dim(trait)[2]
  nGeno <- dim(trait)[1]
  TraitEnv <-  c(trait)

  cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
  CrMatEnv <- diag(nEnv) %x% cross.mat

  ### 2.2 parent matrix

  parent.mat <- IncMat_parent(mppData)

  ### 2.4 modify the par.clu object order parents columns and replace monomorphic

  if (Q.eff == "anc") {

    check <- parent_clusterCheck(par.clu = par.clu)
    par.clu <- check$par.clu[, mppData$parents] # order parents columns

  } else {par.clu <- NULL}


  ### 2.5 Formation of the list of cofactors

  if(is.character(cofactors)){

    cof.pos <- which(mppData$map[, 1] %in% cofactors)

  } else {

    cof.pos <- which(mppData$map[, 1] %in% cofactors[, 1])

  }

  cof.list <- lapply(X = cof.pos, FUN = IncMat_QTL, mppData = mppData,
                     cross.mat = cross.mat, par.mat = parent.mat,
                     par.clu = par.clu, Q.eff = Q.eff, order.MAF = TRUE)

  ### 2.6 Formation of the genome-wide and cofactors partition

  vect.pos <- 1:dim(mppData$map)[1]

  # 2.6.1 cofactor partition tells if the cofactor should be included or
  # not in the model at each position.

  if (is.character(cofactors)){

    cofactors2 <- mppData$map[mppData$map[, 1] %in% cofactors, c(2, 4)]

  } else { cofactors2 <- cofactors[, c(2, 4)] }

  test.cof <- function(x, map, window) {

    t1 <- map$chr == as.numeric(x[1])
    t2 <- abs(map$pos.cM - as.numeric(x[2])) < window
    !(t1 & t2)

  }

  cof.part <- apply(X = cofactors2, MARGIN = 1, FUN = test.cof,
                    map = mppData$map, window = window)


  # 3. computation of the CIM profile (genome scan)
  #################################################

  if (parallel) {

    log.pval <- parLapply(cl = cluster, X = vect.pos, fun = QTLModelCIM_GE,
                          mppData = mppData, nEnv = nEnv, TraitEnv = TraitEnv,
                          cross.mat = cross.mat, CrMatEnv = CrMatEnv,
                          par.mat = parent.mat, Q.eff = Q.eff, par.clu = par.clu,
                          VCOV = VCOV, cof.list = cof.list, cof.part = cof.part,
                          plot.gen.eff = plot.gen.eff)

  } else {

    log.pval <- lapply(X = vect.pos, FUN = QTLModelCIM_GE,
                       mppData = mppData, nEnv = nEnv,
                       TraitEnv = TraitEnv, cross.mat = cross.mat,
                       CrMatEnv = CrMatEnv, par.mat = parent.mat, Q.eff = Q.eff,
                       par.clu = par.clu, VCOV = VCOV, cof.list = cof.list,
                       cof.part = cof.part, plot.gen.eff = plot.gen.eff)

  }


  log.pval <- t(data.frame(log.pval))
  if(plot.gen.eff & (VCOV == "h.err")){log.pval[is.na(log.pval)] <- 1}
  log.pval[, 1] <- check.inf(x = log.pval[, 1]) # check if there are -/+ Inf value
  log.pval[is.na(log.pval[, 1]), 1] <- 0

  # 4. form the results
  #####################

  CIM <- data.frame(mppData$map, log.pval)

  if(plot.gen.eff){

    if(Q.eff == "cr"){

      Env_name <- rep(paste0("_E", 1:nEnv), each = mppData$n.cr)
      Qeff_names <- paste0(rep(unique(mppData$cross.ind), 2), Env_name)

      colnames(CIM)[5:dim(CIM)[2]] <- c("log10pval", Qeff_names)

    } else if (Q.eff == "anc") {

      Env_name <- rep(paste0("_E", 1:nEnv), each = mppData$n.par)
      Qeff_names <- paste0(rep(mppData$parents, 2), Env_name)

      colnames(CIM)[5:dim(CIM)[2]] <- c("log10pval", Qeff_names)

    } else { # par and bi-allelic no modif for gen effects names

      colnames(CIM)[5] <- "log10pval"

    }

  } else {colnames(CIM)[5] <- "log10pval"}


  class(CIM) <- c("QTLprof", "data.frame")

  ### 4.1: Verify the positions for which model could not be computed

  if(sum(CIM$log10pval == 0) > 0) {

    if (sum(CIM$log10pval) == 0){

      message("The computation of the QTL models failled for all positions.
              This could be due to problem in asreml function.")

    } else {

      list.pos <- mppData$map[(CIM$log10pval == 0), 1]

      text <- paste("The computation of the QTL model failed for the following",
                    "positions: ", paste(list.pos, collapse = ", "),
                    ". This could be due to singularities or function issues.")

      message(text)

    }

  }

  return(CIM)

}
