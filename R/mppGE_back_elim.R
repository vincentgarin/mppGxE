################
# mppGE_back_elim #
################

#' Backward elimination MPP GxE analysis
#'
#' Perform a backward elimination on a list of QTL using a MPP GxE model.
#'
#' @param mppData An object of class \code{mppData}.
#'
#' @param trait \code{Character vector} specifying which traits should be used.
#'
#' @param Q.eff \code{Character} expression indicating the assumption concerning
#' the QTL effects: 1) "cr" for cross-specific; 2) "par" for parental; 3) "anc"
#' for ancestral; 4) "biall" for a bi-allelic. Default = "cr".
#'
#' @param VCOV VCOV \code{Character} expression defining the type of variance
#' covariance structure used. "ID" for identity, "CSRT" for within environment
#' cross-specific residual term, "CS_CSRT" for compound symmetry with within
#' environment cross-specific residual term. Default = "CS_CSRT".
#'
#' @param QTL Object of class \code{QTLlist} representing a list of
#' selected marker positions obtained with the function QTL_select() or
#' vector of \code{character} marker positions names. Default = NULL.
#'
#' @param alpha \code{Numeric} value indicating the level of significance for
#' the backward elimination. Default = 0.01.
#'
#' @param workspace size of workspace for the REML routines measured in double
#' precision words (groups of 8 bytes). The default is workspace = 8e6.
#'
#' @return Return:
#'
#' \item{QTL }{\code{Data.frame} of class \code{QTLlist} with five columns :
#' 1) QTL marker names; 2) chromosomes;
#' 3) interger position indicators on the chromosome;
#' 4) positions in centi-Morgan; and 5) -log10(p-values).}
#'
#' @author Vincent Garin
#'
#' @examples
#'
#' # Come later
#'
#' @export
#'

# # arguments
#
# library(mppR)
# library(asreml)
# library(mppGxE)
#
# # underlying functions
#
# source('H:/PhD/R/package/mppR/R/IncMat_cross.R')
# source('H:/PhD/R/package/Test_function/QTLModelBack_GE.R')
# source('H:/PhD/R/package/Test_function/back_formula_GE.R')
#
# setwd("F:/Data_mppR/EUNAM_Flint")
#
# load('./data/mpp_data/mppDataGE_toy.RData')
#
# SIM <- mppGE_SIM(trait = c("PH_KWS", "PH_CIAM"), mppData = mppData,
#                  Q.eff = 'par', VCOV = 'ID', plot.gen.eff = TRUE)
#
# QTL <- QTL_select(SIM, threshold = 5, window = 1000)
#
# trait <- c("PH_KWS", "PH_CIAM")
# Q.eff <- "cr"
# VCOV <- "UN"
# alpha <- 0.01

mppGE_back_elim <- function(mppData, trait, Q.eff = "cr", VCOV = "CS_CSRT",
                         QTL = NULL, alpha = 0.01, workspace = 8e6){

  if(is.null(QTL)){stop("No 'QTL' have been provided.")}

  # form the trait value

  nEnv <- length(trait)
  TraitEnv <- c(mppData$pheno[, trait])


  # form the list of QTLs

  if(is.character(QTL)){

    Q.pos <- which(mppData$map[, 1] %in% QTL)

    QTL <- mppData$map[mppData$map[, 1] %in% QTL, ]

  } else {

    Q.pos <- which(mppData$map[, 1] %in% QTL[, 1])

  }

  Q.list <- lapply(X = Q.pos, FUN = inc_mat_QTL, mppData = mppData,
                   Q.eff = Q.eff, order.MAF = TRUE)

  Q.list <- lapply(X = Q.list, FUN =  function(x, nEnv) diag(nEnv) %x% x,
                   nEnv = nEnv)

  names(Q.list) <- paste0("Q", 1:length(Q.list))

  # Compute the model

  ind <- TRUE

  while(ind) {

    ### 3.1 elaboration of model formulas

    model.formulas <- formula_backward_GE(Q.names = names(Q.list), VCOV = VCOV,
                                          type = 'GE')

    ### 3.2 computation of the models

    pvals <- lapply(X = model.formulas, FUN = QTLModelBack_GE, mppData = mppData,
                    trait = TraitEnv, nEnv = nEnv, Q.list = Q.list, VCOV = VCOV,
                    workspace = workspace)


    pvals <- unlist(pvals)


    ### 3.4 test the p-values


    if(sum(pvals > alpha) > 0) {

      # remove the QTL position from the list

      Q.list <- Q.list[!(pvals==max(pvals))]

      # test if there is no more positions

      if(length(Q.list) == 0){ind <- FALSE}

    } else {

      # stop the procedure

      ind <- FALSE

    }

  }

  QTL <- QTL[as.numeric(substr(names(Q.list), 2, nchar(names(Q.list)))), ,
             drop = FALSE]

  if(dim(QTL)[1] == 0){

    QTL <- NULL

    return(QTL)

  } else{

    class(QTL) <- c("QTLlist", "data.frame")

    return(QTL)

  }

}