#############
# mppGE_QTL_R2 #
#############

#' MPP GxE QTL R squared
#'
#' Compute MPP GxE QTL R squared
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
#' @param glb.only \code{Logical} value. If glb.only = TRUE, only the global and
#' global adjusted R squared will be returned. Default = FALSE.
#'
#' @param workspace size of workspace for the REML routines measured in double
#' precision words (groups of 8 bytes). The default is workspace = 8e6.
#'
#'
#' @return Return:
#'
#' \item{}{}
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
# source('H:/PhD/R/package/mppR/R/IncMat_cross.R')
# source('F:/mppGxE/package/mppGxE/R/QTLModelQeff_GE.R')
# source('F:/mppGxE/package/mppGxE/R/sign.star.R')
#
# setwd("F:/Data_mppR/EUNAM_Flint")
#
# load('./data/mpp_data/mppDataGE_toy.RData')
#
# ### SIM
#
# # ID
#
# SIM <- mppGE_SIM(trait = c("PH_KWS", "PH_CIAM"), mppData = mppData,
#                  Q.eff = 'biall', VCOV = 'ID', plot.gen.eff = TRUE)
#
# trait <- c("PH_KWS", "PH_CIAM")
# Q.eff = "par"
# VCOV = "CS_CSRT"
# QTL <- QTL_select(SIM)
# glb.only = FALSE
# workspace = 8e6

mppGE_QTL_R2 <- function(mppData, trait, Q.eff = "cr", VCOV = "ID",
                           QTL = NULL, glb.only = FALSE, workspace = 8e6){

  # Check data format and arguments

  check_mod_mppGE(mppData = mppData, trait = trait, Q.eff = Q.eff, VCOV = VCOV,
                  QTL_ch = TRUE, QTL = QTL)

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

  n.QTL <- length(Q.pos)

  Q.list <- lapply(X = Q.pos, FUN = inc_mat_QTL, mppData = mppData,
                   Q.eff = Q.eff, order.MAF = TRUE)

  Q.list <- lapply(X = Q.list, FUN =  function(x, nEnv) diag(nEnv) %x% x,
                   nEnv = nEnv)

  names(Q.list) <- paste0("Q", 1:length(Q.list))

  # Compute the R2

  if (VCOV == 'ID'){

    R2.all <- R2_lin_GE(mppData = mppData, trait = TraitEnv, nEnv = nEnv,
                     QTL = do.call(cbind, Q.list))

    R2 <- R2.all[[1]]
    R2.adj <- R2.all[[2]]

  } else {


    R2.all <- R2_LR_GE(mppData = mppData, trait = TraitEnv, nEnv = nEnv,
                        QTL = do.call(cbind, Q.list), VCOV = VCOV,
                       workspace = workspace)

    R2 <- R2.all[[1]]
    R2.adj <- R2.all[[2]]

  }

  # Partial R2
  ############

  if(glb.only) {

    QR2Res <- list(glb.R2 = R2, glb.adj.R2 = R2.adj)

  } else {

    if(n.QTL > 1){

      # functions to compute the R squared or all QTL minus 1 or only 1 QTL
      # position

      part.R2.diff.lin <- function(x, QTL, mppData, trait, nEnv) {
        R2_lin_GE(mppData = mppData, trait = trait,
               QTL = do.call(cbind, Q.list[-x]), nEnv = nEnv)
      }

      part.R2.sg.lin <- function(x, QTL, mppData, trait, nEnv) {
        R2_lin_GE(mppData = mppData, trait = trait,
               QTL = do.call(cbind, Q.list[x]), nEnv = nEnv)
      }


      part.R2.diff.LR <- function(x, QTL, mppData, trait, nEnv, VCOV, workspace) {
        R2_LR_GE(mppData = mppData, trait = trait,
                  QTL = do.call(cbind, Q.list[-x]), nEnv = nEnv, VCOV = VCOV,
                 workspace = workspace)
      }

      part.R2.sg.LR <- function(x, QTL, mppData, trait, nEnv, VCOV, workspace) {
        R2_LR_GE(mppData = mppData, trait = trait,
                  QTL = do.call(cbind, Q.list[x]), nEnv = nEnv, VCOV = VCOV,
                 workspace = workspace)
      }


      if(VCOV == 'ID'){

        R2.dif <- lapply(X = 1:n.QTL, FUN = part.R2.diff.lin, QTL = Q.list,
                         mppData = mppData, trait = TraitEnv, nEnv = nEnv)

        R2_i.dif <- lapply(X = R2.dif, FUN = function(x) x[[1]])
        R2_i.dif.adj <- lapply(X = R2.dif, FUN = function(x) x[[2]])

        R2_i.dif <- R2 - unlist(R2_i.dif) # difference full model and model minus i
        R2_i.dif.adj <- R2.adj - unlist(R2_i.dif.adj)

        R2.sg <- lapply(X = 1:n.QTL, FUN = part.R2.sg.lin, QTL = Q.list,
                        mppData = mppData, trait = TraitEnv, nEnv = nEnv)

        R2_i.sg <- unlist(lapply(X = R2.sg, FUN = function(x) x[[1]]))
        R2_i.sg.adj <- unlist(lapply(X = R2.sg, FUN = function(x) x[[2]]))


        names(R2_i.dif) <- names(R2_i.dif.adj) <- paste0("Q", 1:n.QTL)
        names(R2_i.sg) <- names(R2_i.sg.adj) <- paste0("Q", 1:n.QTL)

      } else {

        R2.dif <- lapply(X = 1:n.QTL, FUN = part.R2.diff.LR, QTL = Q.list,
                         mppData = mppData, trait = TraitEnv, nEnv = nEnv,
                         VCOV = VCOV, workspace = workspace)

        R2_i.dif <- lapply(X = R2.dif, FUN = function(x) x[[1]])
        R2_i.dif.adj <- lapply(X = R2.dif, FUN = function(x) x[[2]])

        R2_i.dif <- R2 - unlist(R2_i.dif) # difference full model and model minus i
        R2_i.dif.adj <- R2.adj - unlist(R2_i.dif.adj)

        R2.sg <- lapply(X = 1:n.QTL, FUN = part.R2.sg.LR, QTL = Q.list,
                        mppData = mppData, TraitEnv, nEnv = nEnv, VCOV = VCOV,
                        workspace = workspace)

        R2_i.sg <- unlist(lapply(X = R2.sg, FUN = function(x) x[[1]]))
        R2_i.sg.adj <- unlist(lapply(X = R2.sg, FUN = function(x) x[[2]]))

        names(R2_i.dif) <- names(R2_i.dif.adj) <- paste0("Q", 1:n.QTL)
        names(R2_i.sg) <- names(R2_i.sg.adj) <- paste0("Q", 1:n.QTL)


      }

      QR2Res <- list(glb.R2 = R2,
                     glb.adj.R2 = R2.adj,
                     part.R2.diff = R2_i.dif,
                     part.adj.R2.diff = R2_i.dif.adj,
                     part.R2.sg = R2_i.sg,
                     part.adj.R2.sg = R2_i.sg.adj)

    } else {

      names(R2) <- names(R2.adj) <- "Q1"

      QR2Res <- list(glb.R2 = R2,
                     glb.adj.R2 = R2.adj,
                     part.R2.diff = R2,
                     part.adj.R2.diff = R2.adj,
                     part.R2.sg = R2,
                     part.adj.R2.sg = R2.adj)

    }

  }

  class(QR2Res) <- c("QR2Res", "list")

  return(QR2Res)

}
