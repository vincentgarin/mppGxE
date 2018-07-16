##################
# QTL_effects_GE #
##################

#' MPP GxE QTL genetic effects
#'
#' Compute MPP GxE QTL genetic effects
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
#' covariance structure used. "ID" for identity, "CS" for compound symmetry,
#' "DG" for heterogeneous environmental (residual) variance,
#' "UCH" for uniform covariance with heterogeneous environmental variance,
#' and "UN" for unstructured. Default = "ID".
#'
#' @param QTL Object of class \code{QTLlist} representing a list of
#' selected marker positions obtained with the function QTL_select() or
#' vector of \code{character} marker positions names. Default = NULL.
#'
#' @return Return:
#'
#' \item{Qeff}{\code{List} of \code{data.frame} (one per QTL) containing the
#' following information:
#'
#' \enumerate{
#'
#' \item{QTL genetic effects per cross or parent.}
#' \item{Standard error of the QTL effects.}
#' \item{Test statistics of the effects (t-test or Wald statistic).}
#' \item{P-value of the test statistics.}
#' \item{Significance of the QTL effects.}
#'
#' }
#'
#' }
#'
#' @author Vincent Garin
#'
#' @examples
#'
#' # Come later
#'
#' @export
#'

# arguments

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
# Q.eff = "biall"
# VCOV = "CS"
# QTL <- QTL_select(SIM)

QTL_effects_GE <- function(mppData, trait, Q.eff = "cr", VCOV = "ID",
                         QTL = NULL){

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

  nQTL <- length(Q.pos)

  Q.list <- lapply(X = Q.pos, FUN = inc_mat_QTL, mppData = mppData,
                   Q.eff = Q.eff, order.MAF = TRUE)

  Q.names <- function(x, Q.list, nEnv){
    rep(paste0("Q", x, attr(Q.list[[x]], "dimnames")[[2]]), nEnv)
  }

  names.QTL <- unlist(lapply(X = 1:nQTL, FUN = Q.names, Q.list = Q.list,
                             nEnv = nEnv))

  if(Q.eff == "anc"){

    n_al <- unlist(lapply(X = Q.list, FUN = function(x) dim(x)[2]))

    e_lab <- paste0("E", 1:nEnv)

    Env.names <- lapply(X = n_al, FUN = function(x, e_lab) rep(e_lab, each = x),
                        e_lab = e_lab)

    Env.names <- unlist(Env.names)

  } else {

    n_al <- NULL

    Env.names <- rep(rep(paste0("E", 1:nEnv), each = dim(Q.list[[1]])[2]), nQTL)

  }

  names.QTL <- paste(names.QTL, Env.names, sep = "_")

  Q.list <- lapply(X = Q.list, FUN =  function(x, nEnv) diag(nEnv) %x% x,
                   nEnv = nEnv)

  names(Q.list) <- paste0("Q", 1:length(Q.list))

  # Compute the model

  model <- QTLModelQeff_GE(mppData = mppData, trait = TraitEnv, nEnv = nEnv,
                           Q.list = Q.list, VCOV = VCOV, names.QTL = names.QTL)


  # process the results

  # arguments



  Qeff <- Qeff_res_processing_GE(model = model, mppData = mppData,
                                 Q.eff = Q.eff, VCOV = VCOV,
                                 names.QTL = names.QTL, nQTL = nQTL,
                                 n_al = n_al, nEnv = nEnv)

  return(Qeff)


}
