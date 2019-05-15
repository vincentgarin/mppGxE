################
# mppGE_back_elim #
################

#' Backward elimination MPP GxE analysis
#'
#' Perform a backward elimination on a list of QTL using a MPP GxE model.
#'
#' @param mppData An object of class \code{mppData}.
#'
#' @param trait \code{Character vector} specifying which traits (environments) should be used.
#'
#' @param Q.eff \code{Character} expression indicating the assumption concerning
#' the QTL effects: 1) "cr" for cross-specific; 2) "par" for parental; 3) "anc"
#' for ancestral; 4) "biall" for a bi-allelic. Default = "cr".
#'
#' @param VCOV VCOV \code{Character} expression defining the type of variance
#' covariance structure used. "ID" for identity, "CSRT" for within environment
#' cross-specific residual terms, "CS_CSRT" for compound symmetry with within
#' environment cross-specific residual terms. Default = "CS_CSRT".
#'
#' @param QTL Object of class \code{QTLlist} representing a list of
#' selected marker positions obtained with the function QTL_select() or
#' a vector of \code{character} marker positions names. Default = NULL.
#'
#' @param alpha \code{Numeric} value indicating the level of significance for
#' the backward elimination. Default = 0.01.
#'
#' @param workspace Size of workspace for the REML routines measured in double
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
#' library(asreml)
#'
#' data(mppData_GE)
#'
#' Qpos <- c("PZE.105068880", "PZE.106098900")
#'
#' QTL <- mppGE_back_elim(mppData = mppData_GE, trait = c('DMY_CIAM', 'DMY_TUM'),
#'                        Q.eff = 'par', QTL = Qpos)
#'
#' @export
#'


mppGE_back_elim <- function(mppData, trait, Q.eff = "cr", VCOV = "CS_CSRT",
                         QTL = NULL, alpha = 0.01, workspace = 8e6){

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

  Q.list <- lapply(X = Q.pos, FUN = inc_mat_QTL, mppData = mppData,
                   Q.eff = Q.eff, order.MAF = TRUE)

  Q.list <- lapply(X = Q.list, FUN =  function(x, nEnv) diag(nEnv) %x% x,
                   nEnv = nEnv)

  names(Q.list) <- paste0("Q", 1:length(Q.list))

  # Compute the model

  ind <- TRUE

  while(ind) {

    ### elaboration of model formulas

    model.formulas <- formula_backward_GE(Q.names = names(Q.list), VCOV = VCOV,
                                          type = 'GE')

    ### computation of the models

    pvals <- lapply(X = model.formulas, FUN = QTLModelBack_GE, mppData = mppData,
                    trait = TraitEnv, nEnv = nEnv, Q.list = Q.list, VCOV = VCOV,
                    workspace = workspace)


    pvals <- unlist(pvals)


    ### test the p-values


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
