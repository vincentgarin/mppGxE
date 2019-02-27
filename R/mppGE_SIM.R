#############
# mppGE_SIM #
#############

#' MPP GxE Simple Interval Maping
#'
#' Computes single QTL models along the genome using different models.
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
#' @param plot.gen.eff \code{Logical} value. If \code{plot.gen.eff = TRUE},
#' the function will save the decomposed genetic effects per cross/parent.
#' These results can be ploted with the function \code{\link{plot_genEffects_GE}}
#' to visualize a genome-wide decomposition of the genetic effects.
#' Default value = FALSE.
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
#'
#' @return Return:
#'
#' \item{SIM }{\code{Data.frame} of class \code{QTLprof}. with five columns :
#' 1) QTL marker or in between position names; 2) chromosomes;
#' 3) interger position indicators on the chromosome;
#' 4) positions in centi-Morgan; and 5) -log10(p-val). And if
#' \code{plot.gen.eff = TRUE}, p-values of the cross or parental QTL effects.}
#'
#' @author Vincent Garin
#'
#' @seealso \code{\link{mppGE_CIM}}, \code{\link{plot_genEffects_GE}}
#'
#' @examples
#'
#' # Come later
#'
#' @export
#'


mppGE_SIM <- function(mppData, trait, Q.eff = "cr", VCOV = "CS_CSRT",
                      plot.gen.eff = FALSE, parallel = FALSE, cluster = NULL,
                      workspace = 8e6) {

  # 1. Check data format and arguments
  ####################################

  check_mod_mppGE(mppData = mppData, trait = trait, Q.eff = Q.eff, VCOV = VCOV,
                  QTL_ch = FALSE)



  # 2. Form required elements for the analysis
  ############################################

  nEnv <- length(trait)
  TraitEnv <- c(mppData$pheno[, trait])

  vect.pos <- 1:dim(mppData$map)[1]

  # 3. computation of the SIM profile (genome scan)
  #################################################

  if (parallel) {

    log.pval <- parLapply(cl = cluster, X = vect.pos, fun = QTLModelSIM_GE,
                          mppData = mppData, nEnv = nEnv, TraitEnv = TraitEnv,
                          Q.eff = Q.eff, VCOV = VCOV,
                          plot.gen.eff = plot.gen.eff)

  } else {

    log.pval <- lapply(X = vect.pos, FUN = QTLModelSIM_GE,
                       mppData = mppData, nEnv = nEnv, TraitEnv = TraitEnv,
                       Q.eff = Q.eff, VCOV = VCOV, plot.gen.eff = plot.gen.eff,
                       workspace = workspace)

  }

  log.pval <- t(data.frame(log.pval))
  if(plot.gen.eff & (VCOV == "ID")){log.pval[is.na(log.pval)] <- 1}
  log.pval[, 1] <- check.inf(x = log.pval[, 1]) # check if there are -/+ Inf value
  log.pval[is.na(log.pval[, 1]), 1] <- 0


  # 4. form the results
  #####################

  SIM <- data.frame(mppData$map, log.pval)

  if(plot.gen.eff){

    if(Q.eff == "cr"){

      Env_name <- rep(paste0("_E", 1:nEnv), each = mppData$n.cr)
      Qeff_names <- paste0(rep(unique(mppData$cross.ind), 2), Env_name)

      colnames(SIM)[5:dim(SIM)[2]] <- c("log10pval", Qeff_names)

    } else if (Q.eff == "anc") {

      Env_name <- rep(paste0("_E", 1:nEnv), each = mppData$n.par)
      Qeff_names <- paste0(rep(mppData$parents, 2), Env_name)

      colnames(SIM)[5:dim(SIM)[2]] <- c("log10pval", Qeff_names)

    } else { # par and bi-allelic no modif for gen effects names

      colnames(SIM)[5] <- "log10pval"

    }

  } else {colnames(SIM)[5] <- "log10pval"}


  class(SIM) <- c("QTLprof", "data.frame")

  ### 4.1: Verify the positions for which model could not be computed

  if(sum(SIM$log10pval == 0) > 0) {

    if (sum(SIM$log10pval) == 0){

      message("The computation of the QTL models failled for all positions.
              This could be due to problem in asreml function.")

    } else {

      list.pos <- mppData$map[(SIM$log10pval == 0), 1]

      text <- paste("The computation of the QTL model failed for the following",
                    "positions: ", paste(list.pos, collapse = ", "),
                    ". This could be due to singularities or function issues.")

      message(text)

    }

  }

  return(SIM)

}
