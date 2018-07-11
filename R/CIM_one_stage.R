#################
# CIM_one_stage #
#################

#' Perform a MPP GxE one stage CIM analysis
#'
#' @param plot_data \code{data.frame} containing plot data with the following
#' columns: genotype indicator, cross, environment, experimental design
#' covariate (e.g. replicate, block, etc.), and the traits.
#'
#' @param mppData Object of class \code{mppData} contaning the genotypic
#' information with genotype list corresponding to the one of \code{plot_data}.
#'
#' @param trait \code{Character} expression for the trait.
#'
#' @param Q.eff \code{Character} expression indicating the assumption concerning
#' the QTL effects: 1) "cr" for cross-specific; 2) "par" for parental; 3) "anc"
#' for ancestral; 4) "biall" for a bi-allelic. Default = "cr".
#'
#' @param VCOV VCOV \code{Character} expression defining the type of variance
#' covariance structure used. "CS" for compound symmetry, "DG" for heterogeneous
#' environmental (residual) variance, "UCH" for uniform covariance with
#' heterogeneous environmental variance, and "UN" for unstructured.
#' Default = "DG".
#'
#' @param cofactors Object of class \code{QTLlist} representing a list of
#' selected marker positions obtained with the function QTL_select() or
#' vector of \code{character} marker positions names.
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
#' @seealso \code{\link{plot_genEffects_GE}}
#'
#' @examples
#'
#' # Come later
#'
#' @export
#'

# arguments

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
# trait = "PH"
# Q.eff = "par"
# VCOV = "DG"
# plot.gen.eff = FALSE
#
# SIM <- SIM_one_stage(plot_data = plot_data, mppData = mppData,
#                      trait = "PH", Q.eff = "par", VCOV = "CS",
#                      plot.gen.eff = TRUE)
#
# cofactors <- QTL_select(SIM, threshold = 5, window = 1000)
#
# window <- 20

CIM_one_stage <- function(plot_data, mppData, trait, Q.eff = "cr", VCOV = "DG",
                          cofactors = NULL, window = 20, plot.gen.eff = FALSE){

  if(VCOV == "UN"){stop("This VCOV is not available for the moment.")}

  # 1. Remove the genotype of plot data that do not have genotypic information

  plot_data <- plot_data[plot_data$genotype %in% mppData$geno.id, ]

  # 2. Determine the environments

  EnvNames <- unique(plot_data$env)

  nEnv <- length(EnvNames)

  # 3. form list of cofactors

  if(is.character(cofactors)){

    cof.pos <- which(mppData$map[, 1] %in% cofactors)

  } else {

    cof.pos <- which(mppData$map[, 1] %in% cofactors[, 1])

  }

  cof.list <- lapply(X = cof.pos, FUN = inc_mat_QTL, mppData = mppData,
                     Q.eff = Q.eff, order.MAF = TRUE)

  # form the cofactors partition

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

  vect.pos <- 1:dim(mppData$map)[1]

  ############# form QTLModelCIM_oneS

  log.pval <- lapply(X = vect.pos, FUN = QTLModelCIM_oneS, plot_data = plot_data,
                     mppData = mppData, trait = trait, nEnv = nEnv,
                     EnvNames = EnvNames, Q.eff = Q.eff, VCOV = VCOV,
                     cof.list = cof.list, cof.part = cof.part,
                     plot.gen.eff = plot.gen.eff)

  ############

  log.pval <- t(data.frame(log.pval))
  # if(plot.gen.eff & (VCOV == "ID")){log.pval[is.na(log.pval)] <- 1}
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
