###############
# QTL_R2_oneS #
###############

#' MPP GxE one stage QTL R2
#'
#' Compute MPP GxE one stage QTL R2
#'
#' @param plot_data \code{data.frame} containing the plot data with the following
#' columns: the trait(s), 'genotype' (genotype indicator), 'cross'
#' (cross indicator), 'env' (environment indicator), and all other experimental
#' design covariates (e.g. replicate, blocks, etc.). The column names of the
#' data.frame must be identical to the one specified ('genotype', 'cross',
#' 'env'). The names of the experimental design covariates must be the same as
#' the one used in 'exp_des_form'.
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
#' @param exp_des_form \code{Character} expression for the random experimental
#' design effects in asreml-R format. For example,
#' 'env:replicate + env:replicate:block'. The column variables names used in
#' 'exp_des_form' should strictly match the names used in 'plot_data'.
#'
#' @param QTL Object of class \code{QTLlist} representing a list of
#' selected marker positions obtained with the function QTL_select() or
#' vector of \code{character} marker positions names. Default = NULL.
#'
#' @param glb.only \code{Logical} value. if \code{glb.only = TRUE}, the
#' function only calculate the global R2 for all QTLs. Default = FALSE.
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

# source('F:/mppGxE/package/mppGxE/R/formula_backward_GE.R')
# source('F:/mppGxE/package/mppGxE/R/QTLModelBack_oneS.R')
#
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
# VCOV = "CS_CSRT"
# alpha = 0.01
# exp_des_form = 'env:Rep + env:Rep:Block'
# workspace = 8e6
#
# SIM <- SIM_one_stage(plot_data = plot_data, mppData = mppData,
#                      trait = "PH", Q.eff = "par", VCOV = "CS_CSRT",
#                      plot.gen.eff = TRUE,
#                      exp_des_form = 'env:Rep + env:Rep:Block')
#
# QTL <- QTL_select(SIM, threshold = 4, window = 1000)


QTL_R2_oneS <- function(plot_data, mppData, trait, Q.eff = "cr", VCOV = "ID",
                        exp_des_form, QTL = NULL, glb.only = FALSE,
                        workspace = 8e6){

  if(is.null(QTL)){stop("No 'QTL' have been provided.")}

  plot_data <- plot_data[plot_data$genotype %in% mppData$geno.id, ]

  # Determine the environments

  EnvNames <- unique(plot_data$env)

  nEnv <- length(EnvNames)

  # form the list of QTLs

  if(is.character(QTL)){

    Q.pos <- which(mppData$map[, 1] %in% QTL)

    QTL <- mppData$map[mppData$map[, 1] %in% QTL, ]

  } else {

    Q.pos <- which(mppData$map[, 1] %in% QTL[, 1])

  }

  nQTL <- length(Q.pos)
  nGeno <- length(mppData$geno.id)

  Q.list0 <- lapply(X = Q.pos, FUN = inc_mat_QTL, mppData = mppData,
                    Q.eff = Q.eff, order.MAF = TRUE)

  Q.list0 <- lapply(X = Q.list0, FUN =  function(x, nEnv) diag(nEnv) %x% x,
                    nEnv = nEnv)

  # expand each QTL to match the genotype information of the plot data

  ref_geno <- plot_data[, c("genotype", "env")]

  Q.list <- vector(mode = "list", length = nQTL)

  nObs <- nGeno * nEnv

  ind_row <- split(1:nObs, factor(sort(rank(1:nObs%%nEnv))))

  for(i in 1:nQTL){

    QTLdat_i <- data.frame(genotype = rep(mppData$geno.id, nEnv), Q.list0[[i]],
                           stringsAsFactors = FALSE)
    Q_i <- c()

    for(j in 1:nEnv){

      gen_j <- ref_geno[ref_geno$env == EnvNames[j], ]
      Q_data_ij <- QTLdat_i[ind_row[[j]], ]
      data_j <- merge(gen_j, Q_data_ij, by = c("genotype"))

      Q_i <- rbind(Q_i, data_j)

    }

    Q.list[[i]] <- Q_i[, -c(1, 2)]

  }

  names(Q.list) <- paste0("Q", 1:length(Q.list))
  rm(Q.list0)

  # numeric indicator to match the column of the plot data with the QTL
  # matrices (This part should be made more fluid).

  ref_geno2 <- data.frame(plot_data[, c("genotype", "env")],
                          id = 1:dim(plot_data)[1])

  ref_i <- c()

  QTLdat_i <- data.frame(genotype = rep(mppData$geno.id, nEnv),
                         stringsAsFactors = FALSE)

  for(j in 1:nEnv){

    ref_ij <- ref_geno2[ref_geno2$env == EnvNames[j], ]
    Q_data_ij <- QTLdat_i[ind_row[[j]], , drop = FALSE]
    ref_ij <- merge(ref_ij, Q_data_ij, by = c("genotype"))

    ref_i <- rbind(ref_i, ref_ij)

  }

  plot_data <- plot_data[ref_i$id, ]

  # Compute the R2

    R2.all <- R2_LR_oneS(plot_data = plot_data, mppData = mppData,
                       trait = trait, nEnv = nEnv,
                       QTL = do.call(cbind, Q.list), VCOV = VCOV,
                       exp_des_form = exp_des_form, workspace = workspace)


    R2 <- R2.all[[1]]
    R2.adj <- R2.all[[2]]


  # Partial R2
  ############

  if(glb.only) {

    QR2Res <- list(glb.R2 = R2, glb.adj.R2 = R2.adj)

  } else {

    if(nQTL > 1){

      # functions to compute the R squared or all QTL minus 1 or only 1 QTL
      # position


      part.R2.diff.LR <- function(x, QTL, plot_data, mppData, trait, nEnv, VCOV,
                                  exp_des_form, workspace) {
        R2_LR_oneS(plot_data = plot_data, mppData = mppData, trait = trait,
                 QTL = do.call(cbind, Q.list[-x]),
                 nEnv = nEnv, VCOV = VCOV, exp_des_form = exp_des_form,
                 workspace = workspace)
      }

      part.R2.sg.LR <- function(x, QTL, plot_data, mppData, trait, nEnv, VCOV,
                                exp_des_form, workspace) {
        R2_LR_oneS(plot_data = plot_data, mppData = mppData, trait = trait,
                 QTL = do.call(cbind, Q.list[x]),
                 nEnv = nEnv, VCOV = VCOV, exp_des_form = exp_des_form,
                 workspace = workspace)
      }


        R2.dif <- lapply(X = 1:nQTL, FUN = part.R2.diff.LR, QTL = Q.list,
                         plot_data = plot_data, mppData = mppData,
                         trait = trait, nEnv = nEnv, VCOV = VCOV,
                         exp_des_form = exp_des_form, workspace = workspace)

        R2_i.dif <- lapply(X = R2.dif, FUN = function(x) x[[1]])
        R2_i.dif.adj <- lapply(X = R2.dif, FUN = function(x) x[[2]])

        R2_i.dif <- R2 - unlist(R2_i.dif) # difference full model and model minus i
        R2_i.dif.adj <- R2.adj - unlist(R2_i.dif.adj)

        R2.sg <- lapply(X = 1:nQTL, FUN = part.R2.sg.LR, QTL = Q.list,
                        plot_data = plot_data, mppData = mppData, trait = trait,
                        nEnv = nEnv, VCOV = VCOV, exp_des_form = exp_des_form,
                        workspace = workspace)

        R2_i.sg <- unlist(lapply(X = R2.sg, FUN = function(x) x[[1]]))
        R2_i.sg.adj <- unlist(lapply(X = R2.sg, FUN = function(x) x[[2]]))

        names(R2_i.dif) <- names(R2_i.dif.adj) <- paste0("Q", 1:nQTL)
        names(R2_i.sg) <- names(R2_i.sg.adj) <- paste0("Q", 1:nQTL)


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
