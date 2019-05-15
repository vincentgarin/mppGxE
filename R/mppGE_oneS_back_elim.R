##################
# mppGE_oneS_back_elim #
##################

#' Backward elimination of one stage MPP GxE analysis
#'
#' Backward elimination on a list of QTL using a one stage MPP GxE
#' model.
#'
#' @param plot_data \code{Data.frame} containing the plot data with the following
#' columns: the trait(s), 'genotype' (genotype indicator), 'check'
#' (check indicator), 'cross' (cross indicator), 'env' (environment indicator),
#' and all other experimental design covariates (e.g. replicate, blocks, etc.).
#' The column names must be ('genotype', 'check', 'cross', env'). The names of
#' the experimental design covariates must be the same as the one used in
#' 'exp_des_form'. For more details see \code{\link{plot_data}}.
#'
#' @param mppData Object of class \code{mppData} contaning the the same
#' genotype identifiers as the one in \code{plot_data} ('genotype').
#'
#' @param trait \code{Character} expression for the trait matching the trait
#' column in 'plot_data' argument.
#'
#' @param Q.eff \code{Character} expression indicating the assumption concerning
#' the QTL effects: 1) "cr" for cross-specific; 2) "par" for parental; 3) "anc"
#' for ancestral; 4) "biall" for a bi-allelic. Default = "cr".
#'
#' @param VCOV VCOV \code{Character} expression defining the type of variance
#' covariance structure used: a) "CSRT" for within environment
#' cross-specific residual terms; b) "CS_CSRT" for compound symmetry with within
#' environment cross-specific residual terms. Default = "CS_CSRT".
#'
#' @param exp_des_form \code{Character} expression for the random experimental
#' design effects in asreml-R format. For example,
#' 'env:replicate + env:replicate:block'. The variable names used in
#' 'exp_des_form' should strictly match the column names used in 'plot_data'.
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
#' data(plot_data)
#'
#' Qpos <- c("PZE.105068880", "PZE.106098900")
#'
#' QTL <- mppGE_oneS_back_elim(plot_data = plot_data, mppData = mppData_GE,
#'                             trait = 'DMY', Q.eff = 'par',
#'                             exp_des_form = 'env:Rep + env:Rep:Block',
#'                             QTL = Qpos)
#'
#' @export
#'


mppGE_oneS_back_elim <- function(plot_data, mppData, trait, Q.eff = "cr",
                           VCOV = "CS_CSRT", exp_des_form, QTL = NULL,
                           alpha = 0.01, workspace = 8e6){

  # Check data format and arguments

  check_mod_mppGE(mppData = mppData, trait = trait, Q.eff = Q.eff, VCOV = VCOV,
                  QTL_ch = TRUE, QTL = QTL, GE = FALSE, plot_data = plot_data,
                  exp_des_form = exp_des_form)

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
      data_j <- merge(gen_j, Q_data_ij, by = c("genotype"), all.x = TRUE)

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
    ref_ij <- merge(ref_ij, Q_data_ij, by = c("genotype"), all.x = TRUE)

    ref_i <- rbind(ref_i, ref_ij)

  }

  plot_data <- plot_data[ref_i$id, ]


  # Compute the model

  ind <- TRUE

  while(ind) {

    ### 3.1 elaboration of model formulas

    model.formulas <- formula_backward_GE(Q.names = names(Q.list), VCOV = VCOV,
                                          type = 'oneS')

    ### 3.2 computation of the models

    pvals <- lapply(X = model.formulas, FUN = QTLModelBack_oneS,
                    plot_data = plot_data, mppData = mppData,
                    trait = trait, nEnv = nEnv, Q.list = Q.list, VCOV = VCOV,
                    exp_des_form = exp_des_form, workspace = workspace)

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
