####################
# mppGE_oneS_QTL_effects #
####################

#' MPP GxE one stage QTL genetic effects
#'
#' Compute MPP GxE one stage QTL genetic effects.
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
#' @param workspace Size of workspace for the REML routines measured in double
#' precision words (groups of 8 bytes). The default is workspace = 8e6.
#'
#' @return Return:
#'
#' \item{Qeff}{\code{List} of \code{data.frame} (one per QTL) containing the
#' following information:
#'
#' \enumerate{
#'
#' \item{QTL genetic effects}
#' \item{Standard error of the QTL effects.}
#' \item{Test statistics of the effects (Wald statistic).}
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
#' library(asreml)
#'
#' data(mppData_GE)
#' data(plot_data)
#'
#' Qpos <- c("PZE.105068880", "PZE.106098900")
#'
#' Qeff <- mppGE_oneS_QTL_effects(plot_data = plot_data, mppData = mppData_GE,
#'                               trait = 'DMY', Q.eff = 'par',
#'                               exp_des_form = 'env:Rep + env:Rep:Block',
#'                               QTL = Qpos)
#'
#' Qeff
#'
#' @export
#'


mppGE_oneS_QTL_effects <- function(plot_data, mppData, trait, Q.eff = "cr",
                             VCOV = "CS_CSRT", exp_des_form, QTL = NULL,
                             workspace = 8e6){

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

  Q.names <- function(x, Q.list, nEnv){
    rep(paste0("Q", x, attr(Q.list[[x]], "dimnames")[[2]]), nEnv)
  }

  names.QTL <- unlist(lapply(X = 1:nQTL, FUN = Q.names, Q.list = Q.list0,
                             nEnv = nEnv))

  if(Q.eff == "anc"){

    n_al <- unlist(lapply(X = Q.list0, FUN = function(x) dim(x)[2]))

    e_lab <- paste0("E", 1:nEnv)

    Env.names <- lapply(X = n_al, FUN = function(x, e_lab) rep(e_lab, each = x),
                        e_lab = e_lab)

    Env.names <- unlist(Env.names)

  } else {

    n_al <- NULL

    Env.names <- rep(rep(paste0("E", 1:nEnv), each = dim(Q.list0[[1]])[2]), nQTL)

  }

  names.QTL <- paste(names.QTL, Env.names, sep = "_")



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

  # model computation

  model <- QTLModelQeff_oneS(plot_data = plot_data, mppData = mppData,
                             trait = trait, Q.list = Q.list,
                             VCOV = VCOV, exp_des_form = exp_des_form,
                             names.QTL = names.QTL, workspace = workspace)


  # process the results

  Qeff <- Qeff_res_processing_GE(model = model, mppData = mppData,
                                 Q.eff = Q.eff, VCOV = VCOV,
                                 names.QTL = names.QTL, Q.pos = Q.pos,
                                 nQTL = nQTL, n_al = n_al, nEnv = nEnv)


  return(Qeff)


}
