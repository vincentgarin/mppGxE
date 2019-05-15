#################
# mppGE_oneS_SIM #
#################

#' MPP GxE one stage SIM analysis
#'
#' Perform a MPP GxE one stage SIM analysis.
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
#' @param plot.gen.eff \code{Logical} value. If \code{plot.gen.eff = TRUE},
#' the function will save the significance of the QTL allelic effects per
#' cross/parent along the genome. These results can be visualized with the
#' function \code{\link{plot_genEffects_GE}}. Default value = FALSE.
#'
#' @param workspace Size of workspace for the REML routines measured in double
#' precision words (groups of 8 bytes). The default is workspace = 8e6.
#'
#' @return Return:
#'
#' \item{SIM }{\code{Data.frame} of class \code{QTLprof}. with five columns :
#' 1) QTL marker or in between position names; 2) chromosomes;
#' 3) interger position indicators on the chromosome;
#' 4) positions in centi-Morgan; and 5) -log10(p-val). And if
#' \code{plot.gen.eff = TRUE}, p-values of the cross or parental QTL allelic effects.}
#'
#' @author Vincent Garin
#'
#' @seealso \code{\link{plot_genEffects_GE}}
#'
#' @examples
#'
#'\dontrun{
#'
#' library(asreml)
#' library(mppR)
#'
#' data(mppData_GE)
#' data(plot_data)
#'
#' # take around 5 min
#'
#' SIM <- mppGE_oneS_SIM(plot_data = plot_data, mppData = mppData_GE,
#'                       trait = 'DMY', Q.eff = 'par',
#'                       exp_des_form = 'env:Rep + env:Rep:Block',
#'                       plot.gen.eff = TRUE)
#'
#' Qpos <- QTL_select(Qprof = SIM, threshold = 3, window = 50)
#'
#' plot(x = SIM, QTL = Qpos)
#'
#' plot_genEffects_GE(mppData = mppData_GE, nEnv = 2, EnvNames = c('CIAM', 'TUM'),
#'                    Qprof = SIM, Q.eff = 'par', QTL = Qpos, text.size = 14)
#'
#'}
#'
#' @export
#'



mppGE_oneS_SIM <- function(plot_data, mppData, trait, Q.eff = "cr", VCOV = "CS_CSRT",
                          exp_des_form, plot.gen.eff = FALSE, workspace = 8e6){

  # Check data format and arguments

  check_mod_mppGE(mppData = mppData, trait = trait, Q.eff = Q.eff, VCOV = VCOV,
                  QTL_ch = FALSE, GE = FALSE, plot_data = plot_data,
                  exp_des_form = exp_des_form)

  # 2. Determine the environments

  EnvNames <- unique(plot_data$env)

  nEnv <- length(EnvNames)

  vect.pos <- 1:dim(mppData$map)[1]

  log.pval <- lapply(X = vect.pos, FUN = QTLModelSIM_oneS, plot_data = plot_data,
                     mppData = mppData, trait = trait, nEnv = nEnv,
                     EnvNames = EnvNames, Q.eff = Q.eff, VCOV = VCOV,
                     exp_des_form = exp_des_form, plot.gen.eff = plot.gen.eff,
                     workspace = workspace)


  log.pval <- t(data.frame(log.pval))
  # if(plot.gen.eff & (VCOV == "ID")){log.pval[is.na(log.pval)] <- 1}
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
