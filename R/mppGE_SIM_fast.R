##################
# mppGE_SIM_fast #
##################

#' MPP GxE Simple Interval Mapping using fast MM GLS algorithm
#'
#' Computes single QTL models along the genome using different models. Add extra
#' elements to provide more information about the procedure. The initial mixed
#' are calculated with functions from the nlme package.
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
#' covariance structure used. 'CS' for compound symmetry assuming a unique
#' genetic covariance between environments. 'CSE' for cross-specific within
#' environment error term. 'CS_CSE' for both compound symmetry plus
#' cross-specific within environment error term. 'UN' for unstructured
#' environmental variance covariance structure allowing a specific genotypic
#' covariance for each pair of environments. Default = 'UN'
#'
#' @param n.cores \code{Numeric}. Specify here the number of cores you like to
#' use. Default = 1.
#' 
#' @param maxIter maximum number of iterations for the lme optimization algorithm.
#' Default = 100.
#' 
#' @param msMaxIter maximum number of iterations for the optimization step inside
#' the lme optimization. Default = 100.
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
#' @references
#' 
#' Pinheiro J, Bates D, DebRoy S, Sarkar D, R Core Team (2021). nlme: Linear
#' and Nonlinear Mixed Effects Models_. R package version 3.1-152,
#' <URL: https://CRAN.R-project.org/package=nlme>.
#'
#' @seealso \code{\link{plot_genEffects_GE}}
#'
#' @examples
#'
#' library(mppR)
#'
#' data(mppData_GE)
#'
#' SIM <- mppGE_SIM_fast(mppData = mppData_GE, trait = c('DMY_CIAM', 'DMY_TUM'),
#'                  Q.eff = 'par', plot.gen.eff = TRUE)
#'
#' Qpos <- QTL_select(Qprof = SIM, threshold = 3, window = 50)
#'
#' plot(x = SIM, QTL = Qpos)
#'
#' plot_genEffects_GE(mppData = mppData_GE, nEnv = 2, EnvNames = c('CIAM', 'TUM'),
#'                    Qprof = SIM, Q.eff = 'par', QTL = Qpos, text.size = 14)
#'
#' @export
#'

# library(mppR)
# library(mppGxE)
# library(nlme)
# library(scrime)

# source('D:/hard_drive_back_up/mppGxE/package/mppGxE/R/check_mod_mppGE.R')
# source('D:/hard_drive_back_up/mppGxE/package/mppGxE/R/MM_comp.R')
# source('D:/hard_drive_back_up/mppGxE/package/mppGxE/R/getVCOV.R')
# source('D:/hard_drive_back_up/mppGxE/package/mppGxE/R/pair_index.R')
# source('D:/hard_drive_back_up/mppGxE/package/mppGxE/R/W_QTL.R')
# source('D:/hard_drive_back_up/mppGxE/package/mppGxE/R/check.inf.R')
# # source('D:/hard_drive_back_up/mppGxE/package/mppGxE/R/W_stat.R')
# 
# data(mppData_GE)
# mppData = mppData_GE
# trait = c('DMY_CIAM', 'DMY_TUM')
# Q.eff = 'par'
# VCOV = "CS"
# n.cores = 1
# maxIter = 100
# msMaxIter = 100

# # alternative data (Grinkan BC-NAM)
# 
# load('D:/hard_drive_back_up/ICRISAT/BCNAM/data/mppData/Grinkan_mppData.RData')
# # load('D:/hard_drive_back_up/ICRISAT/BCNAM/data/genotype/Grinkan_geno012.RData')
# # mppData$geno.IBS <- geno012
# 
# trait <- c("GR1_GYLD", 'GR2_GYLD', 'GR3_GYLD', 'GR4_GYLD')
# VCOV <- 'UN'
# Q.eff <- 'par'
# n.cores = 1
# maxIter = 100
# msMaxIter = 100

mppGE_SIM_fast <- function(mppData, trait, Q.eff = "cr", VCOV = "UN",
                           n.cores = 1, maxIter = 100,
                           msMaxIter = 100) {
  
  ### 1. Check data format and arguments
  check_mod_mppGE(mppData = mppData, trait = trait, Q.eff = Q.eff, VCOV = VCOV,
                  QTL_ch = FALSE, fast = TRUE)
  
  
  ### 2. Form required elements for the analysis
  nEnv <- length(trait)
  TraitEnv <- c(mppData$pheno[, trait])
  NA_id <- is.na(TraitEnv)
  vect.pos <- 1:dim(mppData$map)[1]
  
  ### 3. Compute the initial VCOV structure
  m <- MM_comp(mppData = mppData, nEnv = nEnv, y = TraitEnv, VCOV = VCOV,
               maxIter = maxIter, msMaxIter = msMaxIter)
  
  Vi <- getVCOV(mppData = mppData, model = m$model, VCOV = VCOV, data = m$data,
                nEnv = nEnv, inv = TRUE)
  
  ### 4. Computation of the Wald statistic
  
  cross_mat <- model.matrix(~ cross_env -1, data = m$data)
  cross_mat <- cross_mat[!NA_id, ]
  
  
  # 4.1 Optional cluster construction for parallele computation
  if(n.cores > 1){ parallel <- TRUE; cluster <- makeCluster(n.cores)
    
  } else { parallel <- FALSE; cluster <- NULL }
  
  # 4.2 Iteration along the genome
  if (parallel) {
    
    log.pval <- parLapply(cl = cluster, X = vect.pos, fun = W_QTL,
                          y = TraitEnv[!NA_id], Vi = Vi,
                          mppData = mppData,  nEnv = nEnv, 
                          Q.eff = Q.eff, cross_mat = cross_mat, NA_id = NA_id)
    
  } else {
    
    log.pval <- lapply(X = vect.pos, FUN = W_QTL,
                       y = TraitEnv[!NA_id], Vi = Vi,
                       mppData = mppData,  nEnv = nEnv, 
                       Q.eff = Q.eff, cross_mat = cross_mat, NA_id = NA_id)
    
  }
  
  if(n.cores > 1){stopCluster(cluster)}
  
  log.pval <- t(data.frame(log.pval))
  log.pval[, 1] <- check.inf(x = log.pval[, 1]) # check if there are -/+ Inf value
  log.pval[is.na(log.pval[, 1]), 1] <- 0
  
  ### 5. form the results
  
  SIM <- data.frame(mppData$map, log.pval)
  
    if(Q.eff == "cr"){
      
      Env_name <- rep(paste0("_E", 1:nEnv), each = mppData$n.cr)
      Qeff_names <- paste0(rep(unique(mppData$cross.ind), nEnv), Env_name)
      
      colnames(SIM)[5:dim(SIM)[2]] <- c("log10pval", Qeff_names)
      
    } else if (Q.eff == "anc") {
      
      Env_name <- rep(paste0("_E", 1:nEnv), each = mppData$n.par)
      Qeff_names <- paste0(rep(mppData$parents, nEnv), Env_name)
      
      colnames(SIM)[5:dim(SIM)[2]] <- c("log10pval", Qeff_names)
      
    } else { # par and bi-allelic no modif for gen effects names
      
      colnames(SIM)[5] <- "log10pval"
      
    }
  
  class(SIM) <- c("QTLprof", "data.frame")
  
  # 5.1: Verify the positions for which model could not be computed
  
  if(sum(SIM$log10pval == 0) > 0) {
    
    if (sum(SIM$log10pval) == 0){
      
      message("The computation of the QTL models failled for all positions.")
      
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