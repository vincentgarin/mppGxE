##################
# QTL_effects_GE #
##################

#' MPP GxE QTL genetic effects
#'
#' Compute MPP GxE QTL genetic effects.
#'
#' @param mppData An object of class \code{mppData}.
#'
#' @param trait \code{Character vector} specifying which traits (environments) should be used.
#'
#' @param Q.eff \code{Character} expression indicating the assumption concerning
#' the QTL effects: 1) "cr" for cross-specific; 2) "par" for parental; 3) "anc"
#' for ancestral. Default = "par".
#'
#' @param VCOV VCOV \code{Character} expression defining the type of variance
#' covariance structure used. 'CS' for compound symmetry assuming a unique
#' genetic covariance between environments. 'CSE' for cross-specific within
#' environment error term. 'CS_CSE' for both compound symmetry plus
#' cross-specific within environment error term. 'UN' for unstructured
#' environmental variance covariance structure allowing a specific genotypic
#' covariance for each pair of environments. Default = 'UN'
#'
#' @param QTL Object of class \code{QTLlist} representing a list of
#' selected marker positions obtained with the function QTL_select() or
#' a vector of \code{character} marker positions names. Default = NULL.
#'
#' @param maxIter maximum number of iterations for the lme optimization algorithm.
#' Default = 100.
#' 
#' @param msMaxIter maximum number of iterations for the optimization step inside
#' the lme optimization. Default = 100.
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
#' \item{Wald statistics of the effects.}
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
#' data(mppData_GE)
#'
#' Qpos <- c("PZE.105068880", "PZE.106098900")
#'
#' Qeff <- QTL_effects_GE(mppData = mppData_GE, trait = c('DMY_CIAM', 'DMY_TUM'),
#'                           Q.eff = 'par', QTL = Qpos)
#'
#' Qeff
#'
#' @export
#'

# library(mppR)
# library(mppGxE)
# library(nlme)
# library(Matrix)
# 
# source("~/WD/mppGxE/package/mppGxE/R/check_mod_mppGE.R")
# source("~/WD/mppGxE/package/mppGxE/R/getVCOV.R")
# source("~/WD/mppGxE/package/mppGxE/R/MM_comp.R")
# source("~/WD/mppGxE/package/mppGxE/R/W_QTL.R")
# source("~/WD/mppGxE/package/mppGxE/R/check.inf.R")
# source("~/WD/mppGxE/package/mppGxE/R/lme_comp.R")
# source("~/WD/mppGxE/package/mppGxE/R/sign.star.R")
# 
# setwd('C:/Users/user/Documents/WD/ICRISAT/BCNAM/data')
# load('./mppData/KK2013/mppData.RData')
# 
# map <- mppData$map
# 
# set.seed(435436)
# 
# mk_sel <- map$mk.names[map$chr %in% c(1, 3)]
# mk_sel <- mk_sel[sort(sample(x = 1:length(mk_sel), 500))]
# 
# mppData <- subset(x = mppData, mk.list = mk_sel)
# 
# #### SIM computation and cofactors selection ####
# 
# tr_i <- c("KK1_G_YIELD", "KK2_G_YIELD", "KK3_G_YIELD", "KK4_G_YIELD")
# 
# # SIM <- mppGE_SIM_fast(mppData = mppData, trait = tr_i,
# #                       Q.eff = 'par', n.cores = 3)
# # 
# # QTL <- QTL_select(Qprof = SIM, threshold = 3, window = 200)
# 
# QTL <- c('S1_61851589', 'S3_53134126')
# 
# trait = tr_i
# Q.eff = 'par'
# VCOV = 'UN'
# maxIter = 100
# msMaxIter = 100


QTL_effects_GE <- function(mppData, trait, Q.eff = "par", VCOV = "UN",
                           QTL = NULL, maxIter = 100, msMaxIter = 100){
  
  #### 1. Check data format and arguments ####
  check_mod_mppGE(mppData = mppData, trait = trait, Q.eff = Q.eff, VCOV = VCOV,
                  QTL_ch = TRUE, QTL = QTL, fast = TRUE, CIM = FALSE)
  
  #### 2. Form required elements for the analysis ####
  nEnv <- length(trait)
  TraitEnv <- c(mppData$pheno[, trait])
  NA_id <- is.na(TraitEnv)
  
  #### 3. QTL matrices ####
  if(is.character(QTL)){
    
    QTL.pos <- which(mppData$map$mk.names %in% QTL)
    
  } else {
    
    QTL.pos <- which(mppData$map$mk.names %in% QTL$mk.names)
    
  }
  
  nQTL <- length(QTL.pos)
  
  QTL_list <- mapply(FUN = inc_mat_QTL, x = QTL.pos,
                     MoreArgs = list(Q.eff = Q.eff, mppData = mppData),
                     SIMPLIFY = FALSE)
  
  QTL_list <- lapply(QTL_list, function(x) x[, -1])
  
  nAllele <- sapply(QTL_list, function(x) ncol(x))
  
  
  # combined QTL matrices
  QTL_mat <- do.call(cbind, QTL_list)
  
  #### 4. Computation of the mixed model ####
  
  nGeno <- dim(mppData$pheno)[1]
  env <- rep(paste0('E', 1:nEnv), each = nGeno)
  cross <- rep(mppData$cross.ind, nEnv)
  geno <- rep(rownames(mppData$pheno), nEnv)
  cross_env <- paste0(cross, '_', env)
  
  Q_nm <- colnames(QTL_mat)
  Q_id <- paste0('QTL_', rep(1:nQTL, nAllele))
  QTL_nm <- paste0(Q_id, '_', Q_nm)
  QTL_nm <- paste0(rep(QTL_nm, nEnv), rep(paste0('_E', 1:nEnv), each = length(Q_nm)))
  
  # remove problematic character in parents names
  QTL_nm <- gsub(pattern = '-', replacement = "", x = QTL_nm)
  QTL_nm <- gsub(pattern = ' ', replacement = "", x = QTL_nm)
  
  QTL_mat <- diag(nEnv) %x% QTL_mat
  colnames(QTL_mat) <- QTL_nm
  
  d <- data.frame(trait = TraitEnv, env = env, cross_env = cross_env, geno = geno)
  d[, 2:4] <- lapply(d[, 2:4], as.factor)
  d <- data.frame(d, QTL_mat)
  
  fix_form <- paste0('trait~-1 + cross_env+', paste(QTL_nm, collapse = '+'))
  
  m_sg <- lm(as.formula(fix_form), data = d)
  coeff <- coefficients(m_sg)
  if(any(is.na(coeff))){
    d <- d[, -which(colnames(d) %in% names(coeff[is.na(coeff)]))]
  }
  
  Q_id <- colnames(d)[5:ncol(d)]
  
  fix_form <- paste0('trait~-1 + cross_env+', paste(Q_id, collapse = '+'))
  
  m <- lme_comp(fix_form = fix_form, VCOV = VCOV, data = d,
                maxIter = maxIter, msMaxIter = msMaxIter)
  
  Vi <- getVCOV(mppData = mppData, model = m, VCOV = VCOV,
                data = d, nEnv = nEnv, inv = TRUE)
  
  ##### Get QTL effect (Beta), standard error, Wald test #####
  
  Beta <- m$coefficients$fixed
  Q_ind <- grepl(pattern = 'QTL_', x = names(Beta))
  B_QTL <- Beta[Q_ind]
  B_SD <- sqrt(diag(m$varFix)[Q_ind])
  
  # Wald test
  V_QTL_inv <- qr.solve(m$varFix[Q_ind, Q_ind])
  W_Qa <- rep(NA, length(B_QTL))
  for(i in 1:length(W_Qa)) W_Qa[i] <- (B_QTL[i]^2) * diag(V_QTL_inv)[i]
  W_pval <- pchisq(W_Qa, 1, lower.tail = FALSE)
  W_sign <- sapply(W_pval, sign.star) 
  
  res_tab <- data.frame(Effect = round(B_QTL, 3), Std.dev = round(B_SD, 3),
                        Wald = round(W_Qa, 2), df = 1, p.val = W_pval,
                        sign = W_sign)
  
  
  #### process the results ####
  
  if(Q.eff == 'par'){
    
    p_nm <- mppData$parents
    p_nm <- gsub(pattern = '-', replacement = "", x = p_nm)
    p_nm <- gsub(pattern = ' ', replacement = "", x = p_nm)
    
    ref_QTL <- rep(paste0('QTL_', 1:nQTL), each = (mppData$n.par * nEnv))
    ref_QTL <- paste0(ref_QTL, '_', rep(p_nm, nEnv))
    ref_QTL <- paste0(ref_QTL, '_', rep(paste0('E', 1:nEnv), each = mppData$n.par))
    ref_QTL <- data.frame(ref_QTL)
    
    res_tab <- data.frame(ref_QTL = rownames(res_tab), res_tab)
    
    res_tab <- merge(x = ref_QTL, y = res_tab, by = 'ref_QTL', all.x = TRUE)
    rownames(res_tab) <- res_tab$ref_QTL
    res_tab <- res_tab[ref_QTL$ref_QTL, ]
    res_tab <- res_tab[, -1]
    
    Q_f <- strsplit(x = rownames(res_tab), split = '_')
    Q_f_1 <- unlist(lapply(Q_f, `[[`, 1))
    Q_f_2 <- unlist(lapply(Q_f, `[[`, 2))
    Q_f <- paste0(Q_f_1, '_', Q_f_2)
    Q_f <- factor(Q_f, levels = paste0('QTL_', unique(Q_f_2)))
    
    Qeff.mat <- split(x = res_tab, f = Q_f)
    
    
  } else if (Q.eff == 'anc'){
    
    # Later
    
  } else {
    
    # later
    
  }
  
  return(Qeff.mat)
  
}