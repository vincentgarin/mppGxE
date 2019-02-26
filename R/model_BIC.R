#############
# model_BIC #
#############

#' Compute model Bayesian Information Criterion
#'
#' Compute th Bayesian Information Criterion (BIC) of the model using :
#' BIC = -2 x log(L) + log(N) x p (Schwarz, 1978). L is the residual log
#' likelihood, N is the sample size (number of genotypes), and p is the number
#' of variance covariance parameter.
#'
#' @param mppData An object of class \code{mppData}.
#'
#' @param trait \code{Character vector} specifying which traits should be used.
#'
#' @param VCOV \code{Character} expression defining the type of variance
#' covariance structure used. "ID" for identity, "CS" for compound symmetry,
#' "DG" for heterogeneous environmental (residual) variance,
#' "UCH" for uniform covariance with heterogeneous environmental variance,
#' and "UN" for unstructured. Default = "ID".
#'
#' @author Vincent Garin
#'
#' @references
#'
#' Schwarz, G. (1978). Estimating the dimension of a model. The annals of
#' statistics, 6(2), 461-464.
#'
#' @export
#'

# library(mppR)
#
# source('F:/mppGxE/package/mppGxE/R/IncMat_cross.R')
# source('F:/mppGxE/package/mppGxE/R/IncMat_parent.R')
#
# setwd("F:/Data_mppR/EUNAM_Flint")
# load(file = './data/mpp_data/mppDataGE_toy.RData')
#
# trait <- c("PH_KWS", "PH_CIAM")
# VCOV <- "ID"

model_BIC <- function(mppData, trait, VCOV = "ID"){

  nEnv <- length(trait)
  nGeno <- length(mppData$geno.id)

  TraitEnv <- c(mppData$pheno[, trait])

  env_ind <- rep(paste0("Env", 1:nEnv), each = nGeno)
  env_ind <- factor(env_ind, levels = paste0("Env", 1:nEnv))

  cr_ind <- rep(mppData$cross.ind, times = nEnv)
  cr_ind <- factor(cr_ind, levels = unique(mppData$cross.ind))

  geno_id <- rep(mppData$geno.id, times = nEnv)
  geno_id <- factor(x = geno_id, levels = mppData$geno.id)

  # units <- factor(rep(1:nGeno, each = nEnv), levels = 1:nGeno)

  dataset <- data.frame(trait = TraitEnv, env = env_ind, genotype = geno_id,
                        cross = cr_ind)

if(VCOV == "ID"){ # Identity

  model <- asreml::asreml(fixed = trait ~ -1 + env:cross, rcov = ~ units,
                  data = dataset, trace = FALSE,
                  na.method.Y = "include", na.method.X = "omit")
  nParam <- 1


} else if (VCOV == "CS"){ # Compound symmetry

  model <- asreml::asreml(fixed = trait ~ -1 + env:cross, random = ~ genotype,
                  rcov = ~ units, data = dataset, trace = FALSE,
                  na.method.Y = "include", na.method.X = "omit")
  nParam <- 2


} else if (VCOV == "DG"){ # Diagonal - heter. env. res variance

  model <- asreml::asreml(fixed = trait ~ -1 + env:cross, rcov = ~ at(env):units,
                  data = dataset, trace = FALSE, na.method.Y = "include",
                  na.method.X = "omit")
  nParam <- nEnv


} else if (VCOV == "UCH"){

  model <- asreml::asreml(fixed = trait ~ -1 + env:cross, random = ~ genotype,
                  rcov = ~ at(env):units, data = dataset, trace = FALSE,
                  na.method.Y = "include", na.method.X = "omit")
  nParam <- nEnv + 1



} else if (VCOV == "UN"){

  model <- asreml::asreml(fixed = trait ~ -1 + env:cross, rcov = ~ us(env):genotype,
                  data = dataset, trace = FALSE, na.method.Y = "include",
                  na.method.X = "omit")
  nParam <- (nEnv * (nEnv + 1))/2

}

  BIC = (-2*(summary(model)$loglik)) + (log(nGeno) * nParam)

  return(BIC)

}
