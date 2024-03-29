###################
# mppGE_proc_fast #
###################

#' MPP GxE QTL analysis
#'
#' QTL detection in MPP characterized in multiple environments.
#'
#' The procedure is the following:
#'
#' \enumerate{
#'
#' \item{Simple interval mapping (SIM) to select cofactors
#' (\code{\link{mppGE_SIM_fast}}).}
#'
#' \item{Composite interval mapping (CIM) with selected cofactors
#' (\code{\link{mppGE_CIM_fast}}).}
#'
#' \item{Estimation of QTLs additive allelic effect
#' (\code{\link{QTL_effect_GE}}).}
#' 
#' \item{Estimation of the global QTLs effect R squared and individual QTL effect
#' R squared (\code{\link{QTL_R2_GE}}).}
#'
#' }
#'
#'
#' @param pop.name \code{Character} name of the studied population.
#' Default = "MPP".
#'
#' @param trait.name \code{Character} name of the studied trait.
#' Default = "trait1".
#'
#' @param mppData An object of class \code{mppData}.
#'
#' @param trait \code{Character vector} specifying which traits (environments) should be used.
#'
#' @param EnvNames \code{character} expression indicating the environment or trait
#' labels. By default: Env_1, 2, 3, etc.
#'
#' @param Q.eff \code{Character} expression indicating the assumption concerning
#' the QTL effect: 1) "cr" for cross-specific effects; 2) "par" parental
#' effects; 3) "anc" for an ancestral effects; 4) "biall" for a bi-allelic
#' effects. Default = "par".
#'
#' @param VCOV VCOV \code{Character} expression defining the type of variance
#' covariance structure used. 'CS' for compound symmetry assuming a unique
#' genetic covariance between environments. 'CSE' for cross-specific within
#' environment error term. 'CS_CSE' for both compound symmetry plus
#' cross-specific within environment error term. 'UN' for unstructured
#' environmental variance covariance structure allowing a specific genotypic
#' covariance for each pair of environments. Default = 'UN'
#' 
#' @param VCOV_data \code{Character} specifying if the reference VCOV of the
#' CIM profile computation should be formed  taking all cofactors into
#' consideration ("unique") or if different VCOVs should be formed by removing
#' the cofactor information that is too close of a tested cofactor position
#' ("minus_cof"). Default = "unique"
#'
#' @param thre.cof \code{Numeric} value representing the -log10(p-value)
#' threshold above which a position can be selected as cofactor. Default = 4.
#'
#' @param win.cof \code{Numeric} value in centi-Morgan representing the minimum
#' distance between two selected cofactors. Default = 50.
#'
#' @param window \code{Numeric} distance (cM) on the left and the right of a
#' cofactor position where it is not included in the model. Default = 20.
#'
#' @param thre.QTL \code{Numeric} value representing the -log10(p-value)
#' threshold above which a position can be selected as QTL. Default = 4.
#'
#' @param win.QTL \code{Numeric} value in centi-Morgan representing the minimum
#' distance between two selected QTLs. Default = 20.
#'
#' @param text.size \code{Numeric} value specifying the size of graph axis text
#' elements. Default = 18.
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
#' @param verbose \code{Logical} value indicating if the steps of mpp_proc should
#' be printed. Default = TRUE.
#'
#' @param output.loc Path where a folder will be created to save the results.
#' Default = NULL.
#'
#' @return Return:
#'
#' List containing the following items:
#'
#' \item{n.QTL}{Number of detected QTLs.}
#'
#' \item{cofactors}{\code{Data.frame} with cofactors positions.}
#'
#' \item{QTL}{\code{Data.frame} with QTL positions.}
#' 
#' \item{Q_eff}{\code{list} containing the estimated QTL allelic effects.}
#'
#' \item{R2}{\code{List} containing R squared statistics of the QTL effects}
#'
#' Some output files are also saved at the specified location
#' (\code{output.loc}):
#'
#' \enumerate{
#'
#' \item{The SIM and CIM results in a RData file (SIM.RData, CIM.RData).}
#'
#' \item{The list of cofactors (cofactors.RData).}
#'
#' \item{The list of QTL (QTLs.RData).}
#' 
#' \item{The list of QTL allelic effects (QTL_effects.RData).}
#'
#' \item{The QTL R squared statistics (QTL_R2.RData)}
#' 
#' \item{The number of detected QTLs and adjusted R2 (Glb_res.RData)}
#'
#' \item{The plot of the CIM profile (QTL_profile.pdf) with dotted vertical
#' lines representing the cofactors positions and the
#' plot of the genetic effects per cross or parents obtained with
#' \code{\link{plot_allele_eff_GE}} (gen_eff.pdf) with dashed
#' lines representing the QTL positions.}
#'
#' }
#'
#'
#' @author Vincent Garin
#'
#' @seealso
#'
#' \code{\link{mppGE_CIM_fast}},
#' \code{\link{mppGE_SIM_fast}},
#' \code{\link{QTL_effect_GE}},
#' \code{\link{QTL_R2_GE}}
#'
#' @examples
#'
#' \dontrun{
#'
#' data(mppData_GE)
#'
#' MPP_GE_QTL <- mppGE_proc_fast(pop.name = 'EUNAM', trait.name = 'DMY',
#' mppData = mppData_GE, trait = c('DMY_CIAM', 'DMY_TUM'),
#' n.cores = 1, output.loc = tempdir())
#'
#' }
#'
#' @export
#'


mppGE_proc_fast <- function(pop.name = "MPP", trait.name = "trait1", mppData, trait,
                            EnvNames = NULL,  Q.eff = "par", VCOV = "UN",
                            VCOV_data = "unique", thre.cof = 4, win.cof = 50, window = 20,
                            thre.QTL = 4, win.QTL = 20, text.size = 18,
                            n.cores = 1,  maxIter = 100, msMaxIter = 100,
                            verbose = TRUE, output.loc = NULL) {
  
  
  ##### Check the  parameters #####
  
  # CheckModelGE(mppData = mppData, trait = trait, Q.eff = Q.eff, VCOV = VCOV,
  #              EnvNames = EnvNames, fct = "SIM")
  
  ##### Extra parameter + directory #####
  
  if(is.null(EnvNames)){
    
    EnvNames <- paste0("Env_", 1:length(trait))
    
  }
  
  folder.loc <- file.path(output.loc, paste("QTLGE", pop.name, trait.name, Q.eff,
                                            VCOV, sep = "_"))
  
  dir.create(folder.loc)
  
  
  ##### Cofactors selection - SIM ####
  
  if(verbose){
    
    cat("\n")
    cat("Cofactors selection - SIM")
    cat("\n")
    cat("\n")
    
  }
  
  SIM <- mppGE_SIM_fast(mppData = mppData, trait = trait, Q.eff = Q.eff,
                        VCOV = VCOV, n.cores = n.cores, maxIter = maxIter,
                        msMaxIter = msMaxIter)
  
  # save SIM results in output location
  save(SIM, file = file.path(folder.loc, "SIM.RData"))
  
  # cofactors selection
  
  cofactors <- QTL_select(Qprof = SIM, threshold = thre.cof, window = win.cof)
  
  
  if (is.null(cofactors)) { # test if cofactors have been selected
    
    message("No QTL/cofactor position detected based on the SIM profile.")
    
    return(NULL)
    
  }
  
  ##### Multi-QTL model search - CIM #####
  
  if(verbose){
    
    cat("\n")
    cat("Multi-QTL model search - CIM")
    cat("\n")
    cat("\n")
    
  }
  
  CIM <- mppGE_CIM_fast(mppData = mppData, trait = trait, Q.eff = Q.eff,
                        VCOV = VCOV, VCOV_data = VCOV_data, cofactors = cofactors,
                        window = window, n.cores = n.cores, maxIter = maxIter,
                        msMaxIter = msMaxIter)
  
  # save the list of cofactors
  save(cofactors, file = file.path(folder.loc, "cofactors.RData"))
  
  # save CIM results
  save(CIM, file = file.path(folder.loc, "CIM.RData"))
  
  # select QTL candidates
  QTL <- QTL_select(Qprof = CIM, threshold = thre.QTL, window = win.QTL)
  
  if (is.null(QTL)) { # test if QTL have been selected
    
    message("No QTL position detected based on the (last) CIM profile.")
    # return(NULL)
    QTL <- cofactors
    CIM_fail <- TRUE
    
  } else {
    
    CIM_fail <- FALSE
    
    
  }
  
  # save CIM results
  save(QTL, file = file.path(folder.loc, "QTLs.RData"))
  
  
  ##### QTL effects #####
  Q_eff <- QTL_effect_GE(mppData = mppData, trait = trait, Q.eff = Q.eff,
                          QTL = QTL, VCOV = VCOV, maxIter = maxIter,
                          msMaxIter = msMaxIter)
  
  save(Q_eff, file = file.path(folder.loc, "QTL_effects.RData"))
  
  ##### QTL R2 #####
  R2 <- QTL_R2_GE(mppData = mppData, trait = trait, Q.eff = Q.eff, QTL = QTL)
  
  # save R2 results
  QTL.R2 <- data.frame(QTL[, 1:5], round(R2[[3]], 2), round(R2[[4]], 2),
                       stringsAsFactors = FALSE)
  
  colnames(QTL.R2)[6:7] <- c("R2.diff", "adj.R2.diff")
  
  save(QTL.R2, file = file.path(folder.loc, "QTL_R2.RData"))
  
  # save N QTL and adjusted R2
  glb_res <- list(N_QTL = nrow(QTL), adj_R2 = R2[[2]])
  save(glb_res, file = file.path(folder.loc, "Glb_res.RData"))
  
  ##### plot and results processing #####
  
  if(verbose){
    
    cat("\n")
    cat("Results processing")
    cat("\n")
    cat("\n")
    
  }
  
  if(CIM_fail){
    
    Qprof <- SIM
    main_prof <- paste("SIM", pop.name, trait.name, Q.eff, VCOV)
    
  } else {
    
    Qprof <- CIM
    main_prof <- paste("CIM", pop.name, trait.name, Q.eff, VCOV)
    
  }
  
  
  main.Qeff <- paste("QTL gen. effects", pop.name, trait.name, Q.eff, VCOV)
  
  if(Q.eff == "biall"){t_plot <- "h"} else {t_plot <- "l"}
  
  pdf(file.path(folder.loc, "QTL_profile.pdf"), height = 10, width = 16)
  
  print(plot(x = Qprof, QTL = cofactors, type = t_plot, main = main_prof,
             threshold = thre.QTL, text.size = text.size))
  
  dev.off()
  
  pdf(file.path(folder.loc, "gen_eff.pdf"), height = 10, width = 16)
  
  print(plot_allele_eff_GE(mppData = mppData, nEnv = length(trait),
                           EnvNames = EnvNames, Qprof = Qprof, Q.eff = Q.eff,
                           QTL = QTL, main = main.Qeff, text.size = text.size))
  
  dev.off()
  
  #### QTL report ####
  
  QTL_report_GE(out.file = file.path(folder.loc, "QTL_REPORT.txt"),
                main = paste(pop.name, trait.name, Q.eff, VCOV),
                QTL.info = QTL[, c(1, 2, 4, 5)], QTL.effects = Q_eff,
                R2 = R2)
  
  ### Return R object
  
  results <- list(n.QTL = dim(QTL)[1], cofactors = cofactors[, 1:5],
                  QTL = QTL[, 1:5], Q_eff = Q_eff, R2 = R2)
  
  return(results)
  
}