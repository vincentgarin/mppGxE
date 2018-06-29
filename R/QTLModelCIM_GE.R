##################
# QTLModelCIM_GE #
##################

# function to compute a single position MPP GxE QTL CIM model

QTLModelCIM_GE <- function(x, mppData, nEnv, TraitEnv, cross.mat, CrMatEnv,
                           par.mat, Q.eff, par.clu, VCOV, cof.list, cof.part,
                           plot.gen.eff){
  
  # 1. formation of the QTL incidence matrix
  ###########################################
  
  ### 2.2 QTL position
  
  QTL <- IncMat_QTL(x = x, mppData = mppData, cross.mat = cross.mat,
                    par.mat = par.mat, par.clu = par.clu, Q.eff = Q.eff,
                    order.MAF = TRUE)
  
  ref.name <- colnames(QTL)
  Env_name <- rep(paste0("_E", 1:nEnv), each = length(ref.name))
  ref.name <- paste0(c(ref.name, ref.name), Env_name)
  
  QTL <- diag(nEnv) %x% QTL
  colnames(QTL) <- ref.name
  
  QTL.el <- dim(QTL)[2] # number of QTL elements
  
  ref.name <- colnames(QTL)
  
  ### 2.1 cofactors
  
  cof.mat <- do.call(cbind, cof.list[which(cof.part[x, ])])
  
  # test if no cofactors
  
  if(is.null(cof.mat)){ cof.mat <- rep(0, length(TraitEnv)); cof.el <- 1
  
  } else {
    
    cof.mat <- diag(nEnv) %x% cof.mat
    
    cof.el <- dim(cof.mat)[2]
    
  }
  
  
  # 2. model computation
  ######################
  
  ### 2.1 homogeneous residual variance error
  
  if(VCOV == "h.err"){
    
    model <- tryCatch(expr = lm(TraitEnv ~ - 1 + CrMatEnv + cof.mat + QTL),
                      error = function(e) NULL)
    
    if (is.null(model)){ 
      
      if(plot.gen.eff) {
        
        if(Q.eff == "cr"){ results <- c(0, rep(1, nEnv * mppData$n.cr))
        
        } else if (Q.eff == "biall") { results <- c(0, rep(1, nEnv))
        
        } else { results <- c(0, rep(1, nEnv * mppData$n.par)) }
        
      } else { results <- 0 }
      
    } else {
      
      
      if(!("QTL" %in% rownames(anova(model)))){ # QTL effect could not be
        # estimated probably due to singularities.
        
        if(plot.gen.eff) {
          
          if(Q.eff == "cr"){ results <- c(0, rep(1, nEnv * mppData$n.cr))
          
          } else if (Q.eff == "biall") { results <- c(0, rep(1, nEnv))
          
          } else { results <- c(0, rep(1, nEnv * mppData$n.par)) }
          
        } else { results <- 0 }
        
      } else {
        
        results <- -log10(anova(model)$Pr[which(rownames(anova(model))=="QTL")])
        
        if(plot.gen.eff){
          
          gen.eff <- QTL_pval_GE(mppData = mppData, nEnv = nEnv, model = model,
                                 Q.eff = Q.eff, x = x, par.clu = par.clu)
          
          
          results <- c(results, gen.eff)
          
        }
        
      }
      
    }
    
    ### 2.2 HRT REML or cross-specific variance residual terms
    
  } else if (VCOV == "Env.err"){
    
    EnvInd <- rep(paste0("_E", 1:nEnv), each = length(mppData$geno.id))
    CrEnvInd <- paste0(rep(mppData$cross.ind, nEnv), EnvInd)
    
    dataset <- data.frame(cof.mat = cof.mat, QTL = QTL,
                          Env = factor(EnvInd, levels = unique(EnvInd)),
                          cr.mat = factor(CrEnvInd, levels = unique(CrEnvInd)),
                          trait = TraitEnv)
    
    colnames(dataset) <- c(paste0("cof", 1:cof.el), paste0("Q",1:QTL.el),
                           "Env", "cr.mat", "trait")
    
    formula.QTL <- paste("+", paste0("Q", 1:QTL.el), collapse = " ")
    formula.fix <- paste("trait ~ -1 + cr.mat + grp(cof)", formula.QTL)
    formula.R <- "~at(Env):units" 
    
    model <- tryCatch(expr = asreml(fixed = as.formula(formula.fix),
                                    rcov =  as.formula(formula.R),
                                    group = list(cof=1:cof.el),
                                    data=dataset, trace = FALSE,
                                    na.method.Y = "omit",
                                    na.method.X = "omit"),
                      error = function(e) NULL)
    
    ### 2.3 random pedigree + HVRT or + CSRT
    
  }  
  
  # 3. organise the results for the mixed models similar for all models
  #####################################################################
  
  if(VCOV != "h.err"){
    
    if (is.null(model)){ 
      
      if(plot.gen.eff) {
        
        if(Q.eff == "cr"){ results <- c(0, rep(1, nEnv * mppData$n.cr))
        
        } else if (Q.eff == "biall") { results <- c(0, rep(1, nEnv))
        
        } else { results <- c(0, rep(1, nEnv * mppData$n.par)) }
        
      } else { results <- 0 }
      
    } else {
      
      W.stat <- sum(wald(model)[3:(QTL.el+2), 3])
      
      if(W.stat == 0){
        
        if(plot.gen.eff) {
          
          if(Q.eff == "cr"){ results <- c(0, rep(1, nEnv * mppData$n.cr))
          
          } else if (Q.eff == "biall") { results <- c(0, rep(1, nEnv))
          
          } else { results <- c(0, rep(1, nEnv * mppData$n.par)) }
          
        } else { results <- 0 }
        
      } else {
        
        df <- sum(wald(model)[3:(QTL.el+2), 1])
        
        pval <- pchisq(W.stat, df, lower.tail = FALSE)
        
        results <- -log10(pval)
        
        if(plot.gen.eff){
          
          gen.eff  <- QTL_pval_mix_GE(model = model, Q.eff = Q.eff, QTL.el = QTL.el,
                                      x = x, par.clu = par.clu,
                                      ref.name = ref.name,
                                      par.names = mppData$parents, fct = "CIM")
          
          results  <- c(results, gen.eff)
          
        }
        
      }
      
    }
    
  }
  
  return(results)
  
}