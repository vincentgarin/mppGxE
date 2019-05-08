##########################
# Qeff_res_processing_GE #
##########################

# Function to process the result of a MPP GxE QTL effect estimation model

Qeff_res_processing_GE <- function(model, mppData, Q.eff, VCOV, names.QTL,
                                   Q.pos, nQTL, n_al, nEnv){

  nQTLel <- length(names.QTL)

  if(VCOV == "ID"){

    results <- summary(model)$coefficients
    index <- (substr(rownames(results), 1, 1) == "Q")
    results <- subset(x = results, subset = index, drop = FALSE)
    row.names <- substr(rownames(results), 5, nchar(rownames(results)))
    rownames(results) <- row.names


  } else {

    index <- substr(names(rev(model$coefficients$fixed)), 1, 1) == "Q"

    w.table <- asreml::wald(model)
    w.stat <- w.table[substr(rownames(w.table), 1, 1) == "Q", c(3, 4)]

    results <- cbind(rev(model$coefficients$fixed)[index],
                     rev(sqrt(model$vcoeff$fixed))[index], w.stat)
    results <- as.matrix(results)

  }

  # control for singular values and fill missing values

  ref.mat <- matrix(rep(c(0, 0, 0, 1), nQTLel), nrow = nQTLel,
                    byrow = TRUE)

  index <- match(rownames(results), names.QTL)
  ref.mat[index, ] <- results
  rownames(ref.mat) <- names.QTL

  if(Q.eff == "anc"){

    Q_ind <- rep(paste0('Q', 1:nQTL), 2*n_al)

    ref_mat_temp <- c()

    for(j in 1:nQTL){

      # select the ith QTL

      res_j <- ref.mat[Q_ind %in% paste0("Q", j), ]

      # define the prjection matrix

      A.allele <- mppData$par.clu[Q.pos[j], ]
      A <- model.matrix(~ factor(A.allele) - 1)

      Qmat <- inc_mat_QTL(x = Q.pos[j], mppData = mppData, Q.eff = Q.eff,
                          order.MAF = TRUE)

      anc_ref <- as.numeric(substr(x = colnames(Qmat), start = 9,
                                   stop = nchar(colnames(Qmat))))

      row_ord <- data.frame(c1 = anc_ref, c2 = 1:length(anc_ref))
      row_ord <- row_ord[order(row_ord$c1), ]
      row_ord <- row_ord$c2

      # split into each environment

      n_tot <- n_al[j] * nEnv
      s_env <- 1:n_tot
      env_id <- split(s_env, ceiling(seq_along(s_env)/(n_tot/nEnv)))

      for(k in 1:nEnv){

        res_jk <- res_j[env_id[[k]], ]

        res_jk <- A %*% res_jk[row_ord, ]
        rownames(res_jk) <- paste0(paste0('Q', j), mppData$parents, '_E', k)

        ref_mat_temp <- rbind(ref_mat_temp, res_jk[order(res_jk[, 1]), ])

      }

    }

    ref.mat <- ref_mat_temp

  }

  # add sign stars

  Sign <- sapply(ref.mat[, 4], FUN = sign.star)
  results <- data.frame(ref.mat, Sign, stringsAsFactors = FALSE)

  # add column names

  if(VCOV == "ID"){ col.names <- c("Effect", "Std.Err", "t-test", "p-value")

  } else {col.names <- c("Effect", "Std.Err", "W-stat", "p-value")}

  colnames(results)[1:4] <- col.names

  # split per QTL

  if(Q.eff == "anc"){

    Q_lab <- paste0("Q", 1:nQTL)

    Q_f <- rep(Q_lab, each = mppData$n.par*nEnv)

    Q_f <- factor(Q_f, levels = Q_lab)


  } else {

    Q_f <- factor(rep(paste0("Q", 1:nQTL), each = (nQTLel/nQTL)),
                  levels = paste0("Q", 1:nQTL))

  }

  Qeff.mat <- split(x = results, f = Q_f)

  # split within QTL per environments

  Qeff <- vector(mode = 'list', length = nQTL)

  for(z in 1:nQTL){

    Qeff_z <- Qeff.mat[[z]]

    E_lab <- paste0('E', 1:nEnv)
    E_f <- rep(E_lab, each = dim(Qeff_z)[1]/nEnv)
    E_f <- factor(E_f, levels = E_lab)

    Qeff[[z]] <- split(x = Qeff_z, f = E_f)

  }

  names(Qeff) <- paste0('Q', 1:nQTL)

  return(Qeff)


}
