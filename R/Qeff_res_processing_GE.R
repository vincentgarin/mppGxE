##########################
# Qeff_res_processing_GE #
##########################

# Function to process the result of a MPP GxE QTL effect estimation model

Qeff_res_processing_GE <- function(model, mppData, Q.eff, VCOV, names.QTL,
                                   nQTL, n_al, nEnv){

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

  # if(Q.eff == "anc"){
  #
  #   # The results could/should be projected into parents within env
  #   # do it later for nicer presentation. Now focus on CV.
  #
  #   }

  # add sign stars

  Sign <- sapply(ref.mat[, 4], FUN = sign.star)
  results <- data.frame(ref.mat, Sign, stringsAsFactors = FALSE)

  # add column names

  if(VCOV == "ID"){ col.names <- c("Effect", "Std.Err", "t-test", "p-value")

  } else {col.names <- c("Effect", "Std.Err", "W-stat", "p-value")}

  colnames(results)[1:4] <- col.names

  if(Q.eff == "anc"){

    Q_lab <- paste0("Q", 1:nQTL)

    nal_env <- (n_al*nEnv)

    Q_f <- mapply(FUN = function(x, y) rep(x, y), x = Q_lab, y = nal_env)

    Q_f <- factor(unlist(Q_f), levels = Q_lab)


  } else {

    Q_f <- factor(rep(paste0("Q", 1:nQTL), each = (nQTLel/nQTL)),
                  levels = paste0("Q", 1:nQTL))

  }

  Qeff.mat <- split(x = results, f = Q_f)

  return(Qeff.mat)


}
