#######################
# formula_backward_GE #
#######################

# Write the formula for MPP GxE backward elimination

formula_backward_GE <- function(Q.names, VCOV){


  Q.vect <- vector("list", length = length(Q.names))

  for(i in seq_along(Q.vect)) Q.vect[[i]] <- c(Q.names[-i], Q.names[i])

  if(VCOV == "ID"){

    fbegin <- "trait ~ - 1 + CrMatEnv +"

  } else {

    fbegin <- "trait ~ -1 + check + env:cross +"

  }


  formula.fct <- function(x, fbegin, VCOV) {

    if(VCOV == "ID"){

      paste(fbegin, paste(x, collapse = "+"))

    } else {

      QTL.el <- vapply(X = x, FUN = function(x) paste0("grp(", x, ")"),
                       character(1))
      paste(fbegin, paste(QTL.el, collapse = "+"))

    }




  }

  lapply(X = Q.vect, formula.fct, fbegin = fbegin, VCOV = VCOV)

}
