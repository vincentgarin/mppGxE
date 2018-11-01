#####################
# mod_formulas_oneS #
#####################

# Function to select the random and rcov formula for the mixed model one stage
# analysis

mod_formulas_oneS <- function(VCOV, exp_des_form){

  if (VCOV == 'CSRT'){

    formula.random <- paste0('~ ', exp_des_form)
    formula.rcov <- "~ at(cross_env):units"

  } else if (VCOV == 'CS_CSRT'){

    formula.random <- paste0('~ genotype + ', exp_des_form)
    formula.rcov <- "~ at(cross_env):units"

  } else if(VCOV == 'CS_AR1xAR1'){

    formula.random <- paste0('~ genotype + ', exp_des_form)
    formula.rcov <- "~ at(env):ar1(col):ar1(row)"

  } else if(VCOV == 'CS_CSRT_AR1xAR1'){

    formula.random <- paste0('~ genotype + diag(cross_env):units + ', exp_des_form)
    formula.rcov <- "~ at(env):ar1(col):ar1(row)"

  }

  formulas <- c(formula.random, formula.rcov)

  return(formulas)

}
