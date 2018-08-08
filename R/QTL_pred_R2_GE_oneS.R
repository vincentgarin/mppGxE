#######################
# QTL_pred_R2_GE_oneS #
#######################

QTL_pred_R2_GE_oneS <- function(plot_data, mppData.ts, mppData.vs,
                                trait = NULL, cv.ref, nEnv, Q.eff = "cr",
                                VCOV = "ID", QTL = NULL, her = 1) {

  if(is.character(QTL)){ n.QTL <- length(QTL) } else { n.QTL <- dim(QTL)[1] }

  # Determine the environments

  EnvNames <- unique(plot_data$env)

  # form the reference trait values (we assume that the user give differ
  # reference for the particular data analysis)

  t_val <- c(mppData.vs$pheno[, cv.ref])
  t_val <- rep(t_val, nEnv)


  # 2. obtain the genetic effects (Betas)
  #######################################

  effects <- QTL_effects_oneS(plot_data = plot_data,
                              mppData = mppData.ts, trait = trait,
                              QTL = QTL, Q.eff = Q.eff, VCOV = VCOV)


  Qeff_names <- lapply(X = effects, FUN = function(x) rownames(x))
  Qeff_names <- unlist(Qeff_names)
  names(Qeff_names) <- NULL

  B.ts <- lapply(X = seq_along(effects), FUN = function(x, Qeff) Qeff[[x]][, 1],
                 Qeff = effects)

  # 3. obtain the QTL incidence matrices of the positions (X.vs)
  ##############################################################

  ######## add the same code as for M3 because no need anymore to project
  ######## at the plot level.

  if(is.character(QTL)){

    Q.pos <- which(mppData.vs$map[, 1] %in% QTL)

    QTL <- mppData.vs$map[mppData.vs$map[, 1] %in% QTL, ]

  } else {

    Q.pos <- which(mppData.vs$map[, 1] %in% QTL[, 1])

  }

  nQTL <- length(Q.pos)

  Q.list <- lapply(X = Q.pos, FUN = inc_mat_QTL, mppData = mppData.vs,
                   Q.eff = Q.eff, order.MAF = TRUE)

  Q.names <- function(x, Q.list, nEnv){
    rep(paste0("Q", x, attr(Q.list[[x]], "dimnames")[[2]]), nEnv)
  }

  names.QTL <- unlist(lapply(X = 1:nQTL, FUN = Q.names, Q.list = Q.list,
                             nEnv = nEnv))

  if(Q.eff == "anc"){

    n_al <- unlist(lapply(X = Q.list, FUN = function(x) dim(x)[2]))

    e_lab <- paste0("E", 1:nEnv)

    Env.names <- lapply(X = n_al, FUN = function(x, e_lab) rep(e_lab, each = x),
                        e_lab = e_lab)

    Env.names <- unlist(Env.names)

  } else {

    n_al <- NULL

    Env.names <- rep(rep(paste0("E", 1:nEnv), each = dim(Q.list[[1]])[2]), nQTL)

  }

  names.QTL <- paste(names.QTL, Env.names, sep = "_")

  Q.list <- lapply(X = Q.list, FUN =  function(x, nEnv) diag(nEnv) %x% x,
                   nEnv = nEnv)

  names(Q.list) <- paste0("Q", 1:length(Q.list))

  # 4. Predicted R squared computation cor(y.vs, X.vs * B.ts)^2
  ##############################################################

  X.vs <- as.matrix(do.call(cbind, Q.list))
  colnames(X.vs) <- names.QTL

  # order X.vs same order as the gentic predictor B.ts

  X.vs <- X.vs[, Qeff_names]

  # use only complete case informations

  dataset <- cbind(t_val, X.vs)
  cross.ind <- rep(mppData.vs$cross.ind, nEnv)
  Env_ind <- rep((paste0("E", 1:nEnv)), each = length(mppData.vs$geno.id))
  index <- complete.cases(dataset)

  y.vs <- dataset[index, 1, drop = FALSE]
  X.vs <- dataset[index, 2:dim(dataset)[2], drop = FALSE]
  cross.ind <- cross.ind[index]
  Env_ind <- Env_ind[index]

  # n.cr <- table(factor(cross.ind, levels = unique(cross.ind)))

  B.ts <- unlist(B.ts)
  B.ts[is.na(B.ts)] <- 0

  y.vs.hat <- X.vs %*% B.ts

  dataset <- data.frame(y.vs, y.vs.hat, cross.ind, Env_ind)

  with.cross.cor <- function(x){
    if((length(unique(x[, 1])) == 1) || (length(unique(x[, 2])) == 1)){ 0
    } else {100 * ((cor(x[, 1], x[, 2])^2))}
  }

  # iterate over env

  R2_av <- rep(0, nEnv)

  for(i in 1:nEnv){

    dataset_i <- dataset[dataset$Env_ind == paste0("E", i), ]

    cr_ind_i <- dataset_i$cross.ind

    dataset.cr_i <- split(x = dataset_i,
                          f = factor(cr_ind_i, levels = unique(cr_ind_i)))

    R2.cr_i <- unlist(lapply(X = dataset.cr_i, FUN =  with.cross.cor))

    R2_av[i] <- mean(R2.cr_i)

  }

  # R2_av: averaged pred R2 over cross within environments

  return(R2_av)

}
