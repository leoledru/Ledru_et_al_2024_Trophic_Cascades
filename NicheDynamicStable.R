NicheDynamicStable <- function(S, C, tmax, tstep, SelfReg, Omni = TRUE, TriTroph = FALSE){
  Alpha <- 0.9
  # SelfReg = 1
  TP = 3
  cond1 <- FALSE
  cond2 <- FALSE
  cond3 <- FALSE
  while (!cond1 & !cond2 & !cond3){
    # Generate an interaction matrix A from niche-model #
    #####################################################
    # niche <- NicheModelConstraint(S, C, Alpha, SelfReg = SelfReg, TP = TP, SaveNetwork = T)
    niche <- NicheModelSimple(S, C, ConvCoeff = 0.5, SelfReg = SelfReg, MaxTL = 3)
    A <- niche[["direct"]]
    if (Omni == FALSE){ # remove omnivory links
      a <- A
      diag(a) <- 0
      flow <- -1 * a * (a < 0)
      flow[is.na(flow)] <- 0
      Troph <- round(TrophInd(flow)$TL)
      IdxTop <- which(Troph == max(Troph))
      IdxMinus2Troph <- which(Troph == (max(Troph)-2))
      A[IdxMinus2Troph, IdxTop] <- 0
      A[IdxTop, IdxMinus2Troph] <- 0
    }
    if (TriTroph == TRUE){ # when we only want a trophic chain with a fixed attack rate
      A[A > 0] <- Alpha * 0.5 # convcoeff
      A[A < 0] <- - Alpha
      diag(A) <- -1 # reset diag
    }
    TrophInit <- niche[["trophic"]]
    # Generalized Lotka-Volterra run 1 #
    ####################################
    r <- rep(0.1, length(TrophInit)) # Vecteur de taux de croissance pour les espèces basales
    N0 <- rep(0.1, length(TrophInit)) # Vecteur d'état initial
    # Paramètres de la simulation
    parameters <- list(A = A, r = r, Troph = TrophInit)
    # Objet de simulation
    ode_system <- ode(y = N0, times = seq(0, tmax, by = tstep), func = LotkaVolterraGeneralized, parms = parameters)
    # Création d'un dataframe à partir des résultats de la simulation
    result_df <- as.data.frame(ode_system)
    # Supprimer la première colonne (temps)
    time <- result_df[1]
    result_df <- result_df[-1]
    run1 <- result_df # save densities for plot
    # Remove extinct species
    final_densities <- result_df[nrow(result_df),]
    IdxExtinct <- which(final_densities < 10^-3)
    A <- A[final_densities>10^-3,final_densities>10^-3]
    a <- A
    diag(a) <- 0
    flow <- -1 * a * (a < 0)
    flow[is.na(flow)] <- 0
    Troph <- round(TrophInd(flow)$TL)
    # Is there still a max troph level = 3 ?
    cond1 <- max(Troph) == 3
    if (cond1){
      # Generalized Lotka-Volterra run 2 #
      ####################################
      r <- rep(0.1, length(Troph))
      N0 <- final_densities[final_densities>10^-3] # Vecteur d'état initial = final state of first run
      # Paramètres de la simulation
      parameters <- list(A = A, r = r, Troph = Troph)
      # Objet de simulation
      ode_system <- ode(y = N0, times = seq(0, tmax, by = tstep), func = LotkaVolterraGeneralized, parms = parameters)
      result_df <- as.data.frame(ode_system)
      time <- result_df[1]
      result_df <- result_df[-1]
      run2 <- result_df # save densities for plot
      for (idx in IdxExtinct){ # add columns with zeros for extinct species in run 1
        name = as.character(idx)
        run2 <- add_column(run2, .before = idx, .value = 0)
      }
      colnames(run2) <- seq(1, S, 1) # reset correct names
      final_densities <- result_df[nrow(result_df),]
      cond2 <- length(final_densities[final_densities>10^-3]) == length(N0) # check that no new extinction -> stable system
      A <- A[final_densities>10^-3,final_densities>10^-3]
      a <- A
      diag(a) <- 0
      flow <- -1 * a * (a < 0)
      flow[is.na(flow)] <- 0
      Troph <- round(TrophInd(flow)$TL)
      # Is there still a max troph level = 3 ?
      cond3 <- max(Troph) == 3
    }
  }
  collect <- spectralRadius(A + diag(1, nrow = nrow(A)))
  out <- list(A = A, N0 = final_densities, run1 = run1, run2 = run2, 
              TrophInit = TrophInit, Troph = Troph, Collect = collect)
  return(out)
}









