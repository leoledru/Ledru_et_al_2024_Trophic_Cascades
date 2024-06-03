GenerateStableFwbs <- function(Svec, Cvec, Rep, tmax, tstep, Omni = FALSE, type = "niche"){
  library(fdrtool)
  library(NetIndices)
  library(sparsevar)
  library(deSolve)
  
  Alpha <- 0.5
  SelfReg = 1
  TP = 3
  outputs <- list()
  count <- 0
  for (S in Svec){
    for (C in Cvec){
      for (i in 1:Rep){
        count <- count + 1
        cond1 <- FALSE
        cond2 <- FALSE
        cond3 <- FALSE
        while (!cond1 & !cond2 & !cond3){
          # Generate an interaction matrix A from niche-model #
          ####################################################
          # NICHE MODEL
          if (type == "niche"){
            # condTP <- FALSE
            # while (!condTP){
            # niche <- NicheModelConstraint(S, C, Alpha, SelfReg = SelfReg, TP = TP, SaveNetwork = T, Omni = FALSE)
            niche <- NicheModelSimple(S, C, ConvCoeff = 0.5, SelfReg, MaxTL = 4)
            # condTP <- max(niche[["trophic"]]) == 3
            # }
            A <- niche[["direct"]]
            TrophInit <- niche[["trophic"]]
          }else if (type == "structured"){
            # STRUTURED MODEL
            niche <- StructuredNetwork(S = S, C = C, MaxTroph = 3)
            A <- niche[["NetworkMatrix"]]
            TrophInit <- niche[["NetworkList"]]$TrophicLevel
          }else{
            print("choisir un type de modèle valide : niche ou structured")
            break
          }
          
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
          # Remove extinct species
          final_densities <- result_df[nrow(result_df),]
          IdxExtinct <- which(final_densities < 10^-3)
          A <- A[final_densities>10^-3,final_densities>10^-3]
          a <- A
          if (length(a) > 1){
            diag(a) <- 0
            flow <- -1 * a * (a < 0)
            flow[is.na(flow)] <- 0
            Troph <- round(TrophInd(flow)$TL)
            # Is there still a max troph level >= 3 ?
            cond1 <- max(Troph) >= 3
          }else{
            cond1 <- FALSE
          }
          if (cond1 & length(IdxExtinct)>0){
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
            final_densities <- result_df[nrow(result_df),]
            cond2 <- length(final_densities[final_densities>10^-3]) == length(N0) # check that no new extinction -> stable system
            A <- A[final_densities>10^-3,final_densities>10^-3]
            a <- A
            diag(a) <- 0
            flow <- -1 * a * (a < 0)
            flow[is.na(flow)] <- 0
            Troph <- round(TrophInd(flow)$TL)
            # Is there still a max troph level = 3 ?
            cond3 <- max(Troph) >= 3
          }else if(cond1 & length(IdxExtinct)==0){
            cond2 <- cond1
            cond3 <- cond1
          }
        }
        Abis <- A + diag(1, nrow = nrow(A)) # diag = 0
        connectance <- sum(Abis != 0) / nrow(Abis)^2
        collect <- spectralRadius(Abis)
        outputs[[count]] <- list(S = nrow(A), C = connectance, A = A, 
                                Troph = Troph, Collect = collect)
      }
      print(paste0("S = ", S, " C = ", C))
      # if (type == "niche"){
      #   save(outputs, file = "StableFwbs.Rdata")
      # }else{
      #   save(outputs, file = "StableFwbsStructured.Rdata")
      # }
    }
  }
  return(outputs)
}
