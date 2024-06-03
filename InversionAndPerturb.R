InversionAndPerturb <- function(tmax, tstep, S, C, Omni = FALSE){
  
  # Generate a stable system with max trophic level = 3 #
  out <- NicheDynamicStable(S = S, C = C, tmax = tmax, tstep = tstep, SelfReg = 1, 
                            Omni = Omni, TriTroph = FALSE)
  A <- out[["A"]]
  Collect <- round(out[["Collect"]],1)
  S <- nrow(A)
  C <- sum(A != 0) / S^2
  StableRuns <- rbind(out[["run1"]], out[["run2"]])
  Equilibrium <- StableRuns[(dim(StableRuns)[[1]] - 2000):dim(StableRuns)[[1]],]
  Equilibrium <- Equilibrium[, sapply(Equilibrium, function(x) tail(x, 1) >= 10^-3)]
  # Identify the trophic chain(s) with inversion #
  I <- diag(1, nrow = nrow(A))
  Abis <- A
  diag(Abis) <- 0
  A2 <- Abis%^%2
  AInv <- solve(I - Abis)
  outInv <- InversionEachChain(A)
  # If there is chain(s) with inversion, perturb them and vizualise
  Inversion <- outInv[["Inversion"]]
  Idx <- which(!sapply(Inversion, is.null))
  if (length(Idx) > 0){ # if there at least on chain with inversion
    EqN <- out[["run2"]]
    EqN <- EqN[nrow(EqN),]
    EqTroph <- out[["Troph"]] # Trophy des species non-null
    EqN <- EqN[EqN >= 10^-3] # Densités à l'eq pour état initial
    r <- rep(0.1, length(EqTroph))
    IdxBottom <- unlist(outInv["IdxBottom"]) # index of species at the bottom of chain
    plot_all_species <- vector("list", length(Idx))
    plot_bottom_species <- vector("list", length(Idx))
    plot_variation_eq <- vector("list", length(Idx))
    plot_variation_eq_1 <- vector("list", length(Idx))
    FullRuns <- list()
    IdxBottoms <- list()
    for (i in 1:length(Idx)){ # for each top species of chain with inversion
      # Perturbation
      idx_perturb <- Inversion[[Idx[[i]]]][["top"]] # top node à perturber
      idx_bottom <- Inversion[[Idx[[i]]]][["bottom"]] # bottom nodes qui répondent
      # calcul le r_perturb maximal possible sans extinction et divise par deux
      k <- which(AInv[,idx_perturb] < 0)
      r_perturb <- min(EqN[k]/-AInv[k,idx_perturb]) / 2 
      # find middle node
      edges <- which(Abis < 0, arr.ind = TRUE)
      preys <- edges[,1] # nodes which suffer a negative effect
      preds <- edges[,2] # nodes which add a negative effect on others
      preys_of_top <- preys[preds == idx_perturb] # nodes which suffer a negative effect from the top of the focus chain
      pred_of_bottom <- preds[preys %in% idx_bottom] # nodes which add a negative effect on the bottom(s) of the focus chain
      middle_node <- unique(preds[preds %in% preys_of_top & preds %in% pred_of_bottom])
      # Compute the Eq if only direct cascade and the Eq with net effects : ONLY FOR THE PERTURBED CHAIN
      # EqDirectBottom <- EqN[[Bottom]] + DirectCascade*r_perturb
      EqNChain <- c(EqN[idx_perturb], sum(EqN[middle_node]), sum(EqN[idx_bottom]))
      EqDirect <- EqN + (I[,idx_perturb] + Abis[,idx_perturb] + A2[,idx_perturb])*r_perturb
      EqDirectTop <- EqDirect[idx_perturb]
      EqDirectTop_1 <- EqN[idx_perturb] + (I[idx_perturb,idx_perturb] + Abis[idx_perturb,idx_perturb])*r_perturb
      EqDirectMiddle <- sum(EqDirect[middle_node])
      EqDirectMiddle_1 <- sum(EqN[middle_node] + (I[middle_node,idx_perturb] + Abis[middle_node,idx_perturb])*r_perturb)
      EqDirectBottom <- sum(EqDirect[idx_bottom])
      EqDirectChain <- c(EqDirectTop, EqDirectMiddle, EqDirectBottom)
      EqDirectChain_1 <- c(EqDirectTop_1, EqDirectMiddle_1, EqDirectBottom)
      EqNet <- EqN + AInv[,idx_perturb]*r_perturb
      EqNetTop <- EqNet[idx_perturb]
      EqNetMiddle <- sum(EqNet[middle_node])
      EqNetBottom <- sum(EqNet[idx_bottom])
      EqNetChain <- c(EqNetTop, EqNetMiddle, EqNetBottom)

      plot_bottom_species[[i]] <- vector("list", length(idx_bottom))
      parameters <- list(A = A, r = r, Troph = EqTroph, r_perturb = r_perturb, idx_perturb = idx_perturb)
      tmax = 2000
      out_perturb <- ode(y = EqN, times = seq(0, tmax, by = tstep), 
                         func = LotkaVolterraGeneralizedPerturb, parms = parameters)
      out_perturb <- as.data.frame(out_perturb)
      out_perturb <- out_perturb[-1]
      new_df <- data.frame(Time = numeric(0), Species = numeric(0), Troph = numeric(0),
                           Type = character(0), Density = numeric(0))
      # time <- seq(1, nrow(out_perturb) + nrow(Equilibrium), 1)
      time <- seq(0, tmax + dim(Equilibrium)[[1]]*tstep, tstep)
      # Print the cascade values
      NetCascade <- outInv[["CascadeNet"]][[Idx[[i]]]][which(IdxBottom %in% idx_bottom)]
      DirectCascade <- outInv[["CascadeDirect"]][[Idx[[i]]]][which(IdxBottom %in% idx_bottom)]
      print(paste0("Direct = ", DirectCascade, " Net = ", NetCascade))
      # Réorganiser les données
      for (j in 1:dim(out_perturb)[[2]]) {
        species_num <- j
        if (j %in% idx_perturb){
          species_type <- "Top"
        }else if (j %in% idx_bottom){
          species_type <- "Bottom"
        }else if (j %in% middle_node){
          species_type <- "Middle"
        }
        else{
          species_type <- "Not in perturbed chain"
        }
        species_density <- c(Equilibrium[, j], out_perturb[, j])
        new_data <- data.frame(Time = time, Species = species_num, Troph = EqTroph[j], Type = species_type, Density = species_density)
        new_df <- rbind(new_df, new_data)
      }
      FullRuns[[i]] <- new_df
      IdxBottoms[[i]] <- idx_bottom
      # Créer le graphique avec toutes les espèces
      # plot_all_species[[i]] <- ggplot(new_df, aes(x = Time, y = Density, color = Type, group = Species)) +
      #   geom_line() +
      #   labs(title = paste0("Response to perturbation, S = ", S, " C = ", C, " Collectivity = ", Collect),
      #        x = "Time", y = "Density") +
      #   theme_minimal()
      
      # Intercepts_df <- data.frame(Intercept = EqDirectBottom, Type = "Bottom")
      Intercepts_df <- data.frame(Intercept = EqDirect[idx_bottom], Type = "Bottom")
      plot_all_species[[i]] <- ggplot(new_df) +
        geom_line(aes(x = Time, y = Density, color = Type, 
                      linetype = "Real dynamics", group = Species), linewidth = 1) +
        geom_segment(data = Intercepts_df, aes(y = Intercept, yend = Intercept, x = 200, xend = Inf, color = Type, 
                                               linetype = "Second-order\nprediction")) +
        labs(title = paste0("Response to top-species perturbation\nwith collectivity = ", Collect), 
             x = "Time", y = "Density") +
        scale_color_discrete(name = "Trophic Levels") +
        scale_color_manual(values = c("red", "green", "black", "blue")) +
        scale_linetype_manual(name = "", values = c("Real dynamics" = "solid", "Second-order\nprediction" = "dashed")) +
        geom_vline(xintercept = 200) +
        theme_minimal()
      
      # Variation de densités en mode pyramide de biomasse
      data <- data.frame(
        direct = (EqDirectChain - EqNChain) / EqNChain * 100,
        net = (EqNetChain - EqNChain) / EqNChain * 100,
        y = c(3,2,1)
      )
      plot_variation_eq[[i]] <- ggplot(data) +
        geom_rect(aes(ymin = y[[1]]-0.5, ymax = y[[1]]+0.5, xmin = 0, xmax = direct[[1]]), 
                  fill = "red", alpha = 0.2) +
        geom_rect_pattern(aes(ymin = y[[1]]-0.5, ymax = y[[1]]+0.5, xmin = 0, xmax = net[[1]]), 
                          fill = NA, colour = "black", pattern = "stripe", pattern_size = 0.05,
                          pattern_density = 0.1, pattern_spacing = 0.02, pattern_fill = "black") +
        geom_rect(aes(ymin = y[[2]]-0.5, ymax = y[[2]]+0.5, xmin = 0, xmax = direct[[2]]), 
                  fill = "red", alpha = 0.2) +
        geom_rect_pattern(aes(ymin = y[[2]]-0.5, ymax = y[[2]]+0.5, xmin = 0, xmax = net[[2]]), 
                          fill = NA, colour = "black", pattern = "stripe", pattern_size = 0.05,
                          pattern_density = 0.1, pattern_spacing = 0.02, pattern_fill = "black") +
        geom_rect(aes(ymin = y[[3]]-0.5, ymax = y[[3]]+0.5, xmin = 0, xmax = direct[[3]]), 
                  fill = "red", alpha = 0.2) +
        geom_rect_pattern(aes(ymin = y[[3]]-0.5, ymax = y[[3]]+0.5, xmin = 0, xmax = net[[3]]), 
                          fill = NA, colour = "black", pattern = "stripe", pattern_size = 0.05,
                          pattern_density = 0.1, pattern_spacing = 0.02, pattern_fill = "black") +
        labs(title = paste0("Collectivity = ", Collect),
             x = "% Variation between after and before\nperturbation equilibriums",
             y = "Trophic level") +
        theme(axis.title.x = element_text(size = 14, face = "bold", family = "LM Roman 10"),
              axis.title.y = element_text(size = 14, face = "bold", family = "LM Roman 10"),
              plot.title = element_text(size = 14, face = "bold", family = "LM Roman 10"),
              panel.background = element_rect(fill = 'white', color = 'black'))
      
      data <- data.frame(
        direct = (EqDirectChain_1 - EqNChain) / EqNChain * 100,
        net = (EqNetChain - EqNChain) / EqNChain * 100,
        y = c(3,2,1)
      )
      plot_variation_eq_1[[i]] <- ggplot(data) +
        geom_rect(aes(ymin = y[[1]]-0.5, ymax = y[[1]]+0.5, xmin = 0, xmax = direct[[1]]), 
                  fill = "blue", alpha = 0.2) +
        geom_rect_pattern(aes(ymin = y[[1]]-0.5, ymax = y[[1]]+0.5, xmin = 0, xmax = net[[1]]), 
                          fill = NA, colour = "black", pattern = "stripe", pattern_size = 0.05,
                          pattern_density = 0.1, pattern_spacing = 0.02, pattern_fill = "black") +
        geom_rect(aes(ymin = y[[2]]-0.5, ymax = y[[2]]+0.5, xmin = 0, xmax = direct[[2]]), 
                  fill = "blue", alpha = 0.2) +
        geom_rect_pattern(aes(ymin = y[[2]]-0.5, ymax = y[[2]]+0.5, xmin = 0, xmax = net[[2]]), 
                          fill = NA, colour = "black", pattern = "stripe", pattern_size = 0.05,
                          pattern_density = 0.1, pattern_spacing = 0.02, pattern_fill = "black") +
        geom_rect(aes(ymin = y[[3]]-0.5, ymax = y[[3]]+0.5, xmin = 0, xmax = direct[[3]]), 
                  fill = "red", alpha = 0.2) +
        geom_rect_pattern(aes(ymin = y[[3]]-0.5, ymax = y[[3]]+0.5, xmin = 0, xmax = net[[3]]), 
                          fill = NA, colour = "black", pattern = "stripe", pattern_size = 0.05,
                          pattern_density = 0.1, pattern_spacing = 0.02, pattern_fill = "black") +
        labs(title = paste0("Collectivity = ", Collect),
             x = "% Variation between after and before\nperturbation equilibriums",
             y = "Trophic level") +
        theme(axis.title.x = element_text(size = 14, face = "bold", family = "LM Roman 10"),
              axis.title.y = element_text(size = 14, face = "bold", family = "LM Roman 10"),
              plot.title = element_text(size = 14, face = "bold", family = "LM Roman 10"),
              panel.background = element_rect(fill = 'white', color = 'black'))
      
      # Créer le graphique avec seulement les espèces de Type "Bottom"
      # for (j in 1:length(idx_bottom)){
      #   Direct <- A[idx_bottom[[j]], idx_perturb]
      #   Order2 <- A
      #   diag(Order2) <- 0
      #   Order2 <- Order2%^%2
      #   Order2 <- Order2[idx_bottom[[j]], idx_perturb]
      #   plot_bottom_species[[i]][[j]] <- ggplot(filter(new_df, Type == "Bottom" & Species == idx_bottom[[j]]), aes(x = Time, y = Density, group = Species)) +
      #     geom_line(color = "red") +
      #     annotate(geom = 'text', label = paste0("direct effect (omnivory) A[bottom, top] = ", round(Direct,3),
      #              "\n", "cascade effect A²[bottom, top] = ", round(Order2,3),
      #              "\n", "long-term effect A^-1[bottom, top} = ", round(NetCascade[[j]],3)),
      #              x = 500, y = -Inf+0.008, hjust = 0, vjust = 0) +
      #     labs(title = "Zoom on species at the perturbed chain bottom", x = "Time", y = "Density") +
      #     theme_minimal()
      # }
      # Disposer les deux graphiques côte à côte
      # grid.arrange(plot_all_species, plot_bottom_species, ncol = 2)
    }
    OutEnd <- list(plot_all_species = plot_all_species, 
                   # plot_bottom_species = plot_bottom_species,
                   plot_variation_eq = plot_variation_eq,
                   plot_variation_eq_1 = plot_variation_eq_1,
                   A = out[["A"]], OutInv = outInv, S = S, C = C, FullRuns = FullRuns, IdxBottoms = IdxBottoms)
    return(OutEnd)
  }else{
    print("No chain with inversion")
  }
}
