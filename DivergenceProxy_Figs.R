DivergenceProxy_Figs <- function(A_list, Collect_list, Troph_list, Omni_list){
  # For a sample of stable food webs:
  # - Cascade divergence: measures the ratio of cascade-classic to cascade-net for each food chain
  # - Proxy 1: measures the network integration of each chain (ratio of interactions within the chain to interactions between the chain and the rest of the network)
  # - Proxy 2: measures the omnivory of each trophic chain (the average omnivory of the chain's constituent species)
  # - Proxy 3: Measures intraguild predation for each trophic chain
  # - Outputs three figures relating the cascade divergence metric to each of the three proxies.

  outputs_norm_all <- data.frame(Collect = numeric(0), RatioChain = numeric(0), Ratio = numeric(0),
                                 OmniChain = numeric(0), IGP = numeric(0))
  
  # Collectivity interval for normalization
  UniqueRoundCollect <- sort(unique(round(unlist(Collect_list), 1)))
  IdxEachCollect <- list()
  for (i in 1:length(UniqueRoundCollect)){
    collect <- UniqueRoundCollect[[i]]
    Idx <- which(round(unlist(Collect_list),1) == collect)
    IdxEachCollect[[i]] <- Idx
  }
  # For each group of food webs within the same Collectivity interval
  for (j in 1:length(UniqueRoundCollect)){
    A_sublist <- A_list[IdxEachCollect[[j]]]
    Collect_sublist <- Collect_list[IdxEachCollect[[j]]]
    Troph_sublist <- Troph_list[IdxEachCollect[[j]]]
    Omni_sublist <- Omni_list[IdxEachCollect[[j]]]
    outputs <- data.frame(Collect = numeric(0), RatioChain = numeric(0), Ratio = numeric(0),
                          OmniChain = numeric(0), IGP = numeric(0))
    # For each food web
    for (i in 1:length(A_sublist)){
      # For each chain with cascade (positive second-order interactions from predator to resource)
      # compute direct and net cascades
      Direct <- A_sublist[[i]]
      diag(Direct) <- 0
      A2 <- Direct%^%2
      Trophic <- Troph_sublist[[i]]
      # compute net effects
      I <- diag(1, nrow = nrow(Direct))
      Net <- solve(I - Direct) # -(-I + Direct)^-1 := (I - Direct)^-1
      MaxTroph <- max(Trophic)
      IdxMaxTroph <- which(Trophic == MaxTroph)
      Minus2Troph <- which(Trophic == (MaxTroph - 2))
      IGP_all <- 0
      for (top in IdxMaxTroph){
        for (bottom in Minus2Troph){
          if ((A2[bottom, top] + Direct[bottom, top]) > 0){ # if trophic cascade (despite possible omnivory)
            DirectCascade <- A2[bottom, top] + Direct[bottom, top]
            NetCascade <- Net[bottom, top]
            # Ratio Net/Direct, which show if there is, and how much, amplification, attenuation, or inversion
            Ratio <- NetCascade/DirectCascade
            # Find middle node(s)
            middle <- which(Direct[top, ] > 0 & Direct[bottom, ] < 0) # consumer = positive effect on pred & negative effect on resource
            # Interaction within chain
            WithinChain <- sum(abs(Direct[c(middle, bottom), top])) + sum(abs(Direct[c(top, bottom), middle])) +
              sum(abs(Direct[c(middle, top), bottom]))
            # Interaction between chain and foodweb
            IdxChain <- c(bottom, middle, top)
            OutOfChain <- sum(abs(Direct[IdxChain, which(!(1:nrow(Direct) %in% IdxChain))])) +
              sum(abs(Direct[which(!(1:nrow(Direct) %in% IdxChain)), IdxChain]))
            # Ratio : integration of the food chain in the food web
            RatioChain <- OutOfChain / WithinChain
            if (RatioChain == Inf){
              RatioChain = NaN
            }
            # Omnivory of the species in the chain
            flow <- -1 * Direct * (Direct < 0)
            flow[is.na(flow)] <- 0
            OmniAll <- TrophInd(flow)$OI
            OmniChain <- (OmniAll[top] + OmniAll[bottom] + sum(OmniAll[middle]))/
              (2 + length(middle)) # divide by the nbr of species in the chain
            # IGP : intraguild predation
            IGP <- 0
            for (j in 1:length(middle)){
              IGP_idx <- which(Direct[top, ] > 0 & Direct[middle[[j]], ] < 0 & Direct[bottom, ] < 0)
              IGP <- IGP + sum(abs(Direct[top, IGP_idx]) + abs(Direct[bottom, IGP_idx]) + abs(Direct[middle[[j]], IGP_idx]))
            }
            # Outputs
            Collect <- Collect_sublist[[i]]
            outputs <- rbind(outputs, list(Collect, RatioChain, Ratio, OmniChain, IGP))
          }
        }
      }
    }
    colnames(outputs) <- c("Collect", "RatioChain", "Ratio", "OmniChain", "IGP")
    # Normalize proxies
    normalize <- function(x) {
      (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    }
    outputs$RatioChain <- normalize(outputs$RatioChain) 
    outputs$OmniChain <- normalize(outputs$OmniChain)
    outputs$IGP <- normalize(outputs$IGP) 
    outputs_norm_all <- rbind(outputs_norm_all, outputs)
  }
  
  # Figures
  
  # Food chain integration
  RatioChain <- ggplot(outputs_norm_all, aes(x = RatioChain, y = Ratio)) +
    geom_point(alpha = 0.05, size = 2) +
    geom_smooth(se = TRUE, method = "lm", linewidth = 0.8) +
    stat_poly_line() + stat_poly_eq() +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_y_continuous(trans = pseudolog10_trans, breaks = c(-10^2, -10^1, -10^0, 0, 10^0, 10^1, 10^2)) +
    labs(x = "Ratio (normalized) \n interactions with the foodweb/interactions within chain", y = "Ratio \n net cascade/n-step causality cascade") +
    theme(axis.title.x = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          axis.title.y = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          legend.text = element_text(size = 12, family = "LM Roman 10"),
          panel.background = element_rect(fill = 'white', color = 'black'),
          # panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
          legend.background = element_rect(fill = 'white', color = 'white'),
          legend.box.background = element_rect(fill = 'white', color = 'white'),
          legend.key = element_rect(fill = "transparent"))
  
  # Omnivory 
  OmniChain <- ggplot(outputs_norm_all, aes(x = OmniChain, y = Ratio)) +
    geom_point(alpha = 0.05, size = 2) +
    geom_smooth(se = TRUE, method = "lm", linewidth = 0.8) +
    stat_poly_line() + stat_poly_eq() +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_y_continuous(trans = pseudolog10_trans, breaks = c(-10^2, -10^1, -10^0, 0, 10^0, 10^1, 10^2)) +
    labs(x = "Omnivory (normalized) of the chain's species", y = "Ratio \n net cascade/n-step causality cascade") +
    theme(axis.title.x = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          axis.title.y = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          legend.text = element_text(size = 12, family = "LM Roman 10"),
          panel.background = element_rect(fill = 'white', color = 'black'),
          # panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
          legend.background = element_rect(fill = 'white', color = 'white'),
          legend.box.background = element_rect(fill = 'white', color = 'white'),
          legend.key = element_rect(fill = "transparent"))
  
  # Intraguild predation
  IGP <- ggplot(outputs_norm_all, aes(x = IGP, y = Ratio)) +
    geom_point(alpha = 0.05, size = 2) +
    geom_smooth(se = TRUE, method = "lm", linewidth = 0.8) +
    stat_poly_line() + stat_poly_eq() +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_y_continuous(trans = pseudolog10_trans, breaks = c(-10^2, -10^1, -10^0, 0, 10^0, 10^1, 10^2)) +
    labs(x = "Intra-guild predation (normalized) in the chain", y = "Ratio \n net cascade/n-step causality cascade") +
    theme(axis.title.x = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          axis.title.y = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          legend.text = element_text(size = 12, family = "LM Roman 10"),
          panel.background = element_rect(fill = 'white', color = 'black'),
          # panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
          legend.background = element_rect(fill = 'white', color = 'white'),
          legend.box.background = element_rect(fill = 'white', color = 'white'),
          legend.key = element_rect(fill = "transparent"))
  
  outputs <- list(RatioChain, OmniChain, IGP)
  return(outputs)
}
