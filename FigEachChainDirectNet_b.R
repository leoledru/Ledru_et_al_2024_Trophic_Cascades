FigEachChainDirectNet_b <- function(A_list, Collect_list, Troph_list, Omni_list){
  ### Version qui travaille directement sur des réseaux déjà générés et stables ###
  
  outputs <- data.frame(Collect = numeric(0), Omni = numeric(0), Type = character(0), Value = numeric(0))
  outputs_b <- data.frame(Collect = numeric(0), Omni = numeric(0), InvTot = numeric(0), PercentInv = numeric(0), 
                          Type = numeric(0), RatioTot = numeric(0))
  percent_inv <- data.frame(Collect = numeric(0), Omni = numeric(0), InvPercent = numeric(0))
  OutPerFwb <- data.frame(Collect = numeric(0), Type = character(0), Percent = numeric(0))

  for (i in 1:length(A_list)){
    # For each chain with cascade (positive order two interaction) compute direct and net cascades
    Direct <- A_list[[i]]
    diag(Direct) <- 0
    A2 <- Direct%^%2
    Trophic <- Troph_list[[i]]
    # Compute Net effects
    I <- diag(1, nrow = nrow(Direct)) # matrice identité
    Net <- solve(I - Direct) # -(-I + Direct)^-1 := (I - A)^-1
    NbrChain <- 0
    NbrInv <- 0
    NbrCascade <- 0
    NbrAtt <- 0
    NbrAmp <- 0
    EachDirect <- list()
    EachNet <- list()
    MaxTroph <- max(Trophic)
    IdxMaxTroph <- which(Trophic == MaxTroph)
    Minus2Troph <- which(Trophic == (MaxTroph - 2))
    for (top in IdxMaxTroph){
      for (bottom in Minus2Troph){
        if ((A2[bottom, top] + Direct[bottom, top]) > 0){ # if trophic cascade (positive order 2 effect, despite possible omnivory)
          NbrChain <- NbrChain + 1
          DirectCascade <- A2[bottom, top] + Direct[bottom, top]
          if (DirectCascade < 0) {
            print("directe négative, anormal")
          }
          NetCascade <- Net[bottom, top]
          
          # DIRECT AND NET
          # outputs <- rbind(outputs, list(Collect_list[[i]], "direct", DirectCascade))
          # outputs <- rbind(outputs, list(Collect_list[[i]], "net", NetCascade))
          # ONLY NET, WIEGHTED BY DIRECT
          Ratio <- NetCascade/DirectCascade
          if (Ratio > 0 & Ratio < 0.8){
            outputs <- rbind(outputs, list(Collect_list[[i]], Omni_list[[i]], "attenuated", Ratio))
            NbrAtt <- NbrAtt + 1
          }else if (Ratio > 1.2){
            outputs <- rbind(outputs, list(Collect_list[[i]], Omni_list[[i]], "amplificated", Ratio))
            NbrAmp <- NbrAmp + 1
          }else if (Ratio < 0){
            outputs <- rbind(outputs, list(Collect_list[[i]], Omni_list[[i]], "inverted", Ratio))
            NbrInv <- NbrInv + 1
          }else{
            outputs <- rbind(outputs, list(Collect_list[[i]], Omni_list[[i]], "classic cascade", Ratio))
            NbrCascade <- NbrCascade + 1
          }
          EachDirect <- append(EachDirect, DirectCascade)
          EachNet <- append(EachNet, NetCascade)
        }
      }
    }
    if (length(EachNet) != 0){
      PercentCascade <- NbrCascade / NbrChain * 100 
      PercentAtt <- NbrAtt / NbrChain * 100 
      PercentAmp <- NbrAmp / NbrChain * 100 
      PercentInv <- NbrInv / NbrChain * 100 
      RatioTot <- sum(unlist(EachNet)) / sum(unlist(EachDirect))
      # Collect <- round(Collect_list[[i]],1)
      # Omni <- round(Omni_list[[i]],1)
      Collect <- Collect_list[[i]]
      Omni <- Omni_list[[i]]
      if (RatioTot > 0 & RatioTot < 0.8){
        outputs_b <- rbind(outputs_b, list(Collect, Omni, 
                                           ifelse(sign(sum(unlist(EachDirect))) != sign(sum(unlist(EachNet))), 1, 0),
                                           sum(NbrInv)/length(NbrInv)*100, "attenuated", RatioTot))
      }else if (RatioTot > 1.2){
        outputs_b <- rbind(outputs_b, list(Collect, Omni, 
                                           ifelse(sign(sum(unlist(EachDirect))) != sign(sum(unlist(EachNet))), 1, 0),
                                           sum(NbrInv)/length(NbrInv)*100, "amplificated", RatioTot))
      }else if (RatioTot < 0){
        outputs_b <- rbind(outputs_b, list(Collect, Omni, 
                                           ifelse(sign(sum(unlist(EachDirect))) != sign(sum(unlist(EachNet))), 1, 0),
                                           sum(NbrInv)/length(NbrInv)*100, "inverted", RatioTot))
      }else{
        outputs_b <- rbind(outputs_b, list(Collect, Omni, 
                                           ifelse(sign(sum(unlist(EachDirect))) != sign(sum(unlist(EachNet))), 1, 0),
                                           sum(NbrInv)/length(NbrInv)*100, "classic cascade", RatioTot))
      }
      percent_inv <- rbind(percent_inv, list(Collect_list[[i]], Omni_list[[i]], PercentInv))
      OutPerFwb <- rbind(OutPerFwb, list(Collect, "classic cascade", PercentCascade))
      OutPerFwb <- rbind(OutPerFwb, list(Collect, "attenuated", PercentAtt))
      OutPerFwb <- rbind(OutPerFwb, list(Collect, "amplificated", PercentAmp))
      OutPerFwb <- rbind(OutPerFwb, list(Collect, "inverted", PercentInv))
    }
  }
  colnames(outputs) <- c("Collect", "Omni", "Type", "Value")
  colnames(outputs_b) <- c("Collect", "Omni", "InvTot", "InvPercent", "Type", "RatioTot")
  colnames(percent_inv) <- c("Collect", "Omni", "InvPercent")
  colnames(OutPerFwb) <- c("Collect", "Type", "Percent")
  # Réorganiser les niveaux de la variable de couleur dans l'ordre souhaité
  outputs$Type <- factor(outputs$Type, levels = c("classic cascade", "amplificated", "inverted", "attenuated"))
  outputs_b$Type <- factor(outputs_b$Type, levels = c("classic cascade", "amplificated", "inverted", "attenuated"))
  
  # Figure Percentage of Types of cascade for each foodweb, aggregated by collect
  ## Définir les intervalles pour "Collect"
  CollectInterval <- seq(0, 2, by = 0.1)
  ## Agréger les données
  OutPerFwbMean <- OutPerFwb %>%
    mutate(CollectInterval = cut(Collect, breaks = CollectInterval)) %>%
    group_by(CollectInterval, Type) %>%
    summarise(Mean = mean(Percent), Sd = sd(Percent)) %>%
    ungroup()
  OutPerFwbMean$Sd[is.na(OutPerFwbMean$Sd)] <- 0
  PerFwbMean <- ggplot(OutPerFwbMean, aes(x = CollectInterval, y = Mean, group = Type, color = Type)) +
    geom_line(linewidth = 1.5) +
    geom_ribbon(aes(x = CollectInterval, ymin = pmax(Mean - Sd, 0), ymax = pmin(Mean + Sd, 100), fill = Type),
                alpha = 0.3, color = NA) +
    labs(x = "Collectivity", y = "Percentage", color = "Type", fill = "Type") +
    scale_color_manual(values = c("cascade" = "black", "attenuation" = "blue", 
                                  "amplification" = "darkgreen", "inversion" = "red")) +
    scale_fill_manual(values = c("cascade" = "black", "attenuation" = "blue", 
                                  "amplification" = "darkgreen", "inversion" = "red")) +
    theme(axis.title.x = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title.y = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          legend.text = element_text(size = 12, family = "LM Roman 10"),
          panel.background = element_rect(fill = 'white', color = 'black'),
          legend.background = element_rect(fill = 'white', color = 'white'),
          legend.position = "top",
          legend.box.background = element_rect(fill = 'white', color = 'white'),
          legend.key = element_rect(fill = "transparent")) +
    guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))  # Enlever le titre de la légende
  
  
  # Figure Each Chain
  EachChain <- ggplot(outputs, aes(x = Collect, y = Value)) +
    geom_point(aes(colour = Type), alpha = 0.2, size = 2) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_color_manual(name = "", values = c("attenuated" = "blue",
                                                 "amplificated" = "darkgreen",
                                                 "inverted" = "red",
                                                 "classic cascade" = "black")) +
    scale_y_continuous(trans = pseudolog10_trans, breaks = c(-10^2, -10^1, -10^0, 0, 10^0, 10^1, 10^2)) +
    labs(x = "Collectivity", y = "Ratio \n net cascade / n-step cascade") +
    guides(fill = guide_legend(title = NULL)) +
    scale_alpha(guide = "legend", range = c(1, 1)) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    theme(axis.title.x = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          axis.title.y = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          legend.text = element_text(size = 12, family = "LM Roman 10"),
          panel.background = element_rect(fill = 'white', color = 'black'),
          # panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
          legend.background = element_rect(fill = 'white', color = 'white'),
          legend.position = "top",
          legend.box.background = element_rect(fill = 'white', color = 'white'),
          legend.key = element_rect(fill = "transparent"))
  
  IntervalCollect <- seq(0, 2, by = 0.1)
  outputs$Intervalle <- cut(outputs$Collect, breaks = IntervalCollect, labels = FALSE)
  BarplotType <- ggplot(outputs, aes(x = as.factor(Intervalle), fill = Type)) +
    geom_bar(position = "fill", stat = "count") +
    scale_x_discrete(labels = paste0("[", IntervalCollect[-length(IntervalCollect)], ", ", IntervalCollect[-1], ")")) +
    scale_y_continuous(labels = scales::percent_format(scale = 100)) +
    scale_fill_manual(values = c("attenuated" = "blue",
                                "amplificated" = "darkgreen",
                                "inverted" = "red",
                                "classic cascade" = "black")) +
    labs(x = "Collectivity", y = "Percentage") +
    guides(fill = guide_legend(title = NULL)) +
    theme(axis.title.x = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.y = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          legend.text = element_text(size = 12, family = "LM Roman 10"),
          panel.background = element_rect(fill = 'white', color = 'black'),
          legend.position = "top",
          legend.key = element_rect(fill = "transparent"))
  
  # scatterplot avec le % de chaînes avec inversion pour chaque réseau de métrique Collectivity
  PercentInv <- ggplot(percent_inv, aes(x = Collect, y = InvPercent)) +
    geom_point(size = 2, alpha = 0.2) +
    geom_smooth(method = "loess", color = "black", fill = "red") +
    labs(x = "Collectivity", y = "Percentage of trophic chains \n with long-term inversion") +
    theme(axis.title.x = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          axis.title.y = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          legend.text = element_text(size = 12, family = "LM Roman 10"),
          panel.background = element_rect(fill = 'white', color = 'black'),
          panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
          legend.key = element_rect(fill = "transparent"))
  
  # scatterplot avec le % de chaînes avec inversion pour chaque réseau de métrique Omnivory
  PercentInvOmni <- ggplot(percent_inv, aes(x = Omni, y = InvPercent)) +
    geom_point(size = 2, alpha = 0.2) +
    geom_smooth(method = "loess", color = "black", fill = "red") +
    labs(x = "Omnivory", y = "Percentage of trophic chains \n with long-term inversion") +
    theme(axis.title.x = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          axis.title.y = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          legend.text = element_text(size = 12, family = "LM Roman 10"),
          panel.background = element_rect(fill = 'white', color = 'black'),
          panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
          legend.key = element_rect(fill = "transparent"))
  
  # Figure Full Foodweb
  colnames(outputs_b) <- c("Collect", "Omni", "InvTot", "InvPercent", "Type", "RatioTot")
  
  EachFoodweb <- ggplot(outputs_b, aes(x = Collect, y = RatioTot)) +
    geom_point(aes(colour = Type), alpha = 0.2, size = 2) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_color_manual(name = "", values = c("attenuated" = "blue",
                                             "amplificated" = "darkgreen",
                                             "inverted" = "red",
                                             "classic cascade" = "black")) +
    scale_y_continuous(trans = pseudolog10_trans, breaks = c(-10^2, -10^1, -10^0, 0, 10^0, 10^1, 10^2)) +
    labs(x = "Collectivity", y = "Ratio \n net cascade / n-step cascade") +
    guides(fill = guide_legend(title = NULL)) +
    scale_alpha(guide = "legend", range = c(1, 1)) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    theme(axis.title.x = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          axis.title.y = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          legend.text = element_text(size = 12, family = "LM Roman 10"),
          panel.background = element_rect(fill = 'white', color = 'black'),
          # panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
          legend.background = element_rect(fill = 'white', color = 'white'),
          legend.box.background = element_rect(fill = 'white', color = 'white'),
          legend.position = "top",
          legend.key = element_rect(fill = "transparent"))
  
  EachFoodwebOmni <- ggplot(outputs_b, aes(x = Omni, y = RatioTot)) +
    geom_point(aes(colour = Type), alpha = 0.2, size = 2) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_color_manual(name = "", values = c("attenuated" = "blue",
                                             "amplificated" = "darkgreen",
                                             "inverted" = "red",
                                             "classic cascade" = "black")) +
    scale_y_continuous(trans = pseudolog10_trans, breaks = c(-10^2, -10^1, -10^0, 0, 10^0, 10^1, 10^2)) +
    labs(x = "Omnivory", y = "Ratio \n net cascade/n-step causality cascade") +
    guides(fill = guide_legend(title = NULL)) +
    scale_alpha(guide = "legend", range = c(1, 1)) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    theme(axis.title.x = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          axis.title.y = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          legend.text = element_text(size = 12, family = "LM Roman 10"),
          panel.background = element_rect(fill = 'white', color = 'black'),
          # panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
          legend.background = element_rect(fill = 'white', color = 'white'),
          legend.box.background = element_rect(fill = 'white', color = 'white'),
          legend.key = element_rect(fill = "transparent"))
  
  IntervalCollect <- seq(0, 2.1, by = 0.1)
  outputs_b$Intervalle <- cut(outputs_b$Collect, breaks = IntervalCollect, labels = FALSE)
  BarplotTypeFoodweb <- ggplot(outputs_b, aes(x = as.factor(Intervalle), fill = Type)) +
    geom_bar(position = "fill", stat = "count") +
    scale_x_discrete(labels = paste0("[", IntervalCollect[-length(IntervalCollect)], ", ", IntervalCollect[-1], ")")) +
    scale_y_continuous(labels = scales::percent_format(scale = 100)) +
    scale_fill_manual(values = c("attenuated" = "blue",
                                 "amplificated" = "darkgreen",
                                 "inverted" = "red",
                                 "classic cascade" = "black")) +
    labs(x = "Collectivity", y = "Percentage") +
    guides(fill = guide_legend(title = NULL)) +
    theme(axis.title.x = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.y = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          legend.text = element_text(size = 12, family = "LM Roman 10"),
          panel.background = element_rect(fill = 'white', color = 'black'),
          legend.position = "top",
          legend.key = element_rect(fill = "transparent"))
  
  # Percent of full web inversion in function of intervals of Collectivity
  # breaks <- c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5)
  breaks <- c(0, 0.5, 1, 1.5, 2)
  # breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2)
  # breaks <- seq(0, 2, 0.1)
  Aggregated_data <- outputs_b %>% mutate(interval = cut(Collect,
                          breaks, 
                          include.lowest = TRUE, 
                          right = FALSE)) %>%
    group_by(interval) %>% 
    summarize(Avg_Inv = mean(InvTot, na.rm = TRUE)*100,
              Avg_InvPercent = mean(InvPercent, na.rm = TRUE),
              Avg_Omni = mean(Omni, na.rm = TRUE)) %>%
    ungroup()
  
  PercentFullInv <- ggplot(Aggregated_data, aes(x = interval, y = Avg_Inv, 
                          color = Avg_Omni, size = as.factor(round(Avg_InvPercent)))) +
    geom_point() +
    scale_colour_gradient(low = "blue", high = "red", name = "Omnivory") +
    scale_size_discrete(name = "% of chains \n with inversion") +
    labs(x = "Collectivity", y = "% of food-webs with full inversion", title = "") +
    theme(axis.title.x = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          axis.title.y = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          legend.title = element_text(size = 12, face = "bold", family = "LM Roman 10"),
          panel.background = element_rect(fill = 'white', color = 'black'),
          panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
          legend.key = element_rect(fill = "transparent"))
  
  BarplotFullInv <- ggplot(outputs_b, aes(x = factor(Collect), fill = factor(InvTot))) +
    geom_bar(stat = "count", position = "stack") +
    labs(x = "Collectivity", y = "Food webs distribution") +
    scale_x_discrete(breaks = sort(unique(outputs_b$Collect))[seq(1, length(unique(outputs_b$Collect)), by = 2)],
                     labels = sort(unique(outputs_b$Collect))[seq(1, length(unique(outputs_b$Collect)), by = 2)]) +
    scale_fill_manual(values = c("blue", "red"), labels = c("No full inversion", "Full inversion")) +
    guides(fill = guide_legend(title = NULL)) +
    theme(axis.title.x = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          axis.title.y = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          legend.text = element_text(size = 12, family = "LM Roman 10"),
          panel.background = element_rect(fill = 'white', color = 'black'),
          legend.key = element_rect(fill = "transparent"))
  
  # Percent of full web inversion in function of intervals of Omnivory
  breaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
  Aggregated_data <- outputs_b %>% mutate(interval = cut(Omni,
                                                         breaks, 
                                                         include.lowest = TRUE, 
                                                         right = FALSE)) %>%
    group_by(interval) %>% 
    summarize(Avg_Inv = mean(InvTot, na.rm = TRUE)*100,
              Avg_InvPercent = mean(InvPercent, na.rm = TRUE)) %>%
    ungroup()
  
  PercentFullInvOmni <- ggplot(Aggregated_data, aes(x = interval, y = Avg_Inv, color = Avg_InvPercent)) +
    geom_point(size = 3) +
    scale_colour_gradient(low = "blue", high = "red", name = "% of chains \n with inversion") +
    labs(x = "Omnivory", y = "% of food-webs with full inversion", title = "") +
    theme(axis.title.x = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          axis.title.y = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          legend.title = element_text(size = 12, face = "bold", family = "LM Roman 10"),
          panel.background = element_rect(fill = 'white', color = 'black'),
          panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
          legend.key = element_rect(fill = "transparent"))
  
  BarplotFullInvOmni <- ggplot(outputs_b, aes(x = factor(Omni), fill = factor(InvTot))) +
    geom_bar(stat = "count", position = "stack") +
    labs(x = "Omnivory", y = "Food webs distribution") +
    scale_x_discrete(breaks = sort(unique(outputs_b$Omni))[seq(1, length(unique(outputs_b$Omni)), by = 2)],
                     labels = sort(unique(outputs_b$Omni))[seq(1, length(unique(outputs_b$Omni)), by = 2)]) +
    scale_fill_manual(values = c("blue", "red"), labels = c("No full inversion", "Full inversion")) +
    guides(fill = guide_legend(title = NULL)) +
    theme(axis.title.x = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          axis.title.y = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          legend.text = element_text(size = 12, family = "LM Roman 10"),
          panel.background = element_rect(fill = 'white', color = 'black'),
          legend.key = element_rect(fill = "transparent"))
  

  out <- list(EachChain, BarplotType, PercentFullInv, BarplotFullInv, 
              PercentFullInvOmni, BarplotFullInvOmni, 
              PercentInv, PercentInvOmni, EachFoodweb, EachFoodwebOmni,
              BarplotTypeFoodweb, PerFwbMean)
  return(out)
}
