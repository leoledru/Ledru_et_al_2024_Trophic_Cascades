FigEachChainDirectNet_c <- function(A_list, Collect_list, Troph_list, Omni_list){
  ### Version qui travaille directement sur des réseaux déjà générés et stables ###
  # VERSION POUR RÉSEAUX AVEC MaxTL = 4, ET CASCADE D'ORDRE 3 DU SOMMET À LA BASE 
  
  outputs <- data.frame(Collect = numeric(0), Omni = numeric(0), Type = character(0), Value = numeric(0))
  outputs_b <- data.frame(Collect = numeric(0), Omni = numeric(0), InvTot = numeric(0), PercentInv = numeric(0))
  percent_inv <- data.frame(Collect = numeric(0), Omni = numeric(0), InvPercent = numeric(0))
  for (i in 1:length(A_list)){
    # For each chain with cascade (positive order two interaction) compute direct and net cascades
    Direct <- A_list[[i]]
    diag(Direct) <- 0
    A3 <- Direct%^%3
    Trophic <- Troph_list[[i]]
    # Compute Net effects
    I <- diag(1, nrow = nrow(Direct)) # matrice identité
    Net <- solve(I - Direct) # -(-I + Direct)^-1 := (I - A)^-1
    NbrInv <- list()
    EachDirect <- list()
    EachNet <- list()
    MaxTroph <- max(Trophic)
    IdxMaxTroph <- which(Trophic == MaxTroph)
    Minus2Troph <- which(Trophic == (MaxTroph - 3))
    for (top in IdxMaxTroph){
      for (bottom in Minus2Troph){
        if ((A3[bottom, top] + Direct[bottom, top]) < 0){ # if trophic cascade (negative order 3 effect)
          DirectCascade <- A3[bottom, top] + Direct[bottom, top]
          if (DirectCascade > 0) {
            print("directe négative, anormal")
          }
          NetCascade <- Net[bottom, top]
          # ONLY NET, WIEGHTED BY DIRECT
          outputs <- rbind(outputs, list(Collect_list[[i]], Omni_list[[i]], "relative", 
                                         NetCascade/DirectCascade))
          
          ifelse(sign(NetCascade) != sign(DirectCascade),
                 NbrInv <- append(NbrInv, 1),
                 NbrInv <- append(NbrInv, 0))
          EachDirect <- append(EachDirect, DirectCascade)
          EachNet <- append(EachNet, NetCascade)
        }
      }
    }
    NbrInv <- unlist(NbrInv)
    percent_inv <- rbind(percent_inv, list(Collect_list[[i]], Omni_list[[i]], sum(NbrInv)/length(NbrInv)*100))
    outputs_b <- rbind(outputs_b, list(round(Collect_list[[i]],1), round(Omni_list[[i]],1), 
                                       ifelse(sign(sum(unlist(EachDirect))) != sign(sum(unlist(EachNet))), 1, 0),
                                       sum(NbrInv)/length(NbrInv)*100))
  }
  colnames(outputs) <- c("Collect", "Omni", "Type", "Value")
  colnames(outputs_b) <- c("Collect", "Omni", "InvTot", "InvPercent")
  colnames(percent_inv) <- c("Collect", "Omni", "InvPercent")
  
  EachChain <- ggplot(outputs, aes(x = Collect, y = Value, 
                                   colour = ifelse(Value < 0, "red", "blue"))) +
    geom_point(alpha = 0.5, size = 2) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_y_continuous(trans = pseudolog10_trans, breaks = c(-10^2, -10^1, -10^0, 0, 10^0, 10^1, 10^2)) +
    labs(x = "Collectivity", y = "Weighted net cascade by second-order cascade") +
    theme(axis.title.x = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          axis.title.y = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          legend.text = element_text(size = 12, family = "LM Roman 10"),
          panel.background = element_rect(fill = 'white', color = 'black'),
          # panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
          legend.background = element_rect(fill = 'white', color = 'white'),
          legend.box.background = element_rect(fill = 'white', color = 'white'),
          legend.key = element_rect(fill = "transparent")) +
    # Configuration de la légende
    scale_color_identity(guide = "legend", 
                         labels = c("blue" = "Positive", "red" = "Negative"),
                         breaks = c("blue", "red"),
                         name = "")
  
  # Percent of full web inversion in function of intervals of Collectivity
  # breaks <- c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5)
  # breaks <- c(0, 0.5, 1, 1.5, 2)
  # breaks <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.4, 1.6, 1.8, 2)
  breaks <- seq(0.5, 7.5, 0.5)
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
  
  # Figure façon scatterplot avec le % de chaînes avec inversion pour chaque réseau de métrique Collectivity
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
  
  # Figure façon scatterplot avec le % de chaînes avec inversion pour chaque réseau de métrique Omnivory
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
  
  out <- list(EachChain, PercentFullInv, BarplotFullInv, 
              PercentFullInvOmni, BarplotFullInvOmni, 
              PercentInv, PercentInvOmni)
  return(out)
}
