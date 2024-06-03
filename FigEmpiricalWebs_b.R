FigEmpiricalWebs_b <- function(AllFwbs, MaxInter, MinInter, NbrStep){
  # Measurement for every empirical network and for a whole range of self-regulation 
  # from MinInter to MaxInter by NbrStep steps :
  # - The collectivity
  # - The community-cascade inversion/attenuation/amplification/classic percentage
  # - Percentage of species-cascade inversion/attenuation/amplification/classic
  
  outputsInv <- data.frame(Fwb = numeric(0), Variable = numeric(0), Value = numeric(0), Omni = numeric(0))
  outputsAtt <- data.frame(Fwb = numeric(0), Variable = numeric(0), Value = numeric(0), Omni = numeric(0))
  outputsAmp <- data.frame(Fwb = numeric(0), Variable = numeric(0), Value = numeric(0), Omni = numeric(0))
  outputsClass <- data.frame(Fwb = numeric(0), Variable = numeric(0), Value = numeric(0), Omni = numeric(0))
  CollectVsInv <- data.frame(Collect = numeric(0), Inv = numeric(0))
  # Range of Self-reg
  UnlistFbws <- unlist(AllFwbs)
  # MaxInter <- 2*max(UnlistFbws)
  # MinInter <- min(UnlistFbws[UnlistFbws > 0])
  SelfReg <- seq(MinInter, MaxInter, length.out = NbrStep)
  # for each food web
  for (i in 1:length(AllFwbs)){
    A <- AllFwbs[[i]]
    Connectance <- sum(A != 0) / nrow(A)^2
    SC <- nrow(A)*Connectance
    # Trophic levels
    flow <- -1 * A * (A < 0) 
    flow[is.na(flow)] <- 0
    Troph <- TrophInd(flow)
    Trophic <- round(Troph$TL)
    MaxTroph <- max(Trophic)
    IdxMaxTroph <- which(Trophic == MaxTroph)
    Minus2Troph <- which(Trophic == (MaxTroph - 2))
    Omni <- mean(Troph$OI)
    FullInvAll <- c()
    FullAttAll <- c()
    FullAmpAll <- c()
    FullClassAll <- c()
    # for each value of self-regulation
    for (d in SelfReg){
      # no-dimensional Direct and Net matrices
      D <- rep.int(-d, nrow(A))
      diag(A) <- D
      I <- diag(1, nrow = dim(A)[1])
      A_norm <- A / -D # A_norm -> (-I + A_ij)
      InteractionDirect <- A_norm + I
      A2 <- InteractionDirect%^%2
      InteractionNet <- solve(I - InteractionDirect)
      Collect <- spectralRadius(InteractionDirect)
      # Cascades
      EachDirect <- list()
      EachNet <- list()
      NbrInv <- list()
      NbrAtt <- list()
      NbrAmp <- list()
      NbrClass <- list()
      for (top in IdxMaxTroph){
        for (bottom in Minus2Troph){
          if ((A2[bottom, top] + InteractionDirect[bottom, top]) > 0){ # if trophic cascade (positive order 2 effect, despite possible omnivory)
            DirectCascade <- A2[bottom, top] + InteractionDirect[bottom, top]
            NetCascade <- InteractionNet[bottom, top]
            # count of inversion
            ifelse(sign(NetCascade) != sign(DirectCascade),
                   NbrInv <- append(NbrInv, 1),
                   NbrInv <- append(NbrInv, 0))
            # count of attenuation
            ifelse(sign(NetCascade) == sign(DirectCascade) & NetCascade/DirectCascade<0.8,
                   NbrAtt <- append(NbrAtt, 1),
                   NbrAtt <- append(NbrAtt, 0))
            # count of amplification
            ifelse(sign(NetCascade) == sign(DirectCascade) & NetCascade/DirectCascade>1.2,
                   NbrAmp <- append(NbrAmp, 1),
                   NbrAmp <- append(NbrAmp, 0))   
            # count of classic cascade
            ifelse(sign(NetCascade) == sign(DirectCascade) & NetCascade/DirectCascade>0.8 & NetCascade/DirectCascade<1.2,
                   NbrClass <- append(NbrClass, 1),
                   NbrClass <- append(NbrClass, 0))   
            
            EachDirect <- append(EachDirect, DirectCascade)
            EachNet <- append(EachNet, NetCascade)
          }
        }
      }
      NbrInv <- unlist(NbrInv)
      NbrAtt <- unlist(NbrAtt)
      NbrAmp <- unlist(NbrAmp)
      NbrClass <- unlist(NbrClass)
      PercentChainsInv <- sum(NbrInv)/length(NbrInv)*100
      PercentChainsAtt <- sum(NbrAtt)/length(NbrAtt)*100
      PercentChainsAmp <- sum(NbrAmp)/length(NbrAmp)*100
      PercentChainsClass <- sum(NbrClass)/length(NbrClass)*100
      FullInv <- ifelse(sign(sum(unlist(EachDirect))) != sign(sum(unlist(EachNet))), 100, 0)
      FullInvAll <- c(FullInvAll, FullInv)
      FullAtt <- ifelse(sign(sum(unlist(EachDirect))) == sign(sum(unlist(EachNet))) &
                          sum(unlist(EachNet))/sum(unlist(EachDirect)) < 0.8, 100, 0)
      FullAttAll <- c(FullAttAll, FullAtt)
      FullAmp <- ifelse(sign(sum(unlist(EachDirect))) == sign(sum(unlist(EachNet))) &
                          sum(unlist(EachNet))/sum(unlist(EachDirect)) > 1.2, 100, 0)
      FullAmpAll <- c(FullAmpAll, FullAmp)
      FullClass <- ifelse(sign(sum(unlist(EachDirect))) == sign(sum(unlist(EachNet))) &
                          sum(unlist(EachNet))/sum(unlist(EachDirect)) > 0.8 &
                          sum(unlist(EachNet))/sum(unlist(EachDirect)) < 1.2, 100, 0)
      FullClassAll <- c(FullClassAll, FullClass)
      outputsInv <- rbind(outputsInv, list(i, "Collect", log10(Collect), Omni))
      outputsInv <- rbind(outputsInv, list(i, "PercentChainsInv", PercentChainsInv, Omni))
      outputsAtt <- rbind(outputsAtt, list(i, "PercentChainsAtt", PercentChainsAtt, Omni))
      outputsAmp <- rbind(outputsAmp, list(i, "PercentChainsAmp", PercentChainsAmp, Omni))
      outputsClass <- rbind(outputsClass, list(i, "PercentChainsClass", PercentChainsClass, Omni))
      CollectVsInv <- rbind(CollectVsInv, list(Collect, PercentChainsInv))
    }
    outputsInv <- rbind(outputsInv, list(i, "FullInv", mean(FullInvAll), Omni))
    outputsAtt <- rbind(outputsAtt, list(i, "FullAtt", mean(FullAttAll), Omni))
    outputsAmp <- rbind(outputsAmp, list(i, "FullAmp", mean(FullAmpAll), Omni))
    outputsClass <- rbind(outputsClass, list(i, "FullClass", mean(FullClassAll), Omni))
  }
  colnames(outputsInv) <- c("Fwb", "Variable", "Value", "Omni")
  colnames(outputsAtt) <- c("Fwb", "Variable", "Value", "Omni")
  colnames(outputsAmp) <- c("Fwb", "Variable", "Value", "Omni")
  colnames(outputsClass) <- c("Fwb", "Variable", "Value", "Omni")
  colnames(CollectVsInv) <- c("Collect", "Inv")
  
  # Scatter plot collectivity vs percent of species-cascade inversion
  CollectVsInvPlot <- ggplot(CollectVsInv, aes(x = log10(Collect), y = log10(Inv))) +
    geom_point() +
    stat_poly_line() + stat_poly_eq() +
    labs(x = "Collectivity (log scale)", y = "Percentage of species-cascade inversion (log scale)", color = "Omni") +
    theme(axis.title.x = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          axis.title.y = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          panel.background = element_rect(fill = 'white', color = 'black'),
          panel.grid.major = element_line(color = 'black', linetype = 'dotted'))
  
  # Boxplots (Figure published)
  
  ## INVERSION
  FacetNamesInv <- list('Collect'="Collectivity (log scale)",
                     'PercentChainsInv'="Percentage of community-cascade inversion",
                      'FullInv'="Percentage of species-cascade inversion")
  FacetLabeller <- function(variable, value){
    return(FacetNamesInv[value])}
  GplotInv <- ggplot(outputsInv, aes(x = as.factor(Fwb), y = Value, fill = Variable)) +
    geom_boxplot(data = subset(outputsInv, Variable %in% c("Collect", "PercentChainsInv")), position = position_dodge(width = 0.75)) +
    geom_point(data = subset(outputsInv, Variable == "FullInv"),
               position = position_dodge(width = 0.75),
               aes(group = Variable),
               shape = 17,
               size = 3) +
    facet_wrap(~Variable, scales = "free_y", ncol = 1, labeller=labeller(Variable=FacetLabeller)) +
    labs(x = "Food webs", y = "Value") +
    guides(fill = "none") +
    theme(axis.title.x = element_text(size = 14, face = "bold", family = "LM Roman 10"),
                  axis.title.y = element_text(size = 14, face = "bold", family = "LM Roman 10"),
                  legend.text = element_text(size = 12, family = "LM Roman 10"),
                  strip.text = element_text(size = 12, family = "LM Roman 10"),
                  panel.background = element_rect(fill = 'white', color = 'black'),
                  panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
                  legend.key = element_rect(fill = "transparent"))
  
  ## AMPLIFICATION
  FacetNamesAmp <- list('PercentChainsAmp'="Percentage of community-cascade amplification",
                     'FullAmp'="Percentage of species-cascade amplification")
  FacetLabeller <- function(variable, value){
    return(FacetNamesAmp[value])}
  GplotAmp <- ggplot(outputsAmp, aes(x = as.factor(Fwb), y = Value, fill = Variable)) +
    geom_boxplot(data = subset(outputsAmp, Variable == "PercentChainsAmp"), 
                 position = position_dodge(width = 0.75), fill = "#0077BB") +
    geom_point(data = subset(outputsAmp, Variable == "FullAmp"),
               position = position_dodge(width = 0.75),
               aes(group = Variable),
               shape = 17,
               size = 3) +
    facet_wrap(~Variable, scales = "free_y", ncol = 1, labeller=labeller(Variable=FacetLabeller)) +
    labs(x = "Food webs", y = "Value") +
    guides(fill = "none") +
    theme(axis.title.x = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          axis.title.y = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          legend.text = element_text(size = 12, family = "LM Roman 10"),
          strip.text = element_text(size = 12, family = "LM Roman 10"),
          panel.background = element_rect(fill = 'white', color = 'black'),
          panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
          legend.key = element_rect(fill = "transparent"))
  
  ## ATTENUATION
  FacetNamesAtt <- list('PercentChainsAtt'="Percentage of community-cascade attenuation",
                     'FullAtt'="Percentage of species-cascade attenuation")
  FacetLabeller <- function(variable, value){
    return(FacetNamesAtt[value])}
  GplotAtt <- ggplot(outputsAtt, aes(x = as.factor(Fwb), y = Value, fill = Variable)) +
    geom_boxplot(data = subset(outputsAtt, Variable == "PercentChainsAtt"), 
                 position = position_dodge(width = 0.75), fill = "#0077BB") +
    geom_point(data = subset(outputsAtt, Variable == "FullAtt"),
               position = position_dodge(width = 0.75),
               aes(group = Variable),
               shape = 17,
               size = 3) +
    facet_wrap(~Variable, scales = "free_y", ncol = 1, labeller=labeller(Variable=FacetLabeller)) +
    labs(x = "Food webs", y = "Value") +
    guides(fill = "none") +
    theme(axis.title.x = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          axis.title.y = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          legend.text = element_text(size = 12, family = "LM Roman 10"),
          strip.text = element_text(size = 12, family = "LM Roman 10"),
          panel.background = element_rect(fill = 'white', color = 'black'),
          panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
          legend.key = element_rect(fill = "transparent"))
  
  ## CLASSIC
  FacetNamesClass <- list('PercentChainsClass'="Percentage of classic community-cascade",
                     'FullClass'="Percentage of classic species-cascade")
  FacetLabeller <- function(variable, value){
    return(FacetNamesClass[value])}
  GplotClass <- ggplot(outputsClass, aes(x = as.factor(Fwb), y = Value, fill = Variable)) +
    geom_boxplot(data = subset(outputsClass, Variable == "PercentChainsClass"), 
                 position = position_dodge(width = 0.75), fill = "#0077BB") +
    geom_point(data = subset(outputsClass, Variable == "FullClass"),
               position = position_dodge(width = 0.75),
               aes(group = Variable),
               shape = 17,
               size = 3) +
    facet_wrap(~Variable, scales = "free_y", ncol = 1, labeller=labeller(Variable=FacetLabeller)) +
    labs(x = "Food webs", y = "Value") +
    guides(fill = "none") +
    theme(axis.title.x = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          axis.title.y = element_text(size = 14, face = "bold", family = "LM Roman 10"),
          legend.text = element_text(size = 12, family = "LM Roman 10"),
          strip.text = element_text(size = 12, family = "LM Roman 10"),
          panel.background = element_rect(fill = 'white', color = 'black'),
          panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
          legend.key = element_rect(fill = "transparent"))
  
  
  return(list(CollectVsInvPlot, GplotInv, GplotAmp, GplotAtt, GplotClass))
}
