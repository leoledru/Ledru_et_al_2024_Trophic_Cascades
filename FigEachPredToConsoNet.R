FigEachPredToConsoNet <- function(A_list, Collect_list, Troph_list){
  #' @title Effect of each predator : Figure A.5
  #' @description
    #' Measure the net mean effect of each predator on its consumer & the net effect of each predator on *each* of its consumer ; to showing the double compensation in function of collectivity
  #' @param A_list is a list of which each element is an interaction matrix
  #' @param Collect_list is a list containing the value of collectivity of each interaction matrix of A_list
  #' @param Troph_list is a list of which each element is a list containing the trophic level of each species of the interaction matrix corresponding in A_list
  #' @returns the two figures in a list
  
  Outputs <- data.frame(Collect = numeric(0), PredToConsoWeigthed = numeric(0))
  OutputsEach <- data.frame(Collect = numeric(0), EachPredEachConso = numeric(0))

  for (i in 1:length(A_list)){
    Direct <- A_list[[i]]
    diag(Direct) <- 0
    Trophic <- Troph_list[[i]]
    # Compute Net effects
    I <- diag(1, nrow = nrow(Direct))
    Net <- solve(I - Direct) # -(-I + Direct)^-1 := (I - A)^-1
    MaxTroph <- max(Trophic)
    IdxMaxTroph <- which(Trophic == MaxTroph)
    
    # For each Pred find its Preys and summed net effects (weighted by number of preys, or not)
    for (top in IdxMaxTroph){
      IdxConsos <- which(Direct[top, ] > 0)
      PredToConsoWeighted <- sum(Net[IdxConsos, top]) / length(IdxConsos)
      Outputs <- rbind(Outputs, list(Collect_list[[i]], PredToConsoWeighted))
      for (conso in IdxConsos){
        OutputsEach <- rbind(OutputsEach, list(Collect_list[[i]], Net[conso, top]))
      }
    }
  }
  colnames(Outputs) <- c("Collect", "PredToConsoWeighted")
  colnames(OutputsEach) <- c("Collect", "EachPredEachConso")
  
  PredToConsoWeighted <- ggplot(Outputs, aes(x = Collect, y = PredToConsoWeighted)) +
    geom_point(alpha = 0.1, size = 2) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
    geom_smooth(color = "black", fill = "red") +
    labs(x = "Collectivity", y = "Net summed effects of each predator on its preys\n(weigthed by the number of preys)") +
    theme(axis.title.x = element_text(size = 20, face = "bold", family = "LM Roman 10"),
          axis.title.y = element_text(size = 20, face = "bold", family = "LM Roman 10"),
          legend.text = element_text(size = 12, family = "LM Roman 10"),
          panel.background = element_rect(fill = 'white', color = 'black'),
          # panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
          legend.background = element_rect(fill = 'white', color = 'white'),
          legend.box.background = element_rect(fill = 'white', color = 'white'),
          legend.key = element_rect(fill = "transparent"))
  
  EachPredEachConso <- ggplot(OutputsEach, aes(x = Collect, y = EachPredEachConso)) +
    geom_point(alpha = 0.05, size = 2) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
    geom_smooth(color = "red", fill = "red", se = FALSE) +
    labs(x = "Collectivity", y = "Net effects of each predator on each of its preys") +
    theme(axis.title.x = element_text(size = 20, face = "bold", family = "LM Roman 10"),
          axis.title.y = element_text(size = 20, face = "bold", family = "LM Roman 10"),
          legend.text = element_text(size = 12, family = "LM Roman 10"),
          panel.background = element_rect(fill = 'white', color = 'black'),
          # panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
          legend.background = element_rect(fill = 'white', color = 'white'),
          legend.box.background = element_rect(fill = 'white', color = 'white'),
          legend.key = element_rect(fill = "transparent"))
  
  return(list(PredToConsoWeighted, EachPredEachConso))
}