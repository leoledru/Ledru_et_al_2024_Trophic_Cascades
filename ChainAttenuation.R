ChainAttenuation <- function(Omni = FALSE){
  #' @title Figure 1
  #' @description
    #' Measurement of the classic trophic cascade and the net trophic cascade on a simple trophic chain as a function of the consumer's attack rate on the resource.
  #' @param Omni FALSE or TRUE for a link between pred and resource.
  #' @returns nothing, directly create and show the plot
  
  par(family = "LM Roman 10") # latex font
  n = 3
  AttackRate <- seq(from = 0.1, to = 3, length.out = 100)
  PredAttackRate <- 0.1
  ConvCoeff <- 0.5
  CascadeNet <- vector(mode = "list", length(AttackRate))
  CascadeConcept <- vector(mode = "list", length(AttackRate))
  DiffCascade <- vector(mode = "list", length(AttackRate))
  DiffCascadePercent <- vector(mode = "list", length(AttackRate))
  Horizon <- vector(mode = "list", length(AttackRate))
  NetEffectOnMiddle <- vector(mode = "list", length(AttackRate))
  Collect <- vector(mode = "list", length(AttackRate))
  for (i in 1:length(AttackRate)){
    A <- matrix(0, nrow = n, ncol = n)
    I <- diag(1, nrow = n) # matrice identité
    value <- AttackRate[[i]]
    A[1,2] <- -value
    A[2,3] <- -PredAttackRate
    A[2,1] <- value*ConvCoeff
    A[3,2] <- PredAttackRate*ConvCoeff
    if (Omni == TRUE){
      A[1,3] <- -value / 10
      A[3,1] <- (value / 10)*ConvCoeff
    }
    Collect[[i]] <- spectralRadius(A)
    A2 <- A%^%2
    AInv <- solve(I - A)
    CascadeNet[[i]] <- AInv[1,3]
    CascadeConcept[[i]] <- A[1,3] + A2[1,3]
    DiffCascade[[i]] <- CascadeNet[[i]] - CascadeConcept[[i]]
    DiffCascadePercent[[i]] <- (CascadeNet[[i]] - CascadeConcept[[i]]) / CascadeConcept[[i]] * 100
    NetEffectOnMiddle[[i]] <- AInv[2,3]
  }
  
  plot(unlist(CascadeConcept) ~ AttackRate, pch = 19,
       xlab = "Attack rate a", ylab = "Effect of top-species on basal-species",
       cex.lab = 1.3, ylim = c(-min(unlist(CascadeNet)), max(unlist(CascadeConcept))))
       legend("topleft", legend = c("Classic cascade", "Net cascade"), pch = c(19, 1))
  points(AttackRate, unlist(CascadeNet))
  # plot(Collect, unlist(NetEffectOnMiddle))
  # print(unlist(CascadeConcept))
  # print(unlist(CascadeNet))
}
