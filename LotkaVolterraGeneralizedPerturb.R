# Définition de la fonction du modèle Lotka-Volterra généralisé
LotkaVolterraGeneralizedPerturb <- function(t, N, parameters) {
  A <- parameters$A
  r <- parameters$r
  Troph <- parameters$Troph
  r_perturb <- parameters$r_perturb
  idx_perturb <- parameters$idx_perturb
  r[Troph > 1] <- - 0.01
  r[idx_perturb] <- r[idx_perturb] + r_perturb
  
  dN <- N * r
  interactions <- A * (N %*% t(N))
  dN <- dN + rowSums(interactions)
  return(list(dN))
}
  
#   dN <- rep(0, length(N))
#   for (i in 1:length(N)) {
#     if (Troph[i] == 1) {
#       dN[i] <- N[i] * r[i]  # Termes de croissance pour les espèces basales
#     }
#     if (i == idx_perturb){
#       dN[i] <- N[i] * r_perturb 
#     }
#     for (j in 1:length(N)) {
#       dN[i] <- dN[i] + N[i] * A[i, j] * N[j]  # Termes d'interaction
#     }
#   }
#   return(list(dN))
# }
