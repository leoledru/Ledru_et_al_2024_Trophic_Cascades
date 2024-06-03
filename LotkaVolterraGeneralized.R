# Définition de la fonction du modèle Lotka-Volterra généralisé
LotkaVolterraGeneralized <- function(t, N, parameters) {
  A <- parameters$A
  r <- parameters$r
  Troph <- parameters$Troph
  r[Troph > 1] <- - 0.01 # set intrinsic death-rate = 0.01 for non-basal species
  
  dN <- N * r 
  interactions <- A * (N %*% t(N))
  dN <- dN + rowSums(interactions)
  return(list(dN))
}


# LotkaVolterraGeneralized <- function(t, N, parameters) {
#   A <- parameters$A
#   r <- parameters$r
#   Troph <- parameters$Troph
#   
#   dN <- rep(0, length(N))
#   for (i in 1:length(N)) {
#     if (Troph[i] == 1) {
#       dN[i] <- N[i] * r[i]  # Termes de croissance pour les espèces basales
#     }
#     for (j in 1:length(N)) {
#       dN[i] <- dN[i] + N[i] * A[i, j] * N[j]  # Termes d'interaction
#     }
#   }
#   return(list(dN))
# }
