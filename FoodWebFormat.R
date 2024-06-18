FoodWebFormat <- function(FoodWebs, ConvCoeff = 0.5) {
  
  # extract names of foodwebs
  NamesWebs <- distinct(FoodWebs, network_name)
  AllFwbs <- list()
  for (k in 1:nrow(NamesWebs)){
    # Create a single list of all species
    OneWeb <- subset(FoodWebs, network_name == NamesWebs$network_name[[k]])
    OneWeb$connection_strength <- as.numeric(OneWeb$connection_strength)
    AllSpecies <- unique(c(OneWeb$species1, OneWeb$species2))
    # Initialize a square matrix with zeros
    OneWebMatrix <- matrix(0, nrow = length(AllSpecies), ncol = length(AllSpecies), dimnames = list(AllSpecies, AllSpecies))
    # Fill the adjacency matrix with connection values
    for (i in 1:nrow(OneWeb)) {
      species1 <- OneWeb[i, "species1"]
      species2 <- OneWeb[i, "species2"]
      connection <- OneWeb[i, "connection_strength"]
      OneWebMatrix[species1, species2] <- connection
    }
    # values in the matrix are only positive, they indicate a positive link between the element
    # in row i to the element in column j, i.e. only the effects of prey on preds.
    # we transpose the matrix so that the values indicate the effect of the element in column j 
    # on the element in row i, to match our way of scoring an interaction matrix A
    OneWebMatrix <- t(OneWebMatrix)
    
    # Normalization in [0, 1]
    Vec <- as.vector(OneWebMatrix)
    OneWebMatrixNorm <- scale(Vec, center = min(Vec), scale = max(Vec) - min(Vec))
    OneWebMatrixNorm <- matrix(OneWebMatrixNorm, nrow = nrow(OneWebMatrix), ncol = ncol(OneWebMatrix), byrow = TRUE)
    
    # we need to create negative antysimetric values, i.e. the effect of pred on prey
    for (i in 1:nrow(OneWebMatrixNorm)){
      for (j in 1:nrow(OneWebMatrixNorm)){
        if (i != j){
          if (OneWebMatrixNorm[i, j] == 0 & OneWebMatrixNorm[j, i] != 0){
            OneWebMatrixNorm[i, j] <- - OneWebMatrixNorm[j, i] / ConvCoeff
          }
        }
      }
    }

    # Checking associated cells
    # for (i in 1:nrow(OneWebMatrix)) {
    #   for (j in 1:nrow(OneWebMatrix)) {
    #     if (OneWebMatrix[i, j] != 0 && OneWebMatrix[j, i] == 0) {
    #       # Display indices of non-conforming cells
    #       cat("Cell (", i, ",", j, ") is non-zero, but cell (", j, ",", i, ") is zero.", "\n")
    #     }
    #   }
    # }
    
    # check if max trophic level is >= 3
    flow <- -1 * OneWebMatrixNorm * (OneWebMatrixNorm < 0)
    flow[is.na(flow)] <- 0
    Troph <- TrophInd(flow)
    Trophic <- round(Troph$TL)
    MaxTroph <- max(Trophic)
    if (MaxTroph >= 3){
      AllFwbs[[k]] <- OneWebMatrixNorm
    }    

  }
  return(AllFwbs)
}
