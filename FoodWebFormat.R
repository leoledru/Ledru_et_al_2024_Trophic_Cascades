FoodWebFormat <- function(FoodWebs, ConvCoeff = 0.5) {
  
  # extract names of foodwebs
  NamesWebs <- distinct(FoodWebs, network_name)
  AllFwbs <- list()
  for (k in 1:nrow(NamesWebs)){
    # Supposons que votre objet initial soit appelé "data" avec trois colonnes nommées "species1", "species2", et "connection"
    # Créez une liste unique de toutes les espèces
    OneWeb <- subset(FoodWebs, network_name == NamesWebs$network_name[[k]])
    OneWeb$connection_strength <- as.numeric(OneWeb$connection_strength)
    AllSpecies <- unique(c(OneWeb$species1, OneWeb$species2))
    # Initialisez une matrice carrée avec des zéros
    OneWebMatrix <- matrix(0, nrow = length(AllSpecies), ncol = length(AllSpecies), dimnames = list(AllSpecies, AllSpecies))
    # Remplissez la matrice d'adjacence avec les valeurs de connexion
    for (i in 1:nrow(OneWeb)) {
      species1 <- OneWeb[i, "species1"]
      species2 <- OneWeb[i, "species2"]
      connection <- OneWeb[i, "connection_strength"]
      # Mettez à jour la valeur de connexion dans la matrice d'adjacence
      OneWebMatrix[species1, species2] <- connection
    }
    # les valeurs dans la matrice sont uniquement positives, elles indiquent un lien positif entre l'élément
    # de la ligne i vers l'élément de la colonne j (on le voit quand on compare le plot avec la matrice)
    # c'est-à-dire uniquement les effets des proies sur les preds.
    # on transpose la matrice pour que les valeurs indiquent l'effet de l'élément de la colonne j 
    # sur l'élément de la ligne i pour correspondre à notre manière de noter une matrice d'interaction A
    OneWebMatrix <- t(OneWebMatrix)
    
    # # create igraph object from one network
    # OneWeb <- subset(FoodWebs, network_name == NamesWebs$network_name[[k]])
    # OneWebGraph <- OneWeb %>% select(species1, species2, connection_strength) %>%
    #   graph_from_data_frame(directed = FALSE)
    # # convert igraph object into adjacency matrix
    # OneWebMatrix <- as_adjacency_matrix(OneWebGraph,
    #                                     type = "lower", # because it is directed !
    #                                     attr = "connection_strength",
    #                                     sparse = FALSE)
    # # convert elements into numeric values
    # class(OneWebMatrix) <- "numeric"
    # # remove NA values (previous empty strings)
    # OneWebMatrix[which(is.na(OneWebMatrix) == TRUE)] <- 0
    # 
    # OneWebMatrix <- t(OneWebMatrix)
    
    # Normalization in [0, 1]
    # Convertir la matrice en vecteur
    Vec <- as.vector(OneWebMatrix)
    # Normaliser l'ensemble de la matrice entre 0 et 1
    OneWebMatrixNorm <- scale(Vec, center = min(Vec), scale = max(Vec) - min(Vec))
    # Remettre la matrice dans sa forme originale
    OneWebMatrixNorm <- matrix(OneWebMatrixNorm, nrow = nrow(OneWebMatrix), ncol = ncol(OneWebMatrix), byrow = TRUE)
    
    # il faut maintenant créer les valeurs antysimétriques négatives, c'est-à-dire l'effet des pred sur les proies
    # rq : les réseaux téléchargés sont censés être des foodwebs mais il y a parfois des liens positif-positif...
    # test <- 0
    for (i in 1:nrow(OneWebMatrixNorm)){
      for (j in 1:nrow(OneWebMatrixNorm)){
        if (i != j){
          if (OneWebMatrixNorm[i, j] == 0 & OneWebMatrixNorm[j, i] != 0){
            OneWebMatrixNorm[i, j] <- - OneWebMatrixNorm[j, i] / ConvCoeff
          }
        }
        # print(c(OneWebMatrix[i, j], OneWebMatrix[j, i]))
        # if (OneWebMatrix[i, j] + OneWebMatrix[j, i] != 0){
          # test <- test + 1
        # }
      }
    }

    # print(sum(test)) # number of positive-positive links
    
    # Vérification des cellules associées
    # for (i in 1:nrow(OneWebMatrix)) {
    #   for (j in 1:nrow(OneWebMatrix)) {
    #     if (OneWebMatrix[i, j] != 0 && OneWebMatrix[j, i] == 0) {
    #       # Afficher les indices des cellules non conformes
    #       cat("La cellule (", i, ",", j, ") est non-nulle, mais la cellule (", j, ",", i, ") est nulle.", "\n")
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
