NicheModelSimple <- function(S, C, ConvCoeff, SelfReg, MaxTL = 4){
  
  while (TRUE){
    n.i <- sort(runif(S), decreasing = F)
    r.i <- rbeta(S, 1 , 1/(2*C)-1)*n.i
    c.i <- runif(S, r.i/2, n.i)
    a <- matrix(0, nrow = S, ncol = S)
    for(i in 1:S){
      for(j in 1:S){
        if (i != j){
          if(n.i[j] > (c.i[i] - (.5 * r.i[i])) & n.i[j] < (c.i[i] + .5 * r.i[i])){
            a[j, i] <- - 1
            a[i, j] <- 1
          }
        }
      }
    }
    values <- rhalfnorm(length(a[a == 1]), theta = sqrt(pi/1)) # half-normal
    if (length(a[a == 1]) != length(values)){
      browser()
    }
    a[a == 1] <- values * ConvCoeff
    if (length(a[a == -1]) != length(values)){
      browser()
    }
    a[a == -1] <- - values
    # no unconnected species
    abis <- a
    diag(abis) <- 0
    IdxUnconnected <- which(rowSums(abs(abis)) == 0 & colSums(abs(abis)) == 0)
    if (length(IdxUnconnected) > 0){ # remove unconnected species
      n.i <- n.i[-IdxUnconnected]
      r.i <- r.i[-IdxUnconnected]
      c.i <- c.i[-IdxUnconnected]
      a <- a[-IdxUnconnected, -IdxUnconnected]
    }
    if (nrow(a) < 3){
      next
    }
    # no trophic level > 4
    # browser()
    flow <- -1 * a * (a < 0)
    flow[is.na(flow)] <- 0
    Troph <- TrophInd(flow)
    Trophic <- round(Troph$TL)
    cond2 <- max(Trophic) <= MaxTL
    if (cond2 == FALSE){ # remove species with too high trophic level
      IdxOverTP <- which(Trophic > MaxTL) # find species with too high trophic level
      n.i <- n.i[-IdxOverTP]
      r.i <- r.i[-IdxOverTP]
      c.i <- c.i[-IdxOverTP]
      Trophic <- Trophic[-IdxOverTP]
      a <- a[-IdxOverTP, -IdxOverTP]
    }
    if (nrow(a) < 3){
      next
    }
    break
  }
  
  Interaction <- a
  Connectance <- sum(Interaction != 0) / nrow(Interaction)^2
  ### Normalisation by self-reg and computation of Net Effects matrix
  diag(Interaction) = rep.int(-SelfReg, nrow(Interaction))
  D <- diag(Interaction) # self-reg matrice
  # Normalisation par la self_reg D pour obtenir dg/dN = D.(-I + A)
  I <- diag(1, nrow = nrow(Interaction)) # matrice identité
  A_norm <- Interaction / -D # A_norm -> (-I + A_ij) avec A_ij = aij/(-aii) (normalisée par self-reg)
  A_ij <- A_norm + I  # A_ij -> non-dimensional direct interactions (avec Aii = 0)
  if (det(I - A_ij ) != 0) { # (-I + A_ij) inversible ssi determinant =! 0 
    InteractionNet <- solve(I - A_ij ) # A_ij_inv -> non-dimensional INdirect interactions -(-I + A_ij)^-1 := (I - A)^-1
  } else {
    print("Matrice A_ij non inversible")
  }
  Collect <- spectralRadius(A_ij)
  InteractionDirect <- A_ij - I
  
  output <- list("direct" = InteractionDirect, "net" = InteractionNet, 
                 "collect" = Collect, "connectance" = Connectance, 
                 "trophic" = Trophic)
  return(output)
}

