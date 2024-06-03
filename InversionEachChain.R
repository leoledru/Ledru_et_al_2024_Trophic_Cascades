InversionEachChain <- function(A, Omni = T){
  
  SelfReg <- 1
  # non-dimensionalization
  D <- rep.int(-SelfReg, nrow(A))
  diag(A) <- D
  I <- diag(1, nrow = dim(A)[1]) # matrice identité
  A_norm <- A / -D # A_norm -> (-I + A_ij) avec A_ij = aij/(-aii) (normalisée par self-reg)
  A_ij <- A_norm + I # A_ij -> non-dimensional direct interactions (avec Aii = 0)
  InteractionNet <- solve(I - A_ij)
  
  # Computation of Cascades in function of horizon
  # calcul niveau trophique (package NetIndices)
  flow <- -1 * A_ij * (A_ij < 0)
  flow[is.na(flow)] <- 0
  Troph <- TrophInd(flow)
  Trophic <- round(Troph$TL)
  MaxTroph <- max(Trophic)
  IdxMaxTroph <- which(Trophic == MaxTroph)
  Minus2Troph <- which(Trophic == (MaxTroph - 2)) # Pour cascade à n-2
  BasalTroph <- which(Trophic == 1)
  if (Omni == F){ # remove omnivory links in cascade interactions
    A_ij[Minus2Troph, IdxMaxTroph] <- 0
    InteractionNet <- solve(I - A_ij)
  }
  
  # Cascade nette agreg
  CascadeNetAgreg <- sum(InteractionNet[Minus2Troph, IdxMaxTroph])
  # Cascade directe agreg
  A2 <- A_ij%^%2
  CascadeDirectAgreg <- sum(A_ij[Minus2Troph, IdxMaxTroph] + A2[Minus2Troph, IdxMaxTroph])
  
  # Cascade nette
  CascadeNet <- vector("list", length(IdxMaxTroph))
  for (i in 1:length(IdxMaxTroph)){
    for (j in 1:length(Minus2Troph)){
      CascadeNet[[i]] <- append(CascadeNet[[i]], InteractionNet[Minus2Troph[[j]], IdxMaxTroph[[i]]])
    }
  }
  # Cascade directe
  A2 <- A_ij%^%2
  CascadeDirect <- vector("list", length(IdxMaxTroph))
  for (i in 1:length(IdxMaxTroph)){
    for (j in 1:length(Minus2Troph)){
      value <- A_ij[Minus2Troph[[j]], IdxMaxTroph[[i]]] + A2[Minus2Troph[[j]], IdxMaxTroph[[i]]]
      CascadeDirect[[i]] <- append(CascadeDirect[[i]], value)
    }
  }
  
  # Which chain(s) have inversion ?
  Inversion <- vector("list", length(IdxMaxTroph))
  for (i in 1:length(IdxMaxTroph)){
    Net <- CascadeNet[[i]]
    Direct <- CascadeDirect[[i]]
    IdxInv <- which(sign(Net) != sign(Direct) & Net != 0 & Direct != 0)
    if (length(IdxInv) > 0){
      Inversion[[i]] <- list(top = IdxMaxTroph[[i]], bottom = Minus2Troph[IdxInv])  
    }else{
      Inversion[[i]] <- NULL
    }
  }
  
  out <- list(CascadeNet = CascadeNet, CascadeDirect = CascadeDirect, 
              CascadeNetAgreg = CascadeNetAgreg, CascadeDirectAgreg = CascadeDirectAgreg,
              Inversion = Inversion, IdxBottom = Minus2Troph)  
  return(out)
}
