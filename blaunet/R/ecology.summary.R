ecology.summary <-
function(blauObj, percent = FALSE){
  uniqueEcologies <- unique(blauObj$ids[,2])
  uniqueEcologies <- uniqueEcologies[!is.na(uniqueEcologies)]

  ecologySum <- matrix(0, nrow = 0, ncol = (ncol(blauObj$memberships) + 2)) #ecology name, niche name, num in org, num in niche, num in org but not niche, num not in any other niche, overlap 2+ other niches ##no boundaries, too cluttered

  colnames(ecologySum) <- c("Ecology", "Org/Niche", colnames(blauObj$memberships))

  for (ecologyId in uniqueEcologies){

    ecologyRows <- which(blauObj$ids[,2] == ecologyId)

    focalNiches <- blauObj$isInNiche[ecologyRows, 1:(ncol(blauObj$isInNiche)-1)]

    sum.mat <- t(as.matrix(focalNiches)) %*% as.matrix(focalNiches)

    orig.diag <- diag(sum.mat)

    #how many individuals are exclusively in a niche?
    #get this by summing up all isInNiche rows with sum = 1
    mat.diagonal <- rep(0, ncol(focalNiches))

    for (node in 1:nrow(focalNiches)){
      if(sum(focalNiches[node,]) == 1){
        mat.diagonal <- mat.diagonal + focalNiches[node,]
      }
    }

    #manually replace because diag function gets confused
    for (elem in 1:ncol(mat.diagonal)){
      sum.mat[elem,elem] <- mat.diagonal[1,elem]
    }


    if (percent == TRUE){
      sum.mat <- sum.mat / orig.diag
    }

    ecologySum <- rbind(ecologySum, cbind(cbind(rep(ecologyId, ncol(blauObj$memberships)), colnames(blauObj$memberships)), sum.mat))
  }
  rownames(ecologySum) <- NULL
  return(ecologySum)
}
