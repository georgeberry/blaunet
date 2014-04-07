##place for user-callable functions
##this is how the user interfaces with the underlying functions that actually do type checking, partitioning, and calculations
##repeated tasks here should ALL be done by a function in calculations.r or helper_functions.r
##functions here generally have two cases: one ecology or many ecologies; this pattern should be followed throughout

##it's important to note that functions here take and return a blau object
##functions (like overlaps) that require the output of previous functions should call those functions if the object is not present already
##this preserves continuity and minimizes the potential for user error
##it also avoids tedious sequential calls


#this function is for when the user only needs to compute the niches
#ecologies.off will automatically treat the dataset as having one ecology, even if ecology names are specified
niches <- function(blauObj, dev.range = 1.5, ecologies.off = FALSE){
  uniqueEcologies <-  unique(blauObj$ids[,2])

  if (length(uniqueEcologies) == 1 || ecologies.off == TRUE){
    blauObj <- calc.niches(blauObj, dev.range)
    rownames(blauObj$isInNiche) <- rownames(blauObj$memberships)
  }
  else if (length(uniqueEcologies)) {
    blauObj <- calc.niches.ecology(blauObj, uniqueEcologies, dev.range)
  }

  presentCases <- which(complete.cases(blauObj$dimensions))
  
  blauObj <- getPresentCases(blauObj, presentCases)

  return(blauObj)
}


#computes nodal statuses requiring a focalNiche
nodal.local <- function(blauObj, focal.niche=NULL, dev.range=1.5, ecologies.off=FALSE){
  
  if (is.null(focal.niche)){
    return("Primary Membership needed for nodal.local")
  }

  blauObj$primMemCol <- correctFormat(focal.niche, blauObj$memberships)

  if (ecologies.off == TRUE){
    blauObj <- niches(blauObj, dev.range, ecologies.off)
  }

  uniqueEcologies <- unique(blauObj$ids[,2])
  if(length(uniqueEcologies) == 1 || ecologies.off == TRUE){ #if there's only one ecology and, we don't have isInNiche

    if(is.null(blauObj$isInNiche)){
      blauObj <- niches(blauObj, dev.range)
    }

    blauObj <- calc.nodal(blauObj, mode = "local") #has isInNiches now; does a bunch of stuff if a primaryMembership is specified
  }

  else if(length(uniqueEcologies) > 1) { #if there's more than one ecology, we need to split up primary membership, weights, dimensions and memberships
    if(is.null(blauObj$isInNiche)){
      blauObj <- niches(blauObj, dev.range)
    }    
    blauObj <- calc.nodal.ecology(blauObj, uniqueEcologies, mode = "local")
  }

  colnames(blauObj$nodalLocal) <- c("FocNicher", "Nicher", "MemNotNiche")
  
  return(blauObj)

}


#computes nodal statudes not requiring primary membership
nodal.global <- function(blauObj, dev.range=1.5, ecologies.off = FALSE){
  
  if (ecologies.off == TRUE){
    blauObj <- niches(blauObj, dev.range, ecologies.off)
  }

  uniqueEcologies <- unique(blauObj$ids[,2])
  if(length(uniqueEcologies) == 1 || ecologies.off == TRUE){ #if there's only one ecology and ,we don't have isInNiche

    if(is.null(blauObj$isInNiche)){
      blauObj <- niches(blauObj, dev.range)
    }

    blauObj <- calc.nodal(blauObj, mode = "global") #has isInNiches now; does a bunch of stuff if a primaryMembership is specified
  }

  else if(length(uniqueEcologies) > 1) { #if there's more than one ecology, we need to split up primary membership, weights, dimensions and memberships
    if(is.null(blauObj$isInNiche)){
      blauObj <- niches(blauObj, dev.range)
    }    
    blauObj <- calc.nodal.ecology(blauObj, uniqueEcologies, mode = "global")
  }
  colnames(blauObj$nodalGlobal) <- c("TotalOrgs", "Nicher", "NicheList")
  
  return(blauObj)
}


#computes nodal statudes requiring network information
nodal.network <- function(blauObj, dev.range = 1.5, ecologies.off = FALSE){

  if (ecologies.off == TRUE){
    blauObj <- niches(blauObj, dev.range, ecologies.off)
  }

  uniqueEcologies <- unique(blauObj$ids[,2])
  if(length(uniqueEcologies) == 1 || ecologies.off == TRUE){ #if there's only one ecology and ,we don't have isInNiche

    if(is.null(blauObj$isInNiche)){
      blauObj <- niches(blauObj, dev.range)
    }

    blauObj <- calc.nodal.network(blauObj) #has isInNiches now; does a bunch of stuff if a primaryMembership is specified
  }

  else if(length(uniqueEcologies) > 1) { #if there's more than one ecology, we need to split up primary membership, weights, dimensions and memberships
    if(is.null(blauObj$isInNiche)){
      blauObj <- niches(blauObj, dev.range)
    }    
    blauObj <- calc.nodal.ecology(blauObj, uniqueEcologies, mode = "network")
  }
  colnames(blauObj$nodalNetwork) <- c("Spanner", "NumSpannedTo")
  
  return(blauObj)

}

#this function computes dyadic measures such as our 4 types of spanners and euclidean/mahalanobis distance
#output takes the form of an edgelist with the computed qualities
dyadic <- function(blauObj, dev.range = 1.5, ecologies.off=FALSE, m.dist = TRUE) {

  if (ecologies.off == TRUE){
    blauObj <- niches(blauObj, dev.range, ecologies.off)
  }
  
  uniqueEcologies <- unique(blauObj$ids[,2])
  if(length(uniqueEcologies) == 1 || ecologies.off == TRUE){ #if there's only one ecology and ,we don't have isInNiche

    if(is.null(blauObj$isInNiche)){
      blauObj <- niches(blauObj, dev.range)
    }

    blauObj <- calc.dyadic(blauObj, m.dist)

  }

  else if(length(uniqueEcologies) > 1){
    if(is.null(blauObj$IsInNiche)){
      blauObj <- niches(blauObj, dev.range)
    }
    blauObj <- calc.dyadic.ecology(blauObj, m.dist)
  }

  if (m.dist == TRUE){
    colnames(blauObj$dyadic) <- c('Ego', 'Alter', 'CoNicher', 'CoOutsider', 'Straddler', 'Spanner', 'EucDist', 'MahalanobisDist')
  }
  else{
    colnames(blauObj$dyadic) <- c('Ego', 'Alter', 'CoNicher', 'CoOutsider', 'Straddler', 'Spanner', 'EucDist')
  }

  return(blauObj)
}

#will return the nodal things blauNet has computed
#when adding new working functions, this can be extended easily
export.nodal <- function(blauObj, niches = TRUE){
  if (is.null(blauObj$isInNiche)){
    print("Nothing to export.")
  }

  if (niches == TRUE){
    if ("ecologyNames" %in% colnames(blauObj$isInNiche)){
      to.export <- cbind(blauObj$ids,blauObj$isInNiche[, 1:(ncol(blauObj$isInNiche)-1)])
    }
    else{
      to.export <- cbind(blauObj$ids,blauObj$isInNiche)
    }
  }
  else{
    to.export <- data.frame(matrix(0, nrow = nrow(blauObj$nodalLocal), ncol=0))
  }

  if (!is.null(blauObj$nodalLocal)){
    to.export <- cbind(to.export, blauObj$nodalLocal)
  }
  if (!is.null(blauObj$nodalGlobal)){
    to.export <- cbind(to.export, blauObj$nodalGlobal)
  }
  if (!is.null(blauObj$nodalNetwork)){
    to.export <- cbind(to.export, blauObj$nodalNetwork)
  }
  rownames(to.export) <- NULL
  return(to.export)
}


export.dyadic <- function(blauObj){
  if (is.null(blauObj$dyadic)){
    print("Nothing to export.")
  }
  else{
    return(blauObj$dyadic)
  }
}


#function counts up the number of members in each membership, nichers in each niche, and total people in each ecology
#it returns a matrix with ecologies as columns and these measures as rows
niche.summary <- function(blauObj){

  uniqueEcologies <- unique(blauObj$ids[,2])
  uniqueEcologies <- uniqueEcologies[!is.na(uniqueEcologies)]

  sum.niche <- matrix(0, nrow = length(uniqueEcologies)*ncol(blauObj$memberships), ncol = 7) #ecology name, niche name, num in org, num in niche, num in org but not niche, num not in any other niche, overlap 2+ other niches ##no boundaries, too cluttered

  colnames(sum.niche) <- c("Ecology", "Org/Niche", "OrgMem", "NicheMem", "NicheExc", "NicheOvr", "MemExc")

  rowCount <- 1

  for (ecologyId in uniqueEcologies){

    nicheNum <- 1

    ecologyRows <- which(blauObj$ids[,2] == ecologyId)

    focalMemberships <- blauObj$memberships[ecologyRows, ]
    focalNiches <- blauObj$isInNiche[ecologyRows, 1:(ncol(blauObj$isInNiche)-1)] #exclude last column, which is ecology index

    for (colCyc in 1:ncol(focalMemberships)){
      #in org but not in niche
      #basically find the 1's
      diff <- focalMemberships[, nicheNum] - focalNiches[, nicheNum]
      numOutside <- length(which(diff == 1))

      #in focal niche but in no other niche
      nicheExcl <- sum(apply(focalNiches, 1, function(x) ifelse(sum(x) == 1 && x[nicheNum] == 1, 1, 0)), na.rm = TRUE)

      #overlaps with at least 2 niches
      numNonExclusive <- sum(focalNiches[, nicheNum], na.rm = TRUE) - nicheExcl

      sum.niche[rowCount,] <- c(ecologyId, colnames(blauObj$memberships)[nicheNum], sum(focalMemberships[, nicheNum], na.rm = TRUE), sum(focalNiches[, nicheNum], na.rm = TRUE), nicheExcl, numNonExclusive, numOutside)

      rowCount <- rowCount + 1
      nicheNum <- nicheNum + 1
    }
  }

  return(sum.niche)
}


ecology.summary <- function(blauObj, percent = FALSE){
  uniqueEcologies <- unique(blauObj$ids[,2])
  uniqueEcologies <- uniqueEcologies[!is.na(uniqueEcologies)]

  sum.ecology <- matrix(0, nrow = 0, ncol = (ncol(blauObj$memberships) + 2)) #ecology name, niche name, num in org, num in niche, num in org but not niche, num not in any other niche, overlap 2+ other niches ##no boundaries, too cluttered

  colnames(sum.ecology) <- c("Ecology", "Org/Niche", colnames(blauObj$memberships))

  for (ecologyId in uniqueEcologies){

    ecologyRows <- which(blauObj$ids[,2] == ecologyId)

    focalNiches <- blauObj$isInNiche[ecologyRows, 1:(ncol(blauObj$isInNiche)-1)]

    sum.mat <- t(as.matrix(focalNiches)) %*% as.matrix(focalNiches)

    orig.diag <- diag(sum.mat)

    #how many individuals are exclusively in a niche?
    #get this by summing up all isInNiche rows with sum = 1
    mat.diagonal <- rep(0, ncol(focalNiches))

    if (length(mat.diagonal) > 0){

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

    sum.ecology <- rbind(sum.ecology, cbind(cbind(rep(ecologyId, ncol(blauObj$memberships)), colnames(blauObj$memberships)), sum.mat))
	}
  }
  rownames(sum.ecology) <- NULL
  return(sum.ecology)
}


#very simple: prints active objects so a user can see a list of things to be accessed
active <- function(blauObj){
  sprintf('These elements are active: %s. Access with object$element.', paste(names(blauObj), collapse = ", "))
}