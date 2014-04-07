##a place for doing blauspace computations
##these functions are NOT user callable
##there are two parts to this: computation functions and splitter functions
##computation functions compute measures for one ecology
##splitter functions split data up by ecologies and pass them to the computation functions one at a time.

##COMPUTATION FUNCTIONS
##these functions actually compute blauspace measures
##they should all be set up to calculate for one ecology


#calc.niches computes the topbounds, lowbounds and isInNiche for a single ecology
#if an ecology has no people in given membership, NAs will be returned for those bounds, and zeroes will be returned for niche membership
#if an ecology has exactly one person in a membership, the bounds are a point (topbounds=lowbounds)
#this uses all available data, even if some observations have NAs in other dimensions than the bound being computed
calc.niches <- function(blauObj, dev.range) {
  
  #initialize data objects
  topbounds <- matrix(0, ncol = ncol(blauObj$dimensions), nrow = ncol(blauObj$memberships))
  lowbounds <- matrix(0, ncol = ncol(blauObj$dimensions), nrow = ncol(blauObj$memberships)) 
  means <- matrix(0, ncol = ncol(blauObj$dimensions), nrow = ncol(blauObj$memberships))

  #calculate top and low boundaries
  for (memCyc in 1:ncol(blauObj$memberships)) {
    for (dimCyc in 1:ncol(blauObj$dimensions)) {

      memRows <- which(blauObj$memberships[,memCyc] == 1)
      dimRows <- blauObj$dimensions[memRows, dimCyc]
      memRows <- memRows[!is.na(dimRows)] #gets rid of the missing values in the relevant dimension
      meanData <- blauObj$dimensions[memRows, dimCyc] #rows for relevant dimension

      if (length(meanData) == 0){ #for when there is no information
        #this can happen in one of two cases:
        #1) no members in the group
        #2) all members of group have NA along the relevant dimension
        means[memCyc, dimCyc] <- NA
        topbounds[memCyc, dimCyc] <- NA
        lowbounds[memCyc, dimCyc] <- NA
      }

      else if (length(meanData) == 1){
        #impute our only information if there's 1 obs for the dimension
        means[memCyc,dimCyc] <- meanData #should be just a number
        topbounds[memCyc,dimCyc] <- meanData
        lowbounds[memCyc,dimCyc] <- meanData
      }

      else if (length(meanData) > 1) {
        meanWeights <- blauObj$weights[memRows,]
        means[memCyc,dimCyc] <- sum(meanData*meanWeights)/sum(meanWeights) 
        # Calculate the standard deviation
        # Information on weighted Standard Deviation found at
        # http://www.sosmath.com/CBB/viewtopic.php?t=2656
        sdDenominator <- ((length(meanWeights) - 1) * sum(meanWeights))/(length(meanWeights))
        sdNumerator <- 0
        for (dataCyc in 1:length(meanData)){
          sdNumerator <- sdNumerator + (meanWeights[dataCyc] * (meanData[dataCyc] - means[memCyc,dimCyc])^2 )
        }
        stdDev <- sqrt(sdNumerator/sdDenominator)
        topbounds[memCyc, dimCyc] <- means[memCyc, dimCyc] + stdDev * dev.range
        lowbounds[memCyc, dimCyc] <- means[memCyc, dimCyc] - stdDev * dev.range
      }
    }
  }
  blauObj$topbounds <- topbounds
  blauObj$lowbounds <- lowbounds
  
  colnames(blauObj$topbounds) <- colnames(blauObj$dimensions)
  rownames(blauObj$topbounds) <- colnames(blauObj$memberships)
  colnames(blauObj$lowbounds) <- colnames(blauObj$dimensions)
  rownames(blauObj$lowbounds) <- colnames(blauObj$memberships)
  

  #calculate if each node is in a given niche
  blauObj$isInNiche <- matrix(0, nrow = nrow(blauObj$memberships), ncol = ncol(blauObj$memberships))

  #the inside 'apply' takes each row in dimensions and checks if it's within the boundaries
  #the outside 'apply' checks if all elements of each row in the matrix are true
  for (memCyc in 1:nrow(blauObj$lowbounds)){
    blauObj$isInNiche[,memCyc] <- apply(t(apply(blauObj$dimensions, 1, function(x) x >= blauObj$lowbounds[memCyc,] & x <= blauObj$topbounds[memCyc,])), 1, all)
  }

  #overwrite NAs with zeroes
  blauObj$isInNiche[is.na(blauObj$isInNiche)] <- 0

  colnames(blauObj$isInNiche) <- colnames(blauObj$memberships)
  rownames(blauObj$isInNiche) <- rownames(blauObj$memberships)
  
  return(blauObj)
}


#takes modes: local, global
#calculates nodal measures not requiring network information
calc.nodal <- function(blauObj, mode){ 

  #requires focalNiche (primMem) specification
  if (mode == 'local'){

    #initialize
    blauObj$nodalLocal <- as.data.frame(matrix(0, nrow = nrow(blauObj$memberships), ncol= 3))

    rownames(blauObj$nodalLocal) <- rownames(blauObj$isInNiche)

    #in focal niche
    blauObj$nodalLocal[,1] <- blauObj$isInNiche[, blauObj$primMemCol]

    #total number of niches individual is in
    blauObj$nodalLocal[,2] <- matrix(apply(blauObj$isInNiche, 1, function(x) sum(x, na.rm=TRUE)), ncol = 1, byrow = TRUE)

    #if individual is in primary org but outside of primary niche
    for(nodeCyc in 1:length(blauObj$isInNiche[, blauObj$primMemCol])){
      if (!is.na(blauObj$isInNiche[nodeCyc, blauObj$primMemCol]) && !is.na(blauObj$memberships[nodeCyc, blauObj$primMemCol])){
        if (blauObj$isInNiche[nodeCyc, blauObj$primMemCol] == 0 && blauObj$memberships[nodeCyc, blauObj$primMemCol] == 1){
          blauObj$nodalLocal[nodeCyc, 3] <- 1
        }
      }
    }
    return(blauObj)
  }

  #does not require focalNiche(primMem)
  else if (mode == 'global'){

    #number of organizations individual is in
    orgs <- matrix(apply(blauObj$memberships, 1, function(x) sum(x, na.rm = TRUE)), ncol = 1, byrow = TRUE)

    #number of niches individual is in
    niches <- matrix(apply(blauObj$isInNiche, 1, function(x) c(sum(x, na.rm = TRUE), c(paste(which(x == 1), collapse=' ')))), ncol = 2, byrow = TRUE)

    blauObj$nodalGlobal <- cbind(orgs, niches)

    rownames(blauObj$nodalGlobal) <- rownames(blauObj$isInNiche)

    return(blauObj)
  }
}


#calculates nodal spanners
#depends on the network package
calc.nodal.network <- function(blauObj){

  #initialize
  blauObj$nodalNetwork <- as.data.frame(matrix(0, nrow = nrow(blauObj$memberships), ncol= 2))
  rownames(blauObj$nodalNetwork) <- rownames(blauObj$isInNiche)

  #gets rid of nodes not in the current ecology
  namelist <- network.vertex.names(blauObj$graph)
  diff_names <- setdiff(namelist, rownames(blauObj$dimensions))
  blauObj$graph <- delete.vertices(blauObj$graph, vapply(diff_names, function(x) which(namelist == x), 1))

  edgelist <- as.matrix(blauObj$graph, matrix.type='edgelist')

  #make a named edgelist, makes our computations easier
  charEL <- charEdgelist(edgelist, attr(edgelist, 'vnames'))

  #if we're given an undirected graph (undirected EL/symmetric adjacency matrix)
  #duplicate the EL with the origin nodes reversed
  if (is.directed(blauObj$graph) == FALSE) {
    charEL <- rbind(charEL, cbind(charEL[,2], charEL[,1]))
  }

  #sort edgelist by first element
  if (nrow(charEL) > 1){
    charEL <- charEL[order(charEL[, 1]), ]
  }

  #this is kind of a confusing piece of code at first
  #it sets a 'current' origin node and cycles through all of that node's neighbors
  #when it hits a new 'current' node, it records all of the information for the previous 'current' node
  #then it resets the list of niches spanned to and begins recording information on the new current node
  currentNode <- charEL[1,1]
  spannedTo <- c()

  #cycle through directed edgelist
  #the origin node is element 1, the destination node is element 2
  for (rowCyc in 1:nrow(charEL)){
    edge <- as.vector(charEL[rowCyc,])

    #since EL is sorted, if we see a different origin node,
    #record changes to nodalNetwork
    #update current node
    #reset spannedTo
    if (edge[1] != currentNode){
      blauObj$nodalNetwork[currentNode,1] <- ifelse(length(spannedTo) > 0, 1, 0)
      blauObj$nodalNetwork[currentNode,2] <- length(spannedTo)

      #start new spanner record
      currentNode <- edge[1]
      spannedTo <- c()
      niches1 <- blauObj$isInNiche[edge[1], ]
      niches2 <- blauObj$isInNiche[edge[2], ]
      spannedTo <- union(spannedTo, (which((niches2 - niches1) == 1)))
    }

    else {
      niches1 <- blauObj$isInNiche[edge[1], ]
      niches2 <- blauObj$isInNiche[edge[2], ] 

      #nodal spanners are defined as:
      #node1 is not in nicheA but has a friend in nicheA
      #node1 is then said to 'span' to nicheA

      #niches spanned to are indicated by 1's
      #we get number spanned to
      spannedTo <- union(spannedTo, (which((niches2 - niches1) == 1)))
    }
  }
  
  #save the last elements when loop stops
  blauObj$nodalNetwork[currentNode,1] <- ifelse(length(spannedTo) > 0, 1, 0)
  blauObj$nodalNetwork[currentNode,2] <- length(spannedTo)

  return(blauObj)
}


#calculates dyadic statuses
#depends on the network package
calc.dyadic <- function(blauObj, m.dist) { 

  #gets rid of extraneous nodes
  nameList <- network.vertex.names(blauObj$graph)
  diff_names <- setdiff(nameList, rownames(blauObj$dimensions))
  blauObj$graph <- delete.vertices(blauObj$graph, vapply(diff_names, function(x) which(nameList == x), 1))

  edgelist <- as.matrix(blauObj$graph, matrix.type='edgelist')

  charEL <- charEdgelist(edgelist, attr(edgelist, 'vnames'))

  #if we're given an undirected graph (undirected EL/symmetric adjacency matrix)
  #duplicate the EL with the origin nodes reversed
  if (is.directed(blauObj$graph) == FALSE) {
    charEL <- rbind(charEL, cbind(charEL[,2], charEL[,1]))
  }

  #sort edgelist by first element
  if (nrow(charEL) > 1){
    charEL <- charEL[order(charEL[, 1]), ]
  }

  if (m.dist == TRUE){
    blauObj$dyadic <- as.data.frame(matrix(0, nrow = nrow(charEL), ncol = 6))
  }
  else{
    blauObj$dyadic <- as.data.frame(matrix(0, nrow = nrow(charEL), ncol = 5))
  }

  edgelistNames <- matrix(0, nrow = 0, ncol = 2)

  #here's where we take advantage of treating the network as directed
  for (rowCyc in 1:nrow(charEL)){
    edge <- as.vector(charEL[rowCyc,])

    edgelistNames <- rbind(edgelistNames, c(edge[1], edge[2]))

    for (niche1 in 1:ncol(blauObj$isInNiche)){

      #straddler
      #ego in niche
      #alter outside of all niches
      if (blauObj$isInNiche[edge[1], niche1] == 1 && sum(blauObj$isInNiche[edge[2],]) == 0){

        blauObj$dyadic[rowCyc, 3] <- blauObj$dyadic[rowCyc, 3] + 1
      }

      for (niche2 in 1:ncol(blauObj$isInNiche)){
        if (niche1 == niche2){
          #co nicher
          if ((blauObj$isInNiche[edge[1], niche1] == 1) && (blauObj$isInNiche[edge[2], niche2] == 1)){
            blauObj$dyadic[rowCyc, 1] <- blauObj$dyadic[rowCyc, 1] + 1
          }
        }
        else if (niche1 != niche2){
          #spanner
          #ego in niche 1, not in niche 2
          #alter in niche 2 (may also be in niche 1)
          #or vice versa with ego/alter
          if (blauObj$isInNiche[edge[1], niche1] == 1 && blauObj$isInNiche[edge[1], niche2] != 1 && blauObj$isInNiche[edge[2], niche2] == 1){
            blauObj$dyadic[rowCyc, 4] <- blauObj$dyadic[rowCyc, 4] + 1
          }
          if (blauObj$isInNiche[edge[2], niche1] == 1 && blauObj$isInNiche[edge[2], niche2] != 1 && blauObj$isInNiche[edge[1], niche2] == 1){
            blauObj$dyadic[rowCyc, 4] <- blauObj$dyadic[rowCyc, 4] + 1
          }

          #alter in a niche
          #ego outside of all niches
          if (sum(blauObj$isInNiche[edge[1],]) == 0 && blauObj$isInNiche[edge[2], niche2] == 1){

            blauObj$dyadic[rowCyc, 3] <- blauObj$dyadic[rowCyc, 3] + 1
          }
          if (sum(blauObj$isInNiche[edge[2],]) == 0 && blauObj$isInNiche[edge[1], niche2] == 1){

            blauObj$dyadic[rowCyc, 3] <- blauObj$dyadic[rowCyc, 3] + 1
          }
        }
      }
    }

    #co-outsider
    if (sum(blauObj$isInNiche[edge[1],]) + sum(blauObj$isInNiche[edge[2],]) == 0 ){
      blauObj$dyadic[rowCyc, 2] <- 1
    }

    #euclidean dist
    blauObj$dyadic[rowCyc,5] <- dist(rbind(blauObj$dimensions[edge[1],], blauObj$dimensions[edge[2],]), method='euclidean')

    #mahalanobis dist
    if (m.dist == TRUE){
      blauObj$dyadic[rowCyc,6] <- sqrt(mahalanobis(blauObj$dimensions[edge[1],], blauObj$dimensions[edge[2],], cov(blauObj$dimensions)))
    }
  }

  blauObj$dyadic <- cbind(edgelistNames, blauObj$dyadic)

  return(blauObj)
}



##SPLITTER FUNCTIONS
##these functions calculate measures for lots of ecologies by splitting up a dataset and calling the basic computation functions on each ecology individually
##use the 'splittify' function to do this


#calculates bounds and isInNiche for each ecology across each dimension
calc.niches.ecology <- function(blauObj, uniqueEcologies, dev.range){
  
  uniqueEcologies <- unique(blauObj$ids[,2])
  uniqueEcologies <- uniqueEcologies[!is.na(uniqueEcologies)]
  
  blauObj$isInNiche <- matrix(0, nrow = nrow(blauObj$dimensions), ncol = (ncol(blauObj$memberships) + 1)) #extra column for ecology names
  colnames(blauObj$isInNiche) <- c(colnames(blauObj$memberships), 'ecologyNames')
  rownames(blauObj$isInNiche) <- blauObj$ids[,1]
  
  for(ecologyId in uniqueEcologies){ #iterate through each ecology: all of the calculations for the ecology happen here and they are appended to $isInNiche, $topbounds, and $lowbounds
    ecologyRows <- which(blauObj$ids[,2] == ecologyId) #pull out ROW identifiers for each row in the ecology
    
    miniBlau <- splittify(blauObj, ecologyId, ecologyRows)

    miniBlau <- calc.niches(miniBlau, dev.range) #memberships, dimensions, primaryMemberships are used by niches
    
    blauObj$isInNiche[ecologyRows,] <- cbind(miniBlau$isInNiche, (rep(ecologyId, nrow(miniBlau$isInNiche))))
    
    topbounds <- cbind(as.data.frame(miniBlau$topbounds), as.data.frame(rep(ecologyId, nrow(miniBlau$topbounds)))) #temp object to add ecology names 
    lowbounds <- cbind(as.data.frame(miniBlau$lowbounds), as.data.frame(rep(ecologyId, nrow(miniBlau$lowbounds)))) #temp object to add ecology names
    
    colnames(topbounds) <- c(colnames(blauObj$dimensions), 'ecologyNames')
    rownames(topbounds) <- colnames(blauObj$memberships)
    colnames(lowbounds) <- c(colnames(blauObj$dimensions), 'ecologyNames')
    rownames(lowbounds) <- colnames(blauObj$memberships)
    
    blauObj$topbounds <- rbind(blauObj$topbounds, topbounds) #add it to the bottom
    blauObj$lowbounds <- rbind(blauObj$lowbounds, lowbounds) #add it to the bottom
    
  }

  tempData <- blauObj$isInNiche[,which(colnames(blauObj$isInNiche) != 'ecologyNames')]
  class(tempData) <- 'numeric'

  blauObj$isInNiche <- cbind(as.data.frame(tempData), blauObj$isInNiche[,which(colnames(blauObj$isInNiche) == 'ecologyNames')])

  colnames(blauObj$isInNiche)[ncol(blauObj$isInNiche)] <- 'ecologyNames'

  return(blauObj)

}


#calculates the nodal measures for each ecology
#does different things if the 'mode' argument is different
calc.nodal.ecology <- function(blauObj, uniqueEcologies, mode){

  uniqueEcologies <- unique(blauObj$ids[,2])
  uniqueEcologies <- uniqueEcologies[!is.na(uniqueEcologies)]
  
  if (mode == 'local'){
    blauObj$nodalLocal <- matrix(0, nrow = 0, ncol = 2)

    for(ecologyId in uniqueEcologies){

      ecologyRows <- which(blauObj$ids[,2] == ecologyId)
      miniBlau <- splittify(blauObj, ecologyId, ecologyRows)
      
      miniBlau <- calc.nodal(miniBlau, mode) #memberships, weights, dimensions, primaryMemberships are used by niches
      blauObj$nodalLocal <- rbind(blauObj$nodalLocal, miniBlau$nodalLocal)
    }
  }

  else if (mode == 'global'){
    blauObj$nodalGlobal <- matrix(0, nrow = 0, ncol = 3)

    for(ecologyId in uniqueEcologies){

      ecologyRows <- which(blauObj$ids[,2] == ecologyId)
      miniBlau <- splittify(blauObj, ecologyId, ecologyRows)
      
      miniBlau <- calc.nodal(miniBlau, mode) #memberships, weights, dimensions, primaryMemberships are used by niches
      blauObj$nodalGlobal <- rbind(blauObj$nodalGlobal, miniBlau$nodalGlobal)
    }
  }

  if (mode == 'network'){
    blauObj$nodalNetwork <- matrix(0, nrow = 0, ncol = 2)

    for(ecologyId in uniqueEcologies){

      ecologyRows <- which(blauObj$ids[,2] == ecologyId)
      miniBlau <- splittify(blauObj, ecologyId, ecologyRows)
      
      miniBlau <- calc.nodal.network(miniBlau) #memberships, weights, dimensions, primaryMemberships are used by niches
      blauObj$nodalNetwork <- rbind(blauObj$nodalNetwork, miniBlau$nodalNetwork)
    }
  }
  return(blauObj)
}


#computes dyadic statuses for ecologies
calc.dyadic.ecology <- function(blauObj, m.dist){ #splitter function
  if (m.dist == TRUE){
    blauObj$dyadic <- as.data.frame(matrix(0, nrow = 0, ncol= 8))
  }
  else{
    blauObj$dyadic <- as.data.frame(matrix(0, nrow = 0, ncol= 7))
  }

  uniqueEcologies <- unique(blauObj$ids[,2])
  uniqueEcologies <- uniqueEcologies[!is.na(uniqueEcologies)]

  for(ecologyId in uniqueEcologies) {
    ecologyRows <- which(blauObj$ids[,2] == ecologyId)
    
    miniBlau <- splittify(blauObj, ecologyId, ecologyRows)

    miniBlau <- calc.dyadic(miniBlau, m.dist)

    blauObj$dyadic <- rbind(blauObj$dyadic, miniBlau$dyadic)

  }
  return(blauObj)
}