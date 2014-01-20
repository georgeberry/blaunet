calc.dyadic <-
function(blauObj, m.dist) { 

  #gets rid of extraneous nodes
  nameList <- network.vertex.names(blauObj$graph)
  diff_names <- setdiff(nameList, rownames(blauObj$dimensions))
  blauObj$graph <- delete.vertices(blauObj$graph, vapply(diff_names, function(x) which(nameList == x), 1))

  edgelist <- as.matrix(blauObj$graph, matrix.type="edgelist")

  charEL <- charEdgelist(edgelist, attr(edgelist, "vnames"))

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
