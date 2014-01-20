##this file contains helper functions; these are not user accessible
##any small repeated tasks needed for other functions should be put here. 
##contains: isBinary, correctFormat, isCorrectLength, charEdgelist, splittify, getPresentCases


#checks whether a given column is binary or not
#allows NA values in binary variables
#returns T/F or an error if it can't treat the object as a matrix
isBinary <- function(arg) {
  binary <- TRUE
  arg <- as.matrix(na.exclude(arg))
  if (is.matrix(arg)) {
    u_mat <- unique(arg)
    for (each in u_mat) { 
      if (each != 0 && each != 1) { 
        binary <- FALSE; break
        } 
      } #breaks and returns FALSE as soon as it hits a non-missing, non (0,1) value
    return(binary) 
  } 
  else { 
    print('Blau dimension columns must be coercable to matrix')
  }
}


#takes an argument and converts to numeric by looking up the column in the square input dataset
#arguments can be length 1 or longer
correctFormat <- function(arg, square.data) {
  if (!is.null(arg)){ #should be a max length of 1
    if (length(arg) <= 1) {
      if (is.character(arg)) { 
        arg_temp <- type.convert(arg, as.is=T)
        arg_temp <- which(colnames(square.data) == arg_temp) ; return(arg_temp)
      }
      else if (is.numeric(arg)) { 
        return(arg)
      }    
      else {
        sprintf('Input for option %s could not be treated as a character or vector', arg)
      } 
    }
    else if (length(arg > 1)){
      argList <- rep(0,length(arg))
      if (is.numeric(arg)) { 
        argList <- arg ; return(unique(argList))
      } 
      else {
        for (each in 1:length(arg)) {
          arg_temp <- type.convert(arg[each], as.is=T)
          if (is.numeric(arg_temp)) {
            argList[each] <- as.numeric(arg_temp)
          } else if (is.character(arg_temp)) {
            argList[each] <- as.numeric(which(colnames(square.data) == arg_temp))
          } else print('errors')
        } 
        return(unique(argList))
      }
    } 
  }
}


#allows specification of correct length for an input element
#for example, the primaryMembershipCol should be only length 1
#This checks the user's input, sees if it's length 1 (default), and returns an error if it isn't
isCorrectLength <- function(arg, argLength = 1) { #spits out an error if argument is of wrong length. returns true if there are no errors (NULL or correct length)
  if (!is.null(arg)){
    if (length(arg) != argLength) {sprintf('Option %s takes argument of at most length %s', arg, argLength)}
    else {TRUE}
  }
  else {TRUE}
}


#takes a numeric edgelist and list of names (like the kind produced with the network package)
#returns a character matrix that's a named edgelist
#assumes names are in order/correspond to numbers of edges
charEdgelist <- function(edgelist, vertexNames) {
  char <- matrix('0', ncol = ncol(edgelist), nrow = nrow(edgelist))

  for (col in 1:ncol(edgelist)){
    for (row in 1:nrow(edgelist)) {
      char[row, col] <- vertexNames[edgelist[row, col]]
    }
  }
  return(char)
}


splittify <- function(blauObj, ecologyId, ecologyRows) {
  miniBlau <- list() #this object holds the ROW parts of the blau object that are in the relevant ecology
  miniBlau$ids <- blauObj$ids[ecologyRows, , drop=FALSE]
  miniBlau$memberships <- blauObj$memberships[ecologyRows, , drop=FALSE]
  miniBlau$weights <- blauObj$weights[ecologyRows, , drop=FALSE]
  miniBlau$dimensions <- blauObj$dimensions[ecologyRows, , drop=FALSE]

  if (!is.null(blauObj$primMemCol)){
    miniBlau$primMemCol <- blauObj$primMemCol
  }
  if (!is.null(blauObj$isInNiche)){
    miniBlau$isInNiche <- blauObj$isInNiche[which(blauObj$isInNiche[,'ecologyNames'] == ecologyId), which(colnames(blauObj$isInNiche) != 'ecologyNames'), drop=FALSE] #this picks the rows in the ecology, and all columns except for the one with the ecology names
  }

  if (!is.null(blauObj$topbounds) && !is.null(blauObj$lowbounds)){
    miniBlau$topbounds <- blauObj$topbounds[which(blauObj$topbounds[,'ecologyNames'] == ecologyId), which(colnames(blauObj$topbounds) != 'ecologyNames'), drop=FALSE] 
    miniBlau$lowbounds <- blauObj$lowbounds[which(blauObj$lowbounds[,'ecologyNames'] == ecologyId), which(colnames(blauObj$lowbounds) != 'ecologyNames'), drop=FALSE]     
  }

  #keep the full network because the networked function will delete extraneous vertices
  miniBlau$graph <- blauObj$graph

  return(miniBlau)
}


#easily cuts down all active elements of object to a specified set of rows
#used for fiddling with how missing cases are handled
getPresentCases <- function(blauObj, presentCases){
    blauObj$ids <- blauObj$ids[presentCases, , drop=FALSE]
    blauObj$memberships <- blauObj$memberships[presentCases, , drop=FALSE]
    blauObj$dimensions <- blauObj$dimensions[presentCases, , drop=FALSE]
    blauObj$weights <- blauObj$weights[presentCases, , drop=FALSE]

    if (!is.null(blauObj$isInNiche)) {
      blauObj$isInNiche <- blauObj$isInNiche[presentCases, , drop=FALSE]
    }
    
    if (!is.null(blauObj$primaryMembership)) {
      blauObj$primaryMembership <- blauObj$primaryMembership[presentCases, , drop=FALSE]
    }

    return(blauObj)

    #we don't cut down connections, we keep the full network object
    #we don't cut down top/lowbounds because they're computed in a manner specified with the input function and aren't the same size
    #we don't cut down 'final' elements such as nodalLocal because they've been computed already according to user specification
}