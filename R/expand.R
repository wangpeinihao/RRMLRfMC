#' expand
#'
#' This function is used to expand the Y(category) to a indicator vector
#'
#' @param pri the prior state
#' @param curr the current state
#' @param I a U by U incidence matrix with elements; I(i,j)=1 if state j can be accessed from state i in one step and 0 otherwise
#'
#' @return ry: a indicator vector
#'
#' @export
#'
#' @examples
#'
#'
#'
#'
#'

expand=function(pri,curr,I){    #expand it into indicate variable
  ruler=I[pri,]  #get the corresponding transition for the given prior state
  ry=rep(0,sum(ruler)-1)   #get a zero vector with length of transitions-1 for the given prior state
  ri=which(ruler!=0)      #get the response location
  ci=match(curr,ri)          #get the location for the current state
  ry[ci-1]=1                 #give the index 1 to the current state-1 since we assume the 1st cat as reference
  return(ry)
}

