#' expand
#'
#' This function is used to expand the Y(category) to a indicator vector
#'
#' @param pri the prior state
#' @param curr the current state
#' @param I a U by U incidence matrix with elements; I(i,j)=1 if state j can be accessed from state i in one step and 0 otherwise
#' @param refE a vector with the reference categories
#'
#' @return ry: a indicator vector
#'
#'
#'
#'
#'
#'

expand=function(pri,curr,I,refE){    #expand it into indicate variable
  ruler=I[pri,]  #get the corresponding transition for the given prior state
  ry=rep(0,sum(ruler)-1)   #get a zero vector with length of transitions-1 for the given prior state
  refc=refE[pri]            #get the reference category for this transition
  ri=which(ruler!=0)      #get the response location
  ci=match(curr,ri)          #get the location for the current state
  if(curr!=refc){
    if(curr<refc){ry[ci]=1}else{ry[ci-1]=1}
  }
  return(ry)
}

