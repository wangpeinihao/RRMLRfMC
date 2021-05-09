#' derivativeB
#'
#' This function is used to calculate the loglikelihood with a given matrix B=AG
#'
#' @param B a numeric coefficient matrix
#' @param I U by U incidence matrix with elements; I(i,j)=1 if state j can be
#'     accessed from state i in one step and 0 otherwise
#' @param zy the variable values for a given observation
#' @param refd a vector of reference categories
#'
#' @return loglikelihood
#'
#'
#'
#'
#'
#'

derivativeB <- function(B,I,zy,refd){   ##pri,curr,pred,fpred,obstrans
  p=nrow(B)
  #zy=unlist(zy)
  rsum=apply(I, 1,sum)
  pri=zy[1]
  z=zy[3:(2+p)]   #the p here =p+q
  ZB=z%*%B
  exppart=exp(t(B)%*%as.matrix(z))
  curr=zy[2]
  td=apply(I, 1, sum)
  pstr=td[pri]
  y=expand(pri,curr,I,refE=refd)
  ptrans=rsum[rsum!=0]
  colind=c(0,cumsum(ptrans-1))

  ZBh=ZB[(colind[pri]+1):(colind[pri+1])]
  eph=exppart[(colind[pri]+1):(colind[pri+1])] #exppart for state si
  loglike=-log(1+sum(eph))+ ZBh%*%y

  return(loglike)
}
