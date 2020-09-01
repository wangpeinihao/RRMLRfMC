#' derivativeB
#'
#' This function is used to calculate the loglikelihood with a given matrix B=AG
#'
#' @param B a numeric coefficient matrix
#' @param ptrans a vector with transition number for each state
#' @param zy the variable values for a given observation
#'
#' @return loglikelihood
#'
#' @export
#'
#' @examples
#'
#'
#'
#'
#'

derivativeB <- function(B,ptrans,zy){   ##pri,curr,pred,fpred,obstrans
  p=nrow(B)
  zy=unlist(zy)
  pri=zy[1]
  z=zy[3:(2+p)]   #the p here =p+q
  ZB=z%*%B
  exppart=exp(t(B)%*%as.matrix(z))
  curr=zy[2]
  pstr=zy[length(zy)]
  y=expand(curr,pstr-1)
  colind=c(0,cumsum(ptrans-1))

  ZBh=ZB[(colind[pri]+1):(colind[pri+1])]
  eph=exppart[(colind[pri]+1):(colind[pri+1])] #exppart for state si
  loglike=-log(1+sum(eph))+ ZBh%*%y

  return(loglike)
}
