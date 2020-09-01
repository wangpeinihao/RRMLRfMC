#' derivatives
#'
#' This function is used calculate the derivative values (first and second derivatives for Newton-Raphson method) and loglikelihood when updating A
#'
#' @param A matrix with value from previous iteration
#' @param Gamma G matrix values
#' @param Dmat the coefficient matrix for the fixed variables,
#' @param ptrans a vector with transition number for each state
#' @param zy the variable values for a given observation
#'
#' @return a list of outputs:
#' \itemize{
#' \item fird: the first derivative value
#' \item secd: the second derivative value
#' \item loglike: the loglikelihood
#' }
#'
#' @export
#'
#' @examples
#'
#'
#'
#'
#'

derivatives <- function(A,Gamma,Dmat,ptrans,zy){   ##zy=rows of Adata=pri,curr,pred,fpred,obstrans     derivative
  p=nrow(A)
  q=nrow(Dmat)
  zy=unlist(zy)
  pri=zy[1]
  curr=zy[2]
  pstr=zy[length(zy)]                  #get the number of transitions for this obs
  y=matrix(expand(curr,pstr-1),ncol = 1)
  Kd=sum(ptrans-1)
  Gamma=matrix(Gamma,ncol = Kd)
  Dmat=matrix(Dmat,ncol = Kd)
  z=zy[3:(2+p)]
  z=matrix(z,nrow = 1)
  if(q==0){
    WD=0
  }else{
    w=zy[(3+p):(2+p+q)]
    w=matrix(w,nrow = 1)
    WD=w%*%Dmat
  }
  ZAG=z%*%A%*%Gamma
  WaZ=WD+ZAG
  GZ=kronecker(t(Gamma),z) #21 by R*p matrix
  A=matrix(as.vector(A),ncol = 1)
  exppart=exp(WaZ)
  colind=c(0,cumsum(ptrans-1))   #column index to work on

  GZh=GZ[(colind[pri]+1):(colind[pri+1]),,drop=FALSE]
  eph=exppart[(colind[pri]+1):(colind[pri+1])] #exppart for state si
  WaZh=WaZ[(colind[pri]+1):(colind[pri+1])]
  fird=-t(GZh)%*%matrix(eph,ncol=1)/(1+sum(eph))+t(GZh)%*%y
  secd=(t(GZh)%*%matrix(eph,ncol=1)%*%t(t(GZh)%*%matrix(eph,ncol=1))-(t(GZh*eph)%*%GZh*(1+sum(eph))))/(1+sum(eph))^2
  secd=as.vector(secd)
  loglike=-log(1+sum(eph))+ WaZh%*%y

  return(c(fird,secd,loglike))
}

