#' derivatives
#'
#' This function is used calculate the derivative values (first and second derivatives for Newton-Raphson method) and loglikelihood when updating A
#'
#' @param A matrix with value from previous iteration
#' @param Gamma G matrix values
#' @param Dmat the coefficient matrix for the fixed variables,
#' @param I a U by U incidence matrix with elements; I(i,j)=1 if state j can be accessed from state i in one step and 0 otherwise
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

derivatives <- function(A,Gamma,Dmat,I,zy){   ##zy=rows of Adata=pri,curr,pred,fpred,obstrans     derivative
  rsum=apply(I, 1,sum)
  ptrans=rsum[rsum!=0]
  p=nrow(A)
  q=nrow(Dmat)
  zy=unlist(zy)
  pri=zy[1]
  curr=zy[2]
  td=apply(I, 1, sum)               #get the transition number for each prior state
  pstr=td[pri]                 #get the number of transitions for this obs
  y=matrix(expand(pri,curr,I),ncol = 1)
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

