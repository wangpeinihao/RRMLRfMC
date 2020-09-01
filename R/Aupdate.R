#' Aupdate
#'
#' This function is used to update A matrix
#'
#' @param Dfix the coefficient matrix for fixed variables
#' @param Gamma the G matrix value
#' @param Adata the dataset
#' @param iniA the initial A matrix; if missing, a constant matrix with elements 0.01
#' @param R the rank of reduced rank model
#' @param ptrans a vector with transition number for each state
#' @param  eps the tolerance for convergence, default is 10^-5
#'
#' @return a list of outputs:
#' \itemize{
#' \item NewA: the updated A matrix
#' \item loglikeA: the loglikelihood when updating A
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

Aupdate=function(Dfix,Gamma,Adata,iniA,R,ptrans,eps){  #Adata=pri,curr,pred,fpred,obstrans
  pred=as.matrix(Adata$pred)
  p=ncol(pred)
  #y=apply(matrix(curr,n,1), 1, expand)
  #zyo=cbind(t(y),pred,fpred, pri)
  zyo=cbind(Adata$pri,Adata$curr,Adata$pred,Adata$fpred,Adata$obstrans)


  if(missing(iniA)) iniA <- matrix(rep(0.01,p*R),p,R)#iniA <- matrix(rnorm(p*R),p,R) #give an initial value to A iniA <- matrix(runif(p*R,-0.7,0.7),p,R)
  loglikeold=100
  iteAu=0
  deltaAu=10
  Avc.old=iniA
  #while(deltaAu > eps & iteAu <=200){
  while(deltaAu > eps & iteAu <=100){
    iteAu=iteAu+1
    #Dmat=matrix(Dfix,ncol=21)
    Dmat=Dfix
    #Gamma=matrix(Gamma,ncol = 21)
    Gamma=as.matrix(Gamma)

    newton=apply(as.matrix(zyo),1, function(x) derivatives(A=Avc.old,Gamma=Gamma,Dmat=Dmat,ptrans=ptrans,x))
    newton=apply(newton, 1,sum)

    if(sum(is.na(newton))>0){
      #iteAu=1
      deltaAu=10
      #Avc.old=matrix(rep(0.001,p*R),p,R)
      Avc.old=matrix(rep(0,p*R),p,R)
    } else {
      firstdir=newton[1:(p*R)]
      seconddir=matrix(newton[(p*R+1):(p*R*p*R+p*R)],p*R,p*R)
      loglikenew=newton[p*R*p*R+p*R+1]
      Avc.new=as.vector(Avc.old)-solve(seconddir)%*%firstdir
      deltaAu=abs(loglikenew-loglikeold)
      Avc.old=matrix(Avc.new,nrow = p)
      loglikeold=loglikenew
      print(iteAu)
      print(loglikeold)
    }
  }

  Af=matrix(Avc.new,p,R)
  return(list(NewA=Af,loglikeA=loglikeold))
}
