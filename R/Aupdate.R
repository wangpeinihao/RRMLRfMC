#' Aupdate
#'
#' This function is used to update A matrix
#'
#' @param Dfix the coefficient matrix for study covariates
#' @param Gamma the G matrix value
#' @param Adata the dataset
#' @param R the rank of reduced rank model
#' @param p the number of covariates in the dimension reduction
#' @param q the numbne of study covariates
#' @param I a U by U incidence matrix with elements; I(i,j)=1 if state j can be accessed from state i in one step and 0 otherwise
#' @param iniA initial value for the iteration
#' @param eps the tolerance for convergence, default is 10^-5
#' @param refA a vector of reference categories
#'
#' @return a list of outputs:
#' \itemize{
#' \item NewA: the updated A matrix
#' \item loglikeA: the loglikelihood when updating A
#' }
#'
#'
#'
#'
#'
#'

Aupdate=function(Dfix,Gamma, Adata,R,p,q,I,iniA, eps,refA){

   zyo=Adata

  if(missing(iniA)) iniA <- matrix(rep(0.01,p*R),p,R)#iniA <- matrix(rnorm(p*R),p,R) #give an initial value to A iniA <- matrix(runif(p*R,-0.7,0.7),p,R)
  loglikeold=100; iteAu=0; deltaAu=10; Avc.old=iniA
  while(deltaAu > eps & iteAu <=100){
    iteAu=iteAu+1
    Dmat=Dfix
    Gamma=as.matrix(Gamma)

    newton=apply(as.matrix(zyo),1, function(x) derivatives(A=Avc.old,Gamma,Dmat,I,x,refA))
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
