#' Gupdate
#'
#' This function is used to update G matrix
#'
#' @param A numeric matrix
#' @param Gdata the dataset used to update G
#' @param p the number of covariates in the dimension reduction
#' @param q the numbne of study covariates
#' @param I a U by U incidence matrix with elements; I(i,j)=1 if state j can be
#'     accessed from state i in one step and 0 otherwise
#' @param refG a vector of reference categories
#'
#' @return a list of outputs:
#' \itemize{
#' \item NewG: the updated G matrix
#' \item loglikeK: the loglikelihood when updating G
#' \item sderr: standard errors for the coefficient matrix
#' }
#'
#'
#'
#'
#'

Gupdate=function(A,Gdata,p,q,I,refG){   #Gdata=pri,curr,pred,fpred,obstrans
  pri=Gdata[,1]
  si=length(unique(pri))
  A=as.matrix(A)
  pred=Gdata[,(3:(2+p)),drop=FALSE]
  T=ncol(A)
  wt=which(apply(I, 1, sum)!=0)

  if(q==0){
    fpred=NULL
  }else{
    fpred=Gdata[,((3+p):(2+p+q)),drop=FALSE]
  }
  curr=Gdata[,2]
  K=sum(I)-sum(apply(I, 1, sum)!=0)
  G=matrix(0,T+q,K)
  sderr=matrix(0,T+q,K) #store the standard error
  loglikeK=rep(0,si)
  cp=0
  for (i in 1:si){
    ti=wt[i]
    cm=sum(I[ti,])-1
    resp=as.factor(curr[pri==i])
    predi=cbind(fpred[pri==i,],(pred[pri==i,,drop=FALSE])%*%A)
    rlevi=refG[i]
    resp2=stats::relevel(resp,ref=as.character(rlevi))
    data=as.data.frame(cbind(resp2, predi))
    fit=nnet::multinom(resp2 ~ 0+predi, data = data)
    G[,(cp+1):(cp+cm)] = t(summary(fit)$coefficients)
    sderr[,(cp+1):(cp+cm)]=t(summary(fit)$standard.errors)
    loglikeK[i]=-fit$value
    cp=cp+cm
  }

  loglikeKt=sum(loglikeK)
  return(list(NewG=G,loglikeK=loglikeKt,sderr=sderr))
}
