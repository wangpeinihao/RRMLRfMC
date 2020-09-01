#' Gupdate
#'
#' This function is used to update G matrix
#'
#' @param A numeric matrix
#' @param Gdata the dataset used to update G
#' @param ptrans a vector with transition number for each state
#'
#' @return a list of outputs:
#' \itemize{
#' \item NewG: the updated G matrix
#' \item loglikeK: the loglikelihood when updating G
#' \item sderr: standard errors for the coefficient matrix
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

Gupdate=function(A,Gdata,ptrans){   #Gdata=pri,curr,pred,fpred,obstrans
  pri=Gdata$pri
  si=length(unique(pri))
  A=as.matrix(A)
  pred=as.matrix(unlist(Gdata$pred))
  T=ncol(A)
  if(is.null(Gdata$fpred)){
    fpred=NULL
  }else{
    fpred=as.matrix(unlist(Gdata$fpred))
  }
  q=ncol(fpred)
  if(is.null(q))q=0
  curr=Gdata$curr
  K=sum(ptrans-1)
  G=matrix(0,T+q,K)
  sderr=matrix(0,T+q,K) #store the standard error
  loglikeK=rep(0,si)
  colind=c(0,cumsum(ptrans-1))
  for (i in 1:si){
    resp=as.factor(curr[pri==i])
    predi=cbind(fpred[pri==i,],(pred[pri==i,,drop=FALSE])%*%A)
    rlevi=min(curr[pri==i])
    resp2=relevel(resp,ref=rlevi)
    data=as.data.frame(cbind(resp2,predi))
    fit=multinom(resp2 ~ 0+predi, data = data)
    G[,(colind[i]+1):(colind[i+1])] = t(summary(fit)$coefficients)
    sderr[,(colind[i]+1):(colind[i+1])]=t(summary(fit)$standard.errors)
    loglikeK[i]=-fit$value
  }

  loglikeKt=sum(loglikeK)
  return(list(NewG=G,loglikeK=loglikeKt,sderr=sderr))
}
