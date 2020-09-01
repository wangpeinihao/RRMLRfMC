#' rrmultinom
#'
#' This function is used to fit the reduced rank multinomial logistic regression for markov chain
#'
#' @param data the dataset
#' @param rrfactor the index for variables involved in reduced rank, it could be a numeric vector or a vector with variable name
#' @param  ffactor the index for fixed variables, it could be a numeric vector or a vector with variable name
#' @param intercept used to indicate whether intercept should be considered, 0 means no intercept and 1 means with intercept
#' \itemize{
#' \item 0 means to exclude the intercept in the model
#' \item 1 means to include the intercept in the model
#' }
#' @param priorstate the column index for priorstate
#' @param currentstate the column index for currentstate
#' @param R the rank
#' @param Gamma.start the initial value for G
#' @param eps the tolerance for convergence; the default is 10^-5
#'
#' @return a list of outputs:
#' \itemize{
#' \item Alpha: the final A matrix
#' \item Gamma: the final G matrix
#' \item Beta: the coefficient matrix for variables involved in reduced rank
#' \item Dcoe: the coefficient matrix for the fixed variables
#' \item Dsderr: the standard error matrix for the fixed variables
#' \item Dpval: the p-value matrix for the fixed variables
#' \item coemat: the overall coefficient matrix
#' \item niter: the iteration number to get converged
#' \item df: the degrees of freedom
#' \item loglik: the final loglikelihood
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

rrmultinom=function(data,rrfactor,ffactor=NULL,intercept=0, priorstate,currentstate, R, Gamma.start, eps = 1e-5){
  #the input data have to start from priorstate
  #rrfactor is a vector of index for predictors, it could be column index or column name
  if(missing(rrfactor)) stop("please indicate the factors apply to reduced rank")

  data = as.data.frame(data)
  if(is.numeric(priorstate)){
    pri =  data[,priorstate]
  }else{
    pri = data[,match(priorstate,colnames(data))]
  }

  if(is.numeric(currentstate)){
    curr = data[,currentstate]
  }else{
    curr = data[,match(currentstate,colnames(data))]
  }

  if(is.numeric(rrfactor)){
    pred = as.matrix(data[,rrfactor])
  } else {
    pred = as.matrix(data[,match(rrfactor,colnames(data))])
  }

  if(is.null(ffactor)){
    fpred=NULL
  }else{
    if(is.numeric(ffactor)){
      fpred= as.matrix(data[,ffactor])
    } else {
      fpred= as.matrix(data[,match(ffactor,colnames(data))])
    }
  }
  n=nrow(data)
  p=ncol(pred)
  if(intercept==1){fpred=cbind(rep(1,n),fpred)}
  q=ncol(fpred)
  if(is.null(q)) q=0
  #y=apply(matrix(curr,n,1), 1, expand)
  #zyo=cbind(t(y),pred,fpred,pri)

  pnum=length(unique(pri))         #get the number of prior states
  ptrans=rep(0,pnum)               #store the transition number for each state
  obstrans=rep(0,n)                #store the transitions for each observation
  for(pi in 1:pnum){
    trannumber=length(unique(curr[pri==pi]))
    ptrans[pi]=trannumber
    obstrans[pri==pi]=trannumber
    curr[pri==pi]=curr[pri==pi]-min(curr[pri==pi])+1        #modified curr
  }

  updata=list()
  updata$pri=pri;updata$curr=curr;updata$pred=pred
  updata$fpred=fpred;updata$obstrans=obstrans     #pri,curr,pred,fpred,obstrans

  zyo=cbind(pri,curr,pred,fpred,obstrans)

  K=sum(ptrans)-pnum
  if(R > min(p,K)) stop("rank should not be larger than the minimal dimension of coefficient matrix")

  library(nnet)
  B0=matrix(0,p,K)
  Bcol=c(0,cumsum(ptrans-1))                               #index for B to update
  for (i in 1:pnum){
    resp0=as.factor(curr[pri==i])
    pred0=as.matrix(pred[pri==i,,drop=FALSE])
    rlev=min(curr[pri==i])
    resp01=relevel(resp0,ref=rlev)
    data0=as.data.frame(cbind(resp01,pred0))
    fit=multinom(resp01 ~ 0+pred0, data = data0)
    B0[,(Bcol[i]+1):(Bcol[i+1])] = t(summary(fit)$coefficients)
  }

  svd.B=svd(B0,nv=K)
  Gamma.start = (diag(sqrt(svd.B$d)) %*% t(svd.B$v)[1:p,])[1:R,]
  Gamma.start = matrix(Gamma.start,R,K)
  Gamma.iter = Gamma.start

  A0 = (svd.B$u %*% diag(sqrt(svd.B$d)))[,1:R]
  if (R>1){
    norm.A0 = apply(A0^2,2, function(x) sqrt(sum(x)))
  } else {norm.A0 = sqrt(sum(A0^2))}
  norm.A0.mat = matrix(norm.A0,p,R,byrow=TRUE)
  iniA = A0/norm.A0.mat

  Diter=matrix(0,q,K)

  iter <- 0
  prev.loglik <- 0
  Delta = 10
  #while(abs(Delta) > eps & iter <=200) {
  while(abs(Delta) > eps & iter <=100) {
    iter <- iter + 1
    Anew = Aupdate(Dfix=Diter,Gamma=Gamma.iter, Adata=updata, R=R,ptrans=ptrans,eps=eps)#iniA=iniA,pri=pri,curr=curr,pred=pred, fpred=fpred,
    A.iter = Anew$NewA
    loglikeA=Anew$loglikeA
    svdA=svd(A.iter%*%Gamma.iter,nv=K)
    Ai <- (svdA$u %*% diag(sqrt(svdA$d)))[,1:R]
    if (R>1){
      norm.Ai <- apply(Ai^2,2, function(x) sqrt(sum(x)))
    } else {norm.Ai <- sqrt(sum(Ai^2))}
    norm.Ai.mat <- matrix(norm.Ai,p,R,byrow=TRUE)
    A.iter <- Ai/norm.Ai.mat

    Gnew=Gupdate(A=A.iter,Gdata=updata,ptrans=ptrans)#pri=pri,curr=curr,pred=pred,fpred=fpred
    if(q==0){
      Diter=matrix(0,q,K)
    }else{
      Diter=Gnew$NewG[1:q,,drop=FALSE]
    }

    Gamma.iter = Gnew$NewG[(q+1):(q+R),,drop=FALSE]
    loglikeK=Gnew$loglikeK

    B.iter=A.iter%*%Gamma.iter
    svd <- svd(B.iter)
    Gamma.final <- (diag(sqrt(svd$d)) %*% t(svd$v))[1:R,,drop=FALSE]
    Gamma.iter <- matrix(Gamma.final,R,K)

    coemat=rbind(B.iter,Diter)

    loglikeBi=apply(zyo, 1, function(x) derivativeB(B=coemat,ptrans=ptrans,x)) ##pri,curr,pred,fpred,obstrans
    loglikeB=sum(loglikeBi)

    Delta=prev.loglik-loglikeB
    prev.loglik=loglikeB

    print(iter)
    print(prev.loglik)
  }

  if(q==0){
    Dcoe=NULL
    Dsderr=NULL
    zwpart=NULL
    pwpart=NULL
  }else{
    Dcoe=Diter
    Dsderr=Gnew$sderr[1:q,]
    zwpart=Diter/Dsderr
    pwpart = (1 - pnorm(abs(zwpart), 0, 1)) * 2
  }

  Alpha = A.iter
  Gamma.final = Gamma.iter
  B=B.iter
  svd.B <- svd(B)
  Gamma.final <- (diag(sqrt(svd.B$d)) %*% t(svd.B$v))[1:R,]
  Gamma.final <- matrix(Gamma.final,R,K)
  Alpha <- (svd.B$u %*% diag(sqrt(svd.B$d)))[,1:R]
  if (R>1){
    norm.Alpha <- apply(Alpha^2,2, function(x) sqrt(sum(x)))
  } else {norm.Alpha <- sqrt(sum(Alpha^2))}
  norm.Alpha.mat <- matrix(norm.Alpha,p,R,byrow=TRUE)
  Alpha <- Alpha/norm.Alpha.mat
  Gamma.final <- Gamma.final * matrix(norm.Alpha,R,K)

  coematfin=rbind(Dcoe,B)

  return(list(Alpha = Alpha, Gamma = Gamma.final, Beta = B, Dcoe=Dcoe, Dsderr=Dsderr, Dpval=pwpart,coemat=coematfin, niter = iter, df = R*(p+K-R)+q*K, loglik = prev.loglik))
}

