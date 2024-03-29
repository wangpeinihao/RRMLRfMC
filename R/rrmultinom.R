#' rrmultinom
#'
#' This function is used to fit the reduced rank multinomial logistic regression
#'     for markov chain
#'
#' @param I a U by U incidence matrix with elements; U is number of states; I(i,j)=1
#'     if state j can be accessed from state i in one step and 0 otherwise
#' @param z1 a n by p matrix with covariates involved in the dimension reduction(DR),
#'     n is the number of subjects, p is the number of covariates involved in DR
#' @param z2 a n by q matrix with study covariates (not in dimension reduction),
#'     q is the number of study covariates
#' @param T a M by 3 state matrix,
#' \itemize{
#' \item the first column is a subject number between 1,..,n;
#' \item the second column is time;
#' \item the third column is the state occupied by subject in column 1 at time indicated
#'     in column 2
#' }
#' @param R the rank
#' @param eps the tolerance for convergence; the default is 10^-5
#' @param ref a vector of reference categories; the default is NULL and if NULL
#'     is used, the function will use the first category as the reference category
#'     for each row
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
#' \item converge: three possible values with 0 means fail to converge, 1 means
#'     converges, and 2 means the maximum iteration is achieved
#' }
#'
#'@export
#'
#'@examples
#'# generate the Markov chain
#'U=7
#'I1=I2=I3=rep(1,7)
#'I4=c(0,0,0,1,1,1,1)
#'I5=I6=I7=rep(0,7)
#'I=rbind(I1,I2,I3,I4,I5,I6,I7)
#'# prepare the data
#'data=cogdat
#'n=length(unique(data[,1]))
#'M=nrow(data)+n
#'Mc=0
#'z=matrix(0,n,9)
#'colnames(z)=colnames(data)[5:13]
#'T=matrix(0,M,3)
#'for(i in 1:n){
#'  subdat=data[which(data[,1]==i),,drop=FALSE]
#'  z[i,]=subdat[1,5:13]
#'  mc=nrow(subdat)
#'  T[(Mc+1):(Mc+mc+1),1]=i
#'  T[(Mc+1):(Mc+mc+1),2]=0:mc
#'  T[(Mc+1):(Mc+mc+1),3]=c(subdat[1,3],subdat[,4])
#'  Mc=Mc+mc+1
#'}
#'#z1=z[,c(1:3),drop=FALSE]
#'z2=z[,4,drop=FALSE]
#'# fit the model with rank 1
#'rrmultinom(I,z1=NULL,z2,T,1,eps=9,ref=c(1,1,1,4))
#'

rrmultinom=function(I,z1=NULL,z2=NULL,T,R,eps = 1e-5,ref=NULL){
  #U=number of states; I=incidence matrix;z1=DR variables;
  #z2=study variables;T=state matrix;

  I=as.matrix(I)
  U=nrow(I)
  tn=sum(apply(I, 1, sum)!=0)     #number of transient states
  #if(is.null(ref)) ref=na.omit(apply(I,1,function(x) which(x!=0)[1]))
  if(is.null(ref)) {
    rind = which(apply(I, 1, sum)!=0)
    ref = c(0, length(rind))
    for(ri in 1:length(rind)){
      (rr = rind[ri])
      (ref[ri] = which(I[rr,]!=0)[1])
    }
  }
  if(!is.null(ref)){ if(length(ref)!=tn) stop("please have a correct length for reference vector")}
  wt=which(apply(I, 1, sum)!=0)   #index of transient states
  an=nrow(I)-tn                   #number of absorbing states
  K=sum(I)-tn                     #columns of the coefficient matrix

  if(is.null(z1) &is.null(z2)) stop("please input the covariates")

  if(is.null(z1)){
    z2=as.matrix(z2)
    n=nrow(z2)  #number of subjects
    q=ncol(z2)  #number of study covariates
    M=nrow(T)
    TN=M-n      #number of total transitions
    mdata=matrix(0,TN,4+q)  #a matrix:transition, observation,prior,current,study covariates
    mp=0        #pointer for create mdata
    for(i in 1:n){
      subdat=T[T[,1]==i,]
      mc=nrow(subdat)-1
      mdata[(mp+1):(mp+mc),1]=(mp+1):(mp+mc)
      mdata[(mp+1):(mp+mc),2]=i
      mdata[(mp+1):(mp+mc),3]=subdat[,3][1:mc]
      mdata[(mp+1):(mp+mc),4]=subdat[,3][2:(mc+1)]
      mdata[(mp+1):(mp+mc),(5:(4+q))]=t(replicate(mc,z2[i,]))
      mp=mp+mc
    }

    B0=matrix(0,q,K)
    SE=matrix(0,q,K)
    bp=0                           #pointer to the coefficient matrix column
    fvalue=0                       #final converged value

    for (i in 1:tn){
      ti=wt[i]                     #transient states index
      resp0=as.factor(mdata[mdata[,3]==ti,4])
      bcl=sum(I[ti,])-1            #number of coefficient columns for ith MLR (with )
      pred0=as.matrix(mdata[mdata[,3]==ti,(5:(4+q)),drop=FALSE])
      rlev=ref[i]
      resp01=stats::relevel(resp0,ref=as.character(rlev))
      data0=as.data.frame(cbind(resp01,pred0))
      fit=nnet::multinom(resp01 ~ 0+pred0, data = data0)
      B0[,(bp+1):(bp+bcl)] = t(summary(fit)$coefficients)
      SE[,(bp+1):(bp+bcl)] = t(summary(fit)$standard.errors)
      bp=bp+bcl
      fvalue=fvalue+summary(fit)$value
    }
    
    zwpart=B0/SE
    pwpart = (1 - stats::pnorm(abs(zwpart), 0, 1)) * 2
    return(list(Alpha=NULL,Gamma=NULL,Beta=B0,Dcoe=B0,Dsderr=SE,Dpval=pwpart,coemat=B0,niter=1,df=NULL,loglik=fvalue))
  } else {

    z1=as.matrix(z1)
    if(is.null(z2)){
      q=0
      }else{
        z2=as.matrix(z2)
        q=ncol(z2)    #number of study covariates
      }
    n=nrow(z1)    #number of subjects
    p=ncol(z1)    #number of covariate in DR
    M=nrow(T)
    TN=M-n        #number of total transitions

    mdata=matrix(0,TN,4+q+p)  #create a matrix with columns: transition, observation,prior,current,study covariates
    mp=0        #pointer for create mdata
    if(q==0){
      for(i in 1:n){
        subdat=T[T[,1]==i,]
        mc=nrow(subdat)-1
        mdata[(mp+1):(mp+mc),1]=(mp+1):(mp+mc)
        mdata[(mp+1):(mp+mc),2]=i
        mdata[(mp+1):(mp+mc),3]=subdat[,3][1:mc]
        mdata[(mp+1):(mp+mc),4]=subdat[,3][2:(mc+1)]
        mdata[(mp+1):(mp+mc),(5:(4+p))]=t(replicate(mc,z1[i,]))
        mp=mp+mc
      }
    }else{
      for(i in 1:n){
        subdat=T[T[,1]==i,]
        mc=nrow(subdat)-1
        mdata[(mp+1):(mp+mc),1]=(mp+1):(mp+mc)
        mdata[(mp+1):(mp+mc),2]=i
        mdata[(mp+1):(mp+mc),3]=subdat[,3][1:mc]
        mdata[(mp+1):(mp+mc),4]=subdat[,3][2:(mc+1)]
        mdata[(mp+1):(mp+mc),(5:(4+p))]=t(replicate(mc,z1[i,]))
        mdata[(mp+1):(mp+mc),((4+p+1):(4+p+q))]=t(replicate(mc,z2[i,]))
        mp=mp+mc
      }
    }
    

    if(R > min(p,K)) stop("rank should not be larger than the minimal dimension of coefficient matrix")

    #get the initial value for iteration
    B0=matrix(0,p,K)
    bp=0
    for (i in 1:tn){
      ti=wt[i]                     #transient states index
      resp0=as.factor(mdata[mdata[,3]==ti,4])
      bcl=sum(I[ti,])-1            #number of coefficient columns for ith MLR (with )
      pred0=as.matrix(mdata[mdata[,3]==ti,(5:(4+p)),drop=FALSE])
      rlev=ref[i]
      resp01=stats::relevel(resp0,ref=as.character(rlev))
      data0=as.data.frame(cbind(resp01,pred0))
      fit=nnet::multinom(resp01 ~ 0+pred0, data = data0)
      B0[,(bp+1):(bp+bcl)] = t(summary(fit)$coefficients)
      bp=bp+bcl
    }
    svd.B=svd(B0,nv=K)
    Gamma.start = (diag(sqrt(svd.B$d))%*%t(svd.B$v)[1:p,])[1:R,]
    Gamma.start = matrix(Gamma.start,R,K)
    Gamma.iter = Gamma.start
    A0 = (svd.B$u %*% diag(sqrt(svd.B$d)))[,1:R]
    if (R>1){
      norm.A0 = apply(A0^2,2, function(x) sqrt(sum(x)))
    } else {norm.A0 = sqrt(sum(A0^2))}
    norm.A0.mat = matrix(norm.A0,p,R,byrow=TRUE)
    iniA = A0/norm.A0.mat
    Diter=matrix(0,q,K)

    #newton raphson algorithm
    iter=0;prev.loglik=0;Delta=100;miter=200
    while(abs(Delta) > eps & iter <=miter) {
      tryCatch({
        Anew = Aupdate(Dfix=Diter,Gamma=Gamma.iter, Adata=mdata[,-c(1,2)], R=R,p=p,q=q,I=I, eps=eps,refA=ref)#iniA=iniA,pri=pri,curr=curr,pred=pred, fpred=fpred,
        A.iter = Anew$NewA
        loglikeA=Anew$loglikeA
        # svdA=svd(A.iter%*%Gamma.iter,nv=K)
        # Ai <- (svdA$u %*% diag(sqrt(svdA$d)))[,1:R]
        # if (R>1){
        #   norm.Ai <- apply(Ai^2,2, function(x) sqrt(sum(x)))
        # } else {norm.Ai <- sqrt(sum(Ai^2))}
        # norm.Ai.mat <- matrix(norm.Ai,p,R,byrow=TRUE)
        # A.iter <- Ai/norm.Ai.mat
      }, error=function(e){return(list(convergent=0, iter=1000))})

      Gnew=Gupdate(A=A.iter,Gdata=mdata[,-c(1,2)],p,q,I,refG=ref)#pri=pri,curr=curr,pred=pred,fpred=fpred
      if(q==0){
        Diter=matrix(0,q,K)
      }else{
        Diter=Gnew$NewG[1:q,,drop=FALSE]
      }

      Gamma.iter = Gnew$NewG[(q+1):(q+R),,drop=FALSE]
      loglikeK=Gnew$loglikeK

      B.iter=A.iter%*%Gamma.iter
      # svd <- svd(B.iter)
      # Gamma.final <- (diag(sqrt(svd$d)) %*% t(svd$v))[1:R,,drop=FALSE]
      # Gamma.iter <- matrix(Gamma.final,R,K)
      coemat=rbind(B.iter,Diter)

      loglikeBi=apply(mdata[,-c(1,2)], 1, function(x) derivativeB(B=coemat,I,x,refd=ref)) ##pri,curr,pred,fpred,obstrans
      loglikeB=sum(loglikeBi)

      Delta=prev.loglik-loglikeB
      prev.loglik=loglikeB
      iter=iter+1

      print(iter)
      print(prev.loglik)
    }
    if(iter > miter){convergent=2}else{convergent=1} #indicate the status of convergent

    if(q==0){
      Dcoe=NULL
      Dsderr=NULL
      zwpart=NULL
      pwpart=NULL
    }else{
      Dcoe=Diter
      Dsderr=Gnew$sderr[1:q,]
      zwpart=Diter/Dsderr
      pwpart = (1 - stats::pnorm(abs(zwpart), 0, 1)) * 2
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

    return(list(Alpha = Alpha, Gamma = Gamma.final, Beta = B, Dcoe=Dcoe, Dsderr=Dsderr, Dpval=pwpart,coemat=coematfin, niter = iter, df = R*(p+K-R)+q*K, loglik = prev.loglik,converge=convergent))
  }
}

