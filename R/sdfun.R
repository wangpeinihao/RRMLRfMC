#' sdfun
#'
#' This function is used get the standard error matrix from bootstrap method
#' It returns the matrices of standard error and p-value for the coefficient matrix
#'
#' @param I a U by U incidence matrix with elements; U is the number of states; I(i,j)=1 if state j can be accessed from state i in one step and 0 otherwise
#' @param z1 a n by p matrix with covariates involved in the dimension reduction(DR), n is the number of subjects, p is the number of covariates involved in DR
#' @param z2 a n by q matrix with study covariates (not in dimension reduction), q is the number of study covariates
#' @param T a M by 3 state matrix,
#' \itemize{
#' \item the first column is a subject number between 1,..,n;
#' \item the second column is time;
#' \item the third column is the state occupied by subject in column 1 at time indicated in column 2
#' }
#' @param R the rank
#' @param eps the tolerance for convergence; the default is 10^-5
#' @param B the bootstrap number
#' @param tpoint a matrix has two columns with the participants' visit information about timeline
#' @param ref a vector of reference categories
#'
#' @return a list of outputs:
#' \itemize{
#' \item coe: the coefficient matrix of the original data
#' \item sd: the standard error matrix
#' \item pvalue: the p-value matrix
#' }
#'
#'
#'
#'
#'
#'
#'

sdfun=function(I,z1=NULL,z2=NULL,T,R,eps = 1e-5,B,tpoint=NULL,ref){  #tpoint is a M-n by 2 matrix with the patient id and visit time points(length 4) as the 1st and 2nd columns

  if(is.null(z1) & is.null(z2)) stop("please input the covariates")
  if(is.null(z1)){p=0}else{n=nrow(z1); p=ncol(z1)}
  if(is.null(z2)){q=0}else{n=nrow(z2); q=ncol(z2)}

  U=nrow(I)
  tn=sum(apply(I, 1, sum)!=0)     #number of transient states
  wt=which(apply(I, 1, sum)!=0)   #index of transient states
  an=nrow(I)-tn                   #number of absorbing states
  K=sum(I)-tn                     #columns of the coefficient matrix
  bootcoe=matrix(0,B,(p+q)*K)
  Bi=1

  while(Bi<=B){
    if((is.null(tpoint)) | isTRUE(nrow(tpoint)!=(nrow(T)-n))){     #bootstrap patients

    ind=sample((1:n),n,replace = TRUE)
    znew=matrix(0,n,p+q)
    Tnew=matrix(0,0,3)
    for (j in 1:n) {
      nj=ind[j]
      znew[j,]=c(z1[nj,],z2[nj,])
      Tsub=T[which(T[,1]==nj),]
      Tsub[,1]=rep(j, nrow(Tsub))
      Tnew=rbind(Tnew,Tsub)
    }
    if(p==0){z1new=NULL}else{z1new=znew[,(1:p)]}
    if(q==0){z2new=NULL}else{z2new=znew[,((p+1):(p+q))]}
    } else {                                            #simulate MC for each patient
      tstate=matrix(0,nrow(T)-n,2)                            #a MX2 matrix with total prior and current states
      Mc=0
      for(ni in 1:n){
        subT=T[which(T[,1]==ni),,drop=FALSE]
        mc=nrow(subT)-1
        tstate[(Mc+1):(Mc+mc),1]=subT[,3][-(mc+1)]
        tstate[(Mc+1):(Mc+mc),2]=subT[,3][-1]
        Mc=Mc+mc
      }
      pstate=sort(unique(tstate[,1]))
      pl=length(pstate)
      cstate=sort(unique(tstate[,2]))
      cl=length(cstate)
      Pmatrix=matrix(0,pl,cl)
      colnames(Pmatrix) = cstate
      rownames(Pmatrix) = pstate
      propt=list()
      for(i in 1:pl){             #to get the frequency table
          subtstate=tstate[which(tstate[,1]==i),,drop=FALSE]
          curstate=sort(unique(subtstate))
        	colind=match(curstate,cstate)
        	Pmatrix[i,colind]=table(subtstate[,2])/sum(table(subtstate[,2]))
        	propt[[i]]=table(subtstate[,2])/sum(table(subtstate[,2]))
      }
      yearmax=max(tpoint[,2])                                                  #get the max timepoint
      abstate=setdiff(unique(tstate[,2]),unique(tstate[,1]))                   #get the absorbing states

       Tnew=rep(0,3)
       for(pi in 1:n){
        vsind=which(tpoint[,1]==pi)                #rows index for pi-th subject
        ivdate=tpoint[vsind,2,drop=FALSE][1,] #get the first visit year
        istate=tstate[vsind,1,drop=FALSE][1,] #get the initial state
        sig=1
        pstate=c()          #store the generated prior state
        cstate=c()          #store the generated current state
        iyear=ivdate
        pstate=istate
        ntran=0               #store the transition counts for T[,2]
        Tdata=c(pi,0,istate)      #store the generated T for ith patient
        while(sig==1){
          cstate=sample(as.numeric(names(propt[[pstate]])),1,prob=propt[[pstate]])   #generate the current state,
          iyear=iyear+1
          pstate=cstate
          ntran=ntran+1
          irow=c(pi,ntran,pstate)
          Tdata=rbind(Tdata,irow)
          if(iyear==(yearmax) | cstate %in% abstate) sig=0   #the iyear and absorbing state may need be modified for other data
        }
        Tnew=rbind(Tnew,Tdata)
       }
      Tnew=Tnew[-1,]
      z1new=z1;z2new=z2
    }

    tryCatch({
    rrmodel=rrmultinom(I,z1=z1new,z2=z2new,T=Tnew,R,eps=1e-5,ref)
    niter=rrmodel$niter
    }, error=function(e){niter<<-300})

    if(niter < 200){ # only record the results when niter < 200
      bootcoe[Bi,]=as.vector(rrmodel$coemat)
      Bi=Bi+1
    }
  }

  oricoe=rrmultinom(I,z1,z2,T,R,eps=1e-5,ref)$coemat
  cVB=apply(bootcoe,2,mean)
  cVB=t(replicate(B, cVB))
  sdVB=1/(B-1)*t(bootcoe-cVB)%*%(bootcoe-cVB)
  fsdVB=sqrt(matrix(diag(sdVB),p+q,K))

  zwpart=oricoe/fsdVB
  pwpart = (1 - stats::pnorm(abs(zwpart), 0, 1))*2
  return(list(coe=oricoe,sd=fsdVB,pvalue=pwpart))
}

