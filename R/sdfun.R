#' sdfun
#'
#' This function is used get the standard error matrix from bootstrap method
#' It returns the matrices of standard error and p-value for the coefficient matrix
#'
#' @param dat the dataset
#' @param M the bootstrap number
#' @param rank the rank of the reduced rank model
#' @param rrfactor the index for variables in reduced rank part, it can be a numeric vector or a vector with the variable name
#' @param ffactor the index for variables in fixed part, it can be a numeric vector or a vector with the variable name
#' @param priorstate the column index number for the priorstate, it is a numeric number
#' @param currentstate the column index number for the priorstate, it is a numeric number
#' @param intercept this parameter indicates whether the intercept is included in the model, it has options 0 or 1 with default as 0
#' \itemize{
#' \item 0 means to exclude the intercept in the model
#' \item 1 means to include the intercept in the model
#' }
#'
#' @return a list contains resultand matrices with \strong{sd} as the standard error matrix and \strong{pvalue} as the p-value matrix
#'
#' @export
#'
#' @references
#'
#' @examples
#'
#'
#'
#'
#'

sdfun=function(dat,M,rank,rrfactor=c(3,4,6:11),ffactor=c(12,5),priorstate=1,currentstate=2,intercept=0){
  np=length(unique(dat$priorstate)) #priorstate number
  nt=0  #transition number
  for(ip in 1:np){
    subdat=dat[dat$priorstate==ip,]
    nt=nt+length(unique(subdat$currentstate))
  }
  nt=nt-np
  if (intercept==0){
    colinB=(length(rrfactor)+length(ffactor))*nt
  } else {colinB=(length(rrfactor)+length(ffactor)+1)*nt}
  VB=matrix(NA,M,colinB)
  time=rep(0,M)
  Mi=0
  i=1
  while (Mi < M) {
    print(Mi)
    Bdata=bootdata(dat)
    sdata=Bdata[,-(1:2)]
    R=rank
    eps=1e-5
    cons=rep(1,nrow(sdata))
    sdata=cbind(sdata,cons)
    start_time <- Sys.time()

    tryCatch({
      rrmodel=rrmultinom(sdata,rrfactor=rrfactor,ffactor=ffactor,intercept=intercept,priorstate=priorstate,currentstate=currentstate, R=R, eps = 1e-5)
      niter=rrmodel$niter
    }, error=function(e){niter<<-300})  #if error occurs, let niter=300

    end_time <- Sys.time()
    if(niter < 200){ # only record the results when niter < 200
      VB[i,]=as.vector(rrmodel$coemat)
      itime=end_time - start_time
      time[i]=itime
      i=i+1
      Mi=Mi+1
    }
  }
  sddim1=nrow(rrmodel$coemat)
  sddim2=ncol(rrmodel$coemat)
  cVB=apply(VB,2,mean)
  cVB=t(replicate(M, cVB))
  sdVB=1/(M-1)*t(VB-cVB)%*%(VB-cVB)
  fsdVB=matrix(diag(sdVB),sddim1,sddim2)

  rrmodelori=rrmultinom(dat[-(1:2)],rrfactor=rrfactor,ffactor=ffactor,intercept=intercept,priorstate=priorstate,currentstate=currentstate, R=R, eps = 1e-5)
  coef=rrmodelori$coemat
  zwpart=coef/fsdVB
  pwpart = (1 - pnorm(abs(zwpart), 0, 1)) * 2
  return(list(sd=fsdVB,pvalue=pwpart))
}
