#' bootdata
#'
#' This function is used to generate the bootstrap data based on the first status and the empirical probability
#' It returns a bootstrap data
#'
#' @param data the original dataset, it has to be a data frame with newid, vdata in format of 1990/9/14, priorstate,currentstate, and covariates
#'
#' @return bootstrap data
#'
#' @export
#'
#' @examples
#'
#'
#'
#'
#'

bootdata=function(data){ #the data have to be newid, vdata in format(1990/9/14), priorstate,currenstate and covariates
  person=length(unique(data$newid))
  bdata=rep(0,ncol(data))
  npri=length(unique(data$priorstate))
  propt=list()  #get the probability for sampling
  for(j in 1:npri){
    pdata=data[data$priorstate==j,]
    pj=table(pdata$currentstate)/sum(table(pdata$currentstate))
    propt[[j]]=pj
  }

  yearmax=max(apply(matrix(data$vdate,ncol = 1),1,function(x)as.numeric(substring(x,1,4))))  #get the maximum year it can reach
  abstate=setdiff(unique(data$currentstate),unique(data$priorstate))

  for(i in 1:person){
    idata=data[data$newid==i,]
    ivdate=as.numeric(substring(idata[1,2],1,4)) #get the first visit year
    istate=idata[1,3] #get the initial state
    icov=unlist(idata[1,-(1:4)])  #get the covariates for ith patient
    rowdata=rep(0,ncol(data))   #store the data for ith patient
    sig=1
    pstate=c()          #store the generated prior state
    cstate=c()          #store the generated current state
    iyear=ivdate
    pstate=istate
    while(sig==1){
      cstate=sample(as.numeric(names(propt[[pstate]])),1,prob=propt[[pstate]])   #generate the current state,
      irow=c(i,iyear,pstate,cstate,icov)
      rowdata=rbind(rowdata,irow)
      iyear=iyear+1
      pstate=cstate
      if(iyear==(yearmax+1) | cstate %in% abstate) sig=0   #the iyear and absorbing state may need be modified for other data
    }
    bdata=rbind(bdata,rowdata[-1,])
  }
  bdata=bdata[-1,]
  return(bdata)
}
