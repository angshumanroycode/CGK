cgk.preset<-function(n,d,parameters=NULL,B.parameters=1000){
  if(is.null(parameters))
    parameters=PARAMETERS(n,d,B.parameters)
  return(CGKPRESET(n,d,parameters))
}

cgk<-function(X,parameters=NULL,preset=NULL,B.parameters=1000){
  if(is.null(preset)){
    xdim=dim(X)
    n=xdim[1]
    d=xdim[2]
    preset=cgk.preset(n,d,parameters,B.parameters)
  }
  return(CGK(X,preset))
}

cgk.test.preset<-function(n,d,parameters=NULL,B.parameters=1000,
                          B=1000,preset=NULL){
  if(is.null(preset))
    preset=cgk.preset(n,d,parameters,B.parameters)
  return(CGKTESTPRESET(n,d,preset,B))
}

cgk.test<-function(X,method=c("sum","max","fdr"),alpha=0.05,
                    B.parameters=1000,B=1000,parameters=NULL,
                    preset=NULL,test.preset=NULL){
  if(is.null(test.preset)){
    xdim=dim(X)
    n=xdim[1]
    d=xdim[2]
    test.preset=cgk.test.preset(n,d,parameters,B.parameters,
                                B,preset)
  }
  tstat=CGK(X,test.preset[[1]])
  ndist=test.preset[[2]]
  result=list()
  l=1
  if("sum" %in% method){
    tstat.sum=sum(tstat)
    ndist.sum=rowSums(ndist)
    pvalue=(sum(tstat.sum<ndist.sum)+1)/(length(ndist.sum)+1)
    reject=pvalue<alpha
    subresult=list(reject,pvalue)
    names(subresult)=c("reject","p.value")
    result[[l]]=subresult
    l=l+1
  }
  if("max" %in% method){
    tstat.max=max(tstat)
    ndist.max=apply(ndist,1,max)
    pvalue=(sum(tstat.max<ndist.max)+1)/(length(ndist.max)+1)
    reject=pvalue<alpha
    subresult=list(reject,pvalue)
    names(subresult)=c("reject","p.value")
    result[[l]]=subresult
    l=l+1
  }
  if("fdr" %in% method){
    np=NCOL(ndist)
    nr=NROW(ndist)
    tstat.fdr=matrix(rep(tstat,nr),nrow=nr,byrow=T)
    pvalsort=sort(colMeans(tstat.fdr<ndist))
    Levels=alpha*(1:np)/np
    reject=sum(pvalsort<Levels)>0
    subresult=list(reject)
    names(subresult)=c("reject")
    result[[l]]=subresult
  }
  names(result)=method
  return(result)
}

