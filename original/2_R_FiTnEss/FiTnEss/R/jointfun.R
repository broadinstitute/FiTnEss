#jointfun4

#2. multiple Nta categories

nblgnfun2<-function(q,lp,sigma,nta){
  lbd1<-1/(rlnorm(100000,meanlog=lp,sdlog=sigma))
  ndist<-rnbinom(100000,nta,lbd1)
  pvec=sapply(q,function(x){sum(x>=ndist,na.rm = TRUE)/length(ndist)})
  return(pvec)
}

lownbfun3<-function(q,nta){
  pvec=pnbinom(q,nta,0.7)
  return(pvec)
}

jointfun<-function(q,lambda,lp,sigma,nta){
  pvec=lambda*nblgnfun2(q,lp,sigma,nta)+(1-lambda)*lownbfun3(q,nta)
  return(pvec)
}
