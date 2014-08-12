#' Get quantiles in the left and right window
#' 
#' Internal function. Not to be called by the user
get.eta<-function(x,win){
  # Getting  eta for a single window
  qleft<-runquantile(x,win,c(0.25,0.5,0.75),align="right")    # quantiles of left window of size win
  qright<-runquantile(x,win,c(0.25,0.5,0.75),align="left")    # quantiles of right window of size win
  qrange<-apply(runquantile(x,2*win+1,c(0,1)),1,diff)         # range of the window of size 2*win+1 centered at the point
  nr<-apply((qleft-qright)^2,1,mean)                          # mean squared difference of the quantiles
  eta<-nr/(qrange^2)                                          # the statistic eta
  return(eta)
}

#' Internal function to get statistics for various windows
#' 
#' Internal function. Not to be called by the user
get.zstat<-function(y,times,w){
  # Getting eta for multiple windows.
  # The function returns a matrix with each column corresponding to a window
  # get eta values for all w values
  mzstat<-sapply(w,function(k){get.eta(y,k)})
  colnames(mzstat)<-w
  return(list(y=y,w=w,times=times,mzstat=mzstat))
}

#' Internal function to get split trace by change point
#' 
#' Internal function. Not to be called by the user
split.vector<-function(y,locs){
  # splitting a trace into segments at the step locations
  locs<-ceiling(locs)
  n<-length(y)
  if(any(locs==1)){
    locs<-locs[locs!=1]
  }
  if(any(locs==n)){
    locs<-locs[locs!=n]
  }
  nn<-diff(c(1,locs,n+1))
  id<-rep(1:(length(locs)+1),nn)
  id<-id[1:length(y)]
  yy<-split(y,id)
  return(yy)
}

#' Internal function to get design matrix to be used in the regression
#' 
#' Internal function. Not to be called by the user
get.xmat<-function(n,locs){
  # n - number of data points
  # locs - location of steps
  locs<-ceiling(locs)
  if(any(locs==1)){
    locs<-locs[locs!=1]
  }
  if(any(locs==n)){
    locs<-locs[locs!=n]
  }
  nn<-diff(c(1,locs,n+1))
  id<-rep(1:(length(locs)+1),nn)
  id<-id[1:n]
  id2<-1:n
  id2<-split(id2,id)
  x<-matrix(0,nrow=n,ncol=length(id2))
  id3<-sapply(id2,max)
  for(i in 1:length(id2)){
    if(i==1){
      x[,i]<-1
    }else{
      x[-c(1:id3[i-1]),i]<-1
    }
  }
  return(x)
}


#' Internal function to get perform efficient linear regression
#' 
#' Internal function. Not to be called by the user
sp.lm<-function(x,y,fx,fy){
  xtx<-crossprod(fx)                          # matrix operation x'x
  xty<-crossprod(fx,as.numeric(fy))           # matrix operation x'y
  betas<-solve(xtx,xty)                       # matrix operation (x'x)^(-1)x'y
  fit<-as.numeric(fx%*%betas)              
  res<-fy-fit
  dof<-length(fy)-dim(betas)[1]
  sigma2<-var(res)
  betas.vcov<-sigma2*solve(xtx)
  tstat<-abs(betas/sqrt(diag(betas.vcov)))
  pval<-2*pt(as.numeric(tstat),df=dof,lower.tail=FALSE)
  fit1<-as.numeric(x%*%betas)
  res1<-y-fit1
  res<-res+mean(res1)
  llhood<-sum(log(dnorm(res,sd=sqrt(sigma2))))
  aic<-(-2*llhood)+(2*(ncol(fx)))
  bic<-(-2*llhood)+(log(nrow(fx))*ncol(fx))
  return(list(betas=betas,betas.vcov=betas.vcov,tstat=tstat,fit=fit1,res=res,sigma2=sigma2,pval=pval,aic=aic,bic=bic))
}

#' Internal function to get perform robust ar estimation
#' 
#' Internal function. Not to be called by the user
robust.ar<-function(x,p){
  # fit robust ar to the residuals
  lma<-rlm(x~1)
  w<-lma$w
  mu<-lma$coef
  x<-x-mu
  z<-embed(x,p+1)
  z<-z[,-1]
  x<-x[-c(1:p)]
  w<-w[-c(1:p)]
  dat1<-cbind(x,z)
  colnames(dat1)<-c("x",paste("x",1:ncol(z),sep=""))
  data1<-as.data.frame(dat1)
  lm1<-lm(x~.-1,data=data1,weights=w)
  ar1<-lm1$coef
  x.mean<-mu
  return(list(ar=ar1,x.mean=x.mean,w=w))
}
