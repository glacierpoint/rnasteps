#' Fits multiple step models for a given trace
#' 
#' Fits multiple step models for a given trace
#' @import
#' @param y RNA unwinding trace
#' @param times Time points at which RNA trace was recorded
#' @param w Sequence of window sizes
#' @param cutoff Upper percentile for choosing plausible step locations. Default is 0.9 i.e. only top 10\% highest values are choosen by default
#' @param cor Default is TRUE indicating that the noise is correlated. If FALSE the noise is assumed uncorrelated
#' @return Returns a list with the following elements
#' \itemize{
#' \item results - is a list of length equal to number of models fitted. Each element in this list comprises of elements returned in get.model
#' \item chpt0 - A list with y, times, w and zstat - the statistics based on which the model is fitted
#' }
#' @examples
#' \dontrun{
#' 
#' w<-c(seq(10,90,by=10),seq(100,1000,by=25))           
#' y<-RNA[,1]
#' times<-RNA[,2]   
#' chpt1<-chpt(y,times,w)
#' plot.step(chpt1)
#' # to save as a pdffile
#' # plot.step(chpt1,pdfname="RNAFIG1.pdf")
#' summary(chpt1) # to get bic fit
#' summary(chpt1,type="aic") # to get aic fit
#' summary(chpt1,type=5) # to get fit summary of model number 5
#' plot(chpt1)
#' mymodel<-get.model(chpt1)
#' names(mymodel)
#' get.location(chpt1)
#' }
#' @seealso get.model, plot.chpt, summary.chpt, get.location
#' @references
#' Arunajadai SG, Cheng W (2013) Step Detection in Single-Molecule Real Time Trajectories Embedded in Correlated Noise. PLoS ONE 8(3): e59279. 
#' 
#' Cheng W, Arunajadai SG, Moffitt JR, Tinoco I Jr, Bustamante C. Single-base pair unwinding and asynchronous RNA release by the hepatitis C virus NS3 helicase. Science. 2011 Sep 23;333(6050):1746-9. doi: 10.1126/science.1206023.
chpt<-function(y,times,w) UseMethod("chpt")

ChPt<-function(...){
  print("The use of ChPt is deprecated. Please use the function chpt")
}

chpt.default<-function(y,times,w,cutoff=0.9,cor=TRUE){
  # cutoff is the upper percentile for choosing plausible step locations i.e. only top 10% highest values are choosen by default
  #-----------------------------------------------------------------------------------------------------------------------------
  # chpt0 computes the zstat that is uses to find the plausible set of change points - stored in mzstat
  chpt0<-get.zstat(y,times,w)
  mzstat<-chpt0$mzstat
  #-----------------------------------------------------------------------------------------------------------------------------
  # Getting the initial set of step locations
  # peaks1 are the initial set of step locations
  MU<-apply(mzstat,1,mean)
  snr<-quantile(MU,cutoff)
  if(length(y)>=1000){
    mpeak1<-msPeak(msSet(mzstat),FUN="search",span.supsmu="cv",snr=snr)
    mpeak2<-msAlign(mpeak1,snr=snr)
    peaks1<-mpeak2$peak.class[,1]
  }else{
    peaks1<-which(MU>snr)
  }
  # peaks1 contains the inital list of plausible step locations
  #-----------------------------------------------------------------------------------------------------------------------------
  split.y<-split.vector(y,peaks1)     # split y at the initial step locations
  xmat<-get.xmat(length(y),peaks1)    # design matrix corresponding to the initial step locations
  ind.mat<-xmat                       # ind.mat is matrix that speeds up matrix assignment
  ar.fit<-ar(diff(y),method="mle")    # Getting the order of the AR noise 
  arp<-ar.fit$order
  len<-cumsum(sapply(split.y,length))+1
  len<-c(1,len[-ncol(xmat)])
  len<-len+arp
  index<-1:length(y)
  for(i in 1:length(len)){
    ind.mat[index>len[i],i]<-2
  }
  fx<-xmat                            # fx mat will be the matrix corrsponding to the Cochrane-Orcutt filtered matrix
  # Set Storage
  RESULTS<-list()
  count<-1
  ind<-0
  if(cor==TRUE){
    while(ind==0){
      lm2<-sp.lm(xmat,y,xmat,y)                         # simple linear regression
      lm2[[length(lm2)+1]]<-arp
      names(lm2)[length(lm2)]<-"AR-order"
      if(arp!=0){
        #ar.fit<-ar(lm2$res,order=arp,aic=FALSE,method="mle")     # Fit AR model to noise for pre-determined order above
        ar.fit<-ar(lm2$res,aic=FALSE,method="mle")     # Fit AR model to noise for pre-determined order above
        phi<-ar.fit$ar              
        mu<-ar.fit$x.mean
        # Begin Filtering: Transform xmat to fx and y to fy for regression acoounting correlation
        y2<-y-mu
        n<-length(y2)
        phi<-c(1,-phi)
        n2<-length(phi)-1
        y2<-c(rep(0,n2),y2)
        fy<-filter(y2,phi,method="convolution",sides=1)
        fy<-fy[-c(1:n2)]
        nc<-length(phi)
        sumcoefs<-cumsum(phi)
        sumcoefs<-sumcoefs-(sumcoefs*mu)
        fx[ind.mat==2]<-sumcoefs[nc]
        fx[ind.mat==1]<-rep(sumcoefs,ncol(fx))
        # End Filtering:  at this stage the fx and fy contain the filtered matrix and vector
        lm2<-sp.lm(xmat,y,fx,fy)                      # regression after taking into account the correlation
        lm2[[length(lm2)+1]]<-arp
        names(lm2)[length(lm2)]<-"AR-order"
        lm2[[length(lm2)+1]]<-phi
        names(lm2)[length(lm2)]<-"AR-Coefs"
      }
      RESULTS[[count]]<-lm2
      pval<-lm2$pval[-1]                              # Always ignore the intercept 
      id<-rev(order(pval))[1]                         # is the column to be removed
      id<-id+1
      if(ncol(fx)!=1){
        # here we remove the corresponding columns in all the relevenat matrices
        xmat<-xmat[,-id,drop=FALSE]
        ind.mat<-ind.mat[,-id,drop=FALSE]
        fx<-fx[,-id,drop=FALSE]
        count<-count+1
      }else{
        ind<-1
      }
      gc()            # this is a command to free up available memory
    }
  }else{
    while(ind==0){
      lm2<-sp.lm(xmat,y,xmat,y)                         # simple linear regression
      RESULTS[[count]]<-lm2
      pval<-lm2$pval[-1]                              # Always ignore the intercept 
      id<-rev(order(pval))[1]                         # is the column to be removed
      id<-id+1
      if(ncol(fx)!=1){
        # here we remove the corresponding columns in all the relevenat matrices
        xmat<-xmat[,-id,drop=FALSE]
        ind.mat<-ind.mat[,-id,drop=FALSE]
        fx<-fx[,-id,drop=FALSE]
        count<-count+1
      }else{
        ind<-1
      }
      gc()            # this is a command to free up available memory
    }
  }
  RESULTS<-list(results=RESULTS,chpt0=chpt0)
  class(RESULTS)<-"chpt"
  return(RESULTS)
}
