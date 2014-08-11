#' Summary of the fitted steps from chpt command
#' 
#' Gives the summary of the fitted steps based on crieria or model number
#' @import
#' @param object Object of class chpt returned by the chpt() command
#' @param type Either "bic" (default), "aic" or and integer between 1 and number models fit.
#' @return coefficients  - A matrix with Step size estimate, standard error, t-statistics and p-value
#' @examples
#' \dontrun{
#' # See example in chpt()
#' }
#' @seealso plot.chpt, get.location
summary.chpt<-function(object,type=c("bic","aic")){
  # type can be "bic", "aic" (in which case the corresponding model is chosen) or
  # integer between 1 and number of models specifying a particular model
  # object is of class chpt  -- output from chpt() command
  rlist<-object$results
  if(length(type)>1){
    type<-type[1]
  }
  if (type %in% c("bic","aic")){
    if(type=="bic"){
      type<-order(sapply(rlist,function(z){z$bic}))[1]
    }else{
      if(type=="aic"){
        type<-order(sapply(rlist,function(z){z$aic}))[1]
      }
    }
  }else{
    if(type<1 | type>length(rlist)){
      stop("type has to be one of : aic, bic or and integer between 1 and ",length(rlist))
    }
    
  }
  model<-rlist[[type]]
  betas<-model$betas
  ses<-sqrt(diag(model$betas.vcov))
  tstat<-model$tstat
  p.value<-model$pval
  coefs<-cbind(betas,ses,tstat,p.value)
  colnames(coefs)<-c("Step Size","Stdandard Error","t-statistic","P-Value")
  rownames(coefs)<-c("Base",paste("Step - ",1:(length(betas)-1),sep=""))
  coefficients<-round(coefs,3)
  return(coefficients)
}

#' Plot the fitted step
#' 
#' Plot the fitted step 
#' @import
#' @param object Object of class chpt returned by the chpt() command
#' @param type Either "bic" (default), "aic" or and integer between 1 and number models fit.
#' @param pdfname Name of pdf file with .pdf extension, Default is NULL i.e. no pdf is created
#' @return Plots the steps from the specified model
#' @examples
#' \dontrun{
#' # See example in chpt()
#' }
#' @seealso summary.chpt, get.location
plot.chpt<-function(object,type=c("bic","aic"),pdfname=NULL){
  rlist<-object$results
  wlist<-object$chpt0
  y<-wlist$y
  times<-wlist$times
  if(length(type)>1){
    type<-type[1]
  }
  if(type=="bic"){
    type<-order(sapply(rlist,function(z){z$bic}))[1]
  }else{
    if(type=="aic"){
      type<-order(sapply(rlist,function(z){z$aic}))[1]
    }
  }
  model<-rlist[[type]]
  fit<-model$fit
  locs<-which(diff(fit)!=0)+1
  fit<-cumsum(model$betas)
  ylims<-range(c(y,fit))
  stepfit<-stepfun(times[locs],fit)
  if(!is.null(pdfname)){
    pdf(file=pdfname,width=14,height=7)
  }
  xlims<-range(times)
  plot(times,y,type="l",col="lightgray",ylim=ylims,xlim=xlims,xlab="Time",ylab="Base Pair",main="Unwinding Steps")
  plot(stepfit,add=TRUE,do.points=FALSE,col=2,ylim=ylims,xlim=xlims,lwd=1.5,xlab="Time",ylab="Base Pair",main="Unwinding Steps")
  grid()
  if(!is.null(pdfname)){
    dev.off()
  }
}

#' Get model characteristics
#' 
#' Get model characteristics
#' @import
#' @param object Object of class chpt returned by the chpt() command
#' @param type Either "bic" (default), "aic" or and integer between 1 and number models fit.
#' @return Returns a list with the following elements
#' \itemize{
#' \item betas - Estimate of the step size
#' \item betas.vcov - Estimate of the variance covariance matrix of the step size estimates
#' \item tstat - The t-statistics for the step estimates
#' \item fit - the fitted  RNA trace
#' \item res - the residuals
#' \item sigma2 - Estimated variance of the underlying noise
#' \item pval - the P-value for the steps
#' \item aic - The AIC statistic for the fit
#' \item bic - The BIC statistic for the fit
#' \item AR-order - Order of the auto-regression model fitted to the noise
#' \item AR-coefs - The coefficients from the autoregression model fitted to the noise
#' }
#' @examples
#' \dontrun{
#' # See example in chpt()
#' }
#' @seealso plot.chpt, summary.chpt, get.location
get.model<-function(object,type=c("bic","aic")){
  rlist<-object$results
  if(length(type)>1){
    type<-type[1]
  }
  if (type %in% c("bic","aic")){
    if(type=="bic"){
      type<-order(sapply(rlist,function(z){z$bic}))[1]
    }else{
      if(type=="aic"){
        type<-order(sapply(rlist,function(z){z$aic}))[1]
      }
    }
  }else{
    if(type<1 | type>length(rlist)){
      stop("type has to be one of : aic, bic or and integer between 1 and ",length(rlist))
    }
    
  }
  model<-rlist[[type]]
  return(model)
}

#' Get step locations 
#' 
#' Get step locations 
#' @import
#' @param object Object of class chpt returned by the chpt() command
#' @param type Either "bic" (default), "aic" or and integer between 1 and number of models fit.
#' @return Returns the step locations of the selected model
#' @examples
#' \dontrun{
#' # See example in chpt()
#' }
get.location<-function(object,type=c("bic","aic")){
  rlist<-object$results
  wlist<-object$chpt0
  y<-wlist$y
  times<-wlist$times
  if(length(type)>1){
    type<-type[1]
  }
  if (type %in% c("bic","aic")){
    if(type=="bic"){
      type<-order(sapply(rlist,function(z){z$bic}))[1]
    }else{
      if(type=="aic"){
        type<-order(sapply(rlist,function(z){z$aic}))[1]
      }
    }
  }else{
    if(type<1 | type>length(rlist)){
      stop("type has to be one of : aic, bic or and integer between 1 and ",length(rlist))
    }
    
  }
  model<-rlist[[type]]
  fit<-model$fit
  locs<-which(diff(fit)!=0)+1
  return(locs)
}

#' Simulated RNA trace
#'
#' A dataset containing simulated RNA unwinding trace
#' for 2932 time points.
#'
#' @format A data frame with 2932 rows and 3 variables
#' 
#' The variables are as follows:
#' \itemize{
#'   \item EXT Optical tweezer extension data
#'   \item Time Time points of measurement
#'   \item Force Force measured
#' }
#' @name RNA
NULL