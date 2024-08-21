wt_kernel_est = function(data, dose_vals, outcome_col, dose_col, kern_wts, bw=NULL) {
  kern_info = tryCatch({
    m1 <- my_locpol(as.formula(paste0(outcome_col, " ~ ", dose_col)), 
                    data=data, weig=kern_wts, 
                    deg=1, kernel=locpol::gaussK, bw=bw)
    est = approx(x=m1$lpFit[,m1$X], y=m1$lpFit[,m1$Y], xout=dose_vals)$y
    list(est=est, bw=m1$bw)
  },
  error=function(cond) {
    list(est=rep(NA, length(dose_vals)), bw=NA)
  })
  return(kern_info)
}

my_locpol <- function(formula, data, weig=rep(1,nrow(data)),
                      bw=NULL, kernel=locpol::gaussK, deg=1, xeval=NULL,xevalLen=100)
  ##	
{
  ##  checking
  stopifnot(nrow(data)==length(weig))
  ## compute result
  res <- list()
  res$bw <- bw
  res$KName <- match.call()
  res$kernel <- kernel
  res$deg <- deg
  res$xeval <- xeval
  ## get info from formula
  res$mf <- model.frame(formula,data)
  datCla <- attr(attr(res$mf, "terms"),"dataClasses") 
  varNames <- names(datCla)[datCla=="numeric"] 
  stopifnot(length(varNames)==2) 
  res$Y <- varNames[1] 
  res$X <- varNames[2] 
  ##	sort x's
  xo <- order(res$mf[,res$X])
  res$mf <- res$mf[xo,]
  res$weig <- weig[xo]
  ## xeval
  if( is.null(xeval) )
    res$xeval <- seq(min(res$mf[,res$X]),max(res$mf[,res$X]),len=xevalLen)
  else 
    res$xeval <- sort(xeval)
  if( is.null(res$bw) ) 
    res$bw <- locpol::regCVBwSelC(data[,res$X],data[,res$Y],res$deg,
                                  res$kernel,res$weig)
  ## regression estimation
  res$lpFit <- locpol::locPolSmootherC(res$mf[,res$X], res$mf[,res$Y], res$xeval, 
                                       res$bw, res$deg, res$kernel, DET = TRUE, res$weig )
  names(res$lpFit)[] <- c(res$X,res$Y,paste(res$Y,1:deg,sep=""),"xDen")
  res$lpFit$xDen <- res$lpFit$xDen^(1/(deg+1))/(nrow(data)*res$bw)
  return(res)
}