"%frt%" <- function(refTrain,testTrain) frt(refTrain,testTrain)



#'Computes Forward Recurrence Times from Two transformedTrain Objects
#'
#'Computes the (transformed) time differences between spikes of a
#'\code{refTrain} and the (next) ones of a \code{testTrain}. Both
#'\code{refTrain} and \code{testTrain} should be \code{transformedTrain}
#'objects.
#'
#'When two spike trains have been time transformed using \emph{the same}
#'procedure, which does make one of the trains (the \code{testTrain}) the
#'realization a homogeneous Poisson process with rate 1, the elapsed time
#'between the spikes of the other train (\code{refTrain}) and the ones of
#'\code{testTrain} should be exponentially distributed with rate 1. These
#'elapsed times are returned by \code{frt}.
#'
#'@aliases frt %frt%
#'@param refTrain a \code{transformedTrain} object.
#'@param testTrain a \code{transformedTrain} object.
#'@return An object of class \code{frt} containing the elapsed times.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{transformedTrain}}, \code{\link{plot.frt}},
#'\code{\link{summary.frt}}, \code{\link{mkGLMdf}}
#'@keywords models
#'@examples
#'
#'\dontrun{
#'## Let us consider neuron 1 of the CAL2S data set
#'data(CAL2S)
#'CAL2S <- lapply(CAL2S,as.spikeTrain)
#'CAL2S[["neuron 1"]]
#'renewalTestPlot(CAL2S[["neuron 1"]])
#'summary(CAL2S[["neuron 1"]])
#'## Make a data frame with a 4 ms time resolution
#'cal2Sdf <- mkGLMdf(CAL2S,0.004,0,60)
#'## keep the part relative to neuron 1, 2 and 3 separately
#'n1.cal2sDF <- cal2Sdf[cal2Sdf$neuron=="1",]
#'n2.cal2sDF <- cal2Sdf[cal2Sdf$neuron=="2",]
#'n3.cal2sDF <- cal2Sdf[cal2Sdf$neuron=="3",]
#'## remove unnecessary data
#'rm(cal2Sdf)
#'## Extract the elapsed time since the second to last and
#'## third to last for neuron 1. Normalise the result. 
#'n1.cal2sDF[c("rlN.1","rsN.1","rtN.1")] <- brt4df(n1.cal2sDF,"lN.1",2,c("rlN.1","rsN.1","rtN.1"))
#'## load mgcv library
#'library(mgcv)
#'## fit a model with a tensorial product involving the last
#'## three spikes and using a cubic spline basis for the last two
#'## To gain time use a fixed df regression spline
#'n1S.fitA <- gam(event ~ te(rlN.1,rsN.1,bs="cr",fx=TRUE) + rtN.1,data=n1.cal2sDF,family=binomial(link="logit"))
#'## transform time
#'N1.Lambda <- transformedTrain(n1S.fitA)
#'## check out the resulting spike train using the fact
#'## that transformedTrain objects inherit from spikeTrain
#'## objects
#'N1.Lambda
#'## Use more formal checks
#'summary(N1.Lambda)
#'plot(N1.Lambda,which=c(1,2,4,5),ask=FALSE)
#'## Transform spike trains of neuron 2 and 3
#'N2.Lambda <- transformedTrain(n1S.fitA,n2.cal2sDF$event)
#'N3.Lambda <- transformedTrain(n1S.fitA,n3.cal2sDF$event)
#'## Check interactions
#'summary(N2.Lambda %frt% N1.Lambda)
#'summary(N3.Lambda %frt% N1.Lambda)
#'plot(N2.Lambda %frt% N1.Lambda,ask=FALSE)
#'plot(N3.Lambda %frt% N1.Lambda,ask=FALSE)
#'}
#'
frt <- function(refTrain,
                testTrain
                ) {

  class(refTrain) <- NULL
  class(testTrain) <- NULL
  refTrain <- refTrain[!is.na(refTrain)]
  testTrain <- testTrain[!is.na(testTrain)]
  result <- sapply(refTrain,
                   function(t) {
                     goodOnes <- testTrain >= t
                     if (sum(goodOnes) == 0) return(NA)
                     min(testTrain[goodOnes])-t
                   }
                   )
  result <- result[!is.na(result)]
  class(result) <- "frt"
  result
}



#'Plots and Summarizes frt Objects.
#'
#'\code{plot.frt} generates interactively (by default) 2 plots, the survivor
#'function with confidence intervals and the Berman's test with confidence
#'bands. \code{summary.frt} generates a concise summary of \code{frt} objects.
#'It is mostly intended for use in batch processing situations where a decision
#'to stop with the current model or go on with a more complicated one must be
#'made automatically.
#'
#'If the reference and test (transformed) spike trains used in the
#'\code{\link{frt}} call which generated \code{x} (or \code{object}) are not
#'correlated (and if the transformed test train is indeed homogeneous Poisson
#'with rate 1), the elements of \code{x} (or \code{object}) should be iid
#'realizations of an exponential with rate 1. Two test plots are generated by
#'\code{plot.frt} in the same way as the corresponding ones (testing the same
#'thing) of \code{\link{plot.transformedTrain}}.
#'
#'The same correspondence holds between \code{summary.frt} and
#'\code{\link{summary.transformedTrain}}.
#'
#'@aliases plot.frt summary.frt
#'@param x a \code{transformedTrain} object.
#'@param object a \code{transformedTrain} object.
#'@param which if a subset of the plots is required, specify a subset of the
#'numbers \code{1:2}.
#'@param main title to appear above the plots, if missing the corresponding
#'element of \code{caption} will be used.
#'@param caption Default caption to appear above the plots or, if \code{main}
#'is given, bellow it
#'@param ask logical; if \code{TRUE}, the user is \emph{ask}ed before each
#'plot, see \code{\link{par}(ask=.)}.
#'@param \dots additional arguments passed to \code{\link{plot}}.
#'@return \code{summary.frt} returns a vector with named elements stating if
#'the Berman's test is passed with a 95\% and a 99\% confidence.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{transformedTrain}}, \code{\link{frt}},
#'\code{\link{mkGLMdf}}
#'@keywords models
#'@examples
#'
#'\dontrun{
#'## Let us consider neuron 1 of the CAL2S data set
#'data(CAL2S)
#'CAL2S <- lapply(CAL2S,as.spikeTrain)
#'CAL2S[["neuron 1"]]
#'renewalTestPlot(CAL2S[["neuron 1"]])
#'summary(CAL2S[["neuron 1"]])
#'## Make a data frame with a 4 ms time resolution
#'cal2Sdf <- mkGLMdf(CAL2S,0.004,0,60)
#'## keep the part relative to neuron 1, 2 and 3 separately
#'n1.cal2sDF <- cal2Sdf[cal2Sdf$neuron=="1",]
#'n2.cal2sDF <- cal2Sdf[cal2Sdf$neuron=="2",]
#'n3.cal2sDF <- cal2Sdf[cal2Sdf$neuron=="3",]
#'## remove unnecessary data
#'rm(cal2Sdf)
#'## Extract the elapsed time since the second to last and
#'## third to last for neuron 1. Normalise the result. 
#'n1.cal2sDF[c("rlN.1","rsN.1","rtN.1")] <- brt4df(n1.cal2sDF,"lN.1",2,c("rlN.1","rsN.1","rtN.1"))
#'## load mgcv library
#'library(mgcv)
#'## fit a model with a tensorial product involving the last
#'## three spikes and using a cubic spline basis for the last two
#'## To gain time use a fixed df regression spline
#'n1S.fitA <- gam(event ~ te(rlN.1,rsN.1,bs="cr",fx=TRUE) + rtN.1,data=n1.cal2sDF,family=binomial(link="logit"))
#'## transform time
#'N1.Lambda <- transformedTrain(n1S.fitA)
#'## check out the resulting spike train using the fact
#'## that transformedTrain objects inherit from spikeTrain
#'## objects
#'N1.Lambda
#'## Use more formal checks
#'summary(N1.Lambda)
#'plot(N1.Lambda,which=c(1,2,4,5),ask=FALSE)
#'## Transform spike trains of neuron 2 and 3
#'N2.Lambda <- transformedTrain(n1S.fitA,n2.cal2sDF$event)
#'N3.Lambda <- transformedTrain(n1S.fitA,n3.cal2sDF$event)
#'## Check interactions
#'summary(N2.Lambda %frt% N1.Lambda)
#'summary(N3.Lambda %frt% N1.Lambda)
#'plot(N2.Lambda %frt% N1.Lambda,ask=FALSE)
#'plot(N3.Lambda %frt% N1.Lambda,ask=FALSE)
#'}
#'
plot.frt <- function(x,
                     which=1:2,
                     main,
                     caption=c("Log Survivor Function",
                       "Berman's Test"),
                     ask=TRUE,
                     ...) {

  show <- logical(2)
  show[which] <- TRUE

  if (ask) {
    op <- par(ask = TRUE)
    on.exit(par(op))
  } else {
    if (sum(show)==2) layout(matrix(1:2,nrow=2))
  }

  X <- sort(x)
  nI <- length(X)
  mainGiven <- !missing(main)
  
  if (show[1]) {
    Y <- (nI:1)/nI
    YTh <- exp(-X)
    p95 <- qbinom(0.975,nI,YTh)/nI
    m95 <- qbinom(0.025,nI,YTh)/nI
    p99 <- qbinom(0.995,nI,YTh)/nI
    m99 <- qbinom(0.005,nI,YTh)/nI
    plot(X,Y,
         type="n",log="y",
         xlab="Rescaled Forward Recurrence Time",
         ylab="Survivor Fct",
         main=ifelse(mainGiven,main,caption[1]),
         sub=ifelse(mainGiven,caption[1],""),
         ...)
    lines(X,YTh)
    lines(X,p95,lty=2)
    lines(X,m95,lty=2)
    lines(X,p99,lty=2)
    lines(X,m99,lty=2)
    lines(X,Y,col=2,lwd=2)
  }

  lambda <- 1-exp(-x)
  
  if (show[2]) {
    plot(c(0,1),c(0,1),type="n",
         xlab=expression(U[(k)]),
         ylab="Cumulative Distribution",
         main=ifelse(mainGiven,main,caption[2]),
         sub=ifelse(mainGiven,caption[2],""),
         ...
         )
    abline(a=0,b=1)
    abline(a=1.36/sqrt(nI),1,lty=2)
    abline(a=-1.36/sqrt(nI),1,lty=2)
    abline(a=1.63/sqrt(nI),1,lty=2)
    abline(a=-1.63/sqrt(nI),1,lty=2)
    lines(sort(lambda),(1:(nI))/nI,col=2,lwd=2)
  }

}

summary.frt <- function(object,...) {

  lambda <- sort(1-exp(-object))
  nI <- length(lambda)
  M <- (1:(nI))/nI - lambda
  in95 <- all(-1.36/sqrt(nI) <= M & M <= 1.36/sqrt(nI))
  in99 <- all(-1.63/sqrt(nI) <= M & M <= 1.63/sqrt(nI))
  c(in95=in95,in99=in99)
}
