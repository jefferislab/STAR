#'Function to Smooth a lockedTrain Object and Related Methods: The Penalized
#'Regression Spline Approach
#'
#'Smooths a \code{lockedTrain} object using a \code{gam} model with the Poisson
#'family after binning the object.
#'
#'\code{gamlockedTrain} essentially generates a smooth version of the histogram
#'obtained by \code{\link{hist.lockedTrain}}. The Idea is to build the
#'histogram first with a "too" small bin width before fitting a regression
#'spline to it with a Poisson distribution of the observed counts.
#'
#'@aliases gamlockedTrain print.gamlockedTrain summary.gamlockedTrain
#'plot.gamlockedTrain
#'@param lockedTrain a \code{\link{lockedTrain}} object.
#'@param bw the bin width (in s) used to generate the observations on which the
#'gam fit will be performed. See details below.
#'@param bs the type of splines used. See \code{\link[mgcv]{s}}.
#'@param k the dimension of the basis used to represent the smooth psth. See
#'\code{\link[mgcv]{s}}.
#'@param x an \code{gamlockedTrain} object.
#'@param object an \code{gamlockedTrain} object.
#'@param xlim a numeric (default value supplied). See \code{\link{plot}}.
#'@param ylim a numeric (default value supplied). See \code{\link{plot}}.
#'@param xlab a character (default value supplied). See \code{\link{plot}}.
#'@param ylab a character (default value supplied). See \code{\link{plot}}.
#'@param main a character (default value supplied). See \code{\link{plot}}.
#'@param lwd line width used to plot the estimated density. See
#'\code{\link{plot}}.
#'@param col color used to plot the estimated density. See \code{\link{plot}}.
#'@param \dots additional arguments passed to \code{\link[mgcv]{gam}} in
#'\code{gamlockedTrain}. Not used in \code{print.gamlockedTrain} and
#'\code{summary.gamlockedTrain}. Passed to \code{\link{plot}} in
#'\code{plot.gamlockedTrain}.
#'@return A list of class \code{gamlockedTrain} is returned by
#'\code{gamlockedTrain}. This list has the following components:
#'
#'\code{print.gamlockedTrain} returns the result of
#'\code{\link[mgcv]{print.gam}} applied to the component \code{gamFit} of its
#'argument.
#'
#'\code{summary.gamlockedTrain} returns the result of
#'\code{\link[mgcv]{summary.gam}} applied to the component \code{gamFit} of its
#'argument.
#'@returnItem gamFit the \code{\link[mgcv]{gamObject}} generated.
#'@returnItem Time the vector of bin centers.
#'@returnItem nRef the number of spikes in the reference train. See
#'\code{\link{hist.lockedTrain}}.
#'@returnItem testFreq the mean frequency of the test neuron. See
#'\code{\link{hist.lockedTrain}}.
#'@returnItem bwV the vector of bin widths used.
#'@returnItem CCH a logical which is \code{TRUE} if a cross-intensity was
#'estimated and \code{FALSE} in the case of an auto-intensity.
#'@returnItem call the matched call.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{lockedTrain}}, \code{\link{plot.lockedTrain}},
#'\code{\link[mgcv]{gam}}
#'@references Wood S.N. (2006) \emph{Generalized Additive Models: An
#'Introduction with R}. Chapman and Hall/CRC Press.
#'@keywords models smooth regression
#'@examples
#'
#'\dontrun{
#'## load e070528spont data set
#'data(e070528spont)
#'## create a lockedTrain object with neuron 1 as reference
#'## and neuron 3 as test up to lags of +/- 250 ms
#'lt1.3 <- lockedTrain(e070528spont[[1]],e070528spont[[3]],laglim=c(-1,1)*0.25)
#'## look at the cross raster plot
#'lt1.3
#'## build a histogram of it using a 10 ms bin width
#'hist(lt1.3,bw=0.01)
#'## do it the smooth way
#'slt1.3 <- gamlockedTrain(lt1.3)
#'plot(slt1.3)
#'## do some check on the gam fit
#'summary(slt1.3)
#'gam.check(gamObj(slt1.3))
#'}
gamlockedTrain <- function(lockedTrain,
                         bw=0.001,
                         bs="cr",
                         k=100,
                         ...)
### smooths a lockedTrain object
### uses the gam model with the Poisson family after
### binning the object.
{
  ## check that lockedTrain is a "lockedTrain" object
  if (!inherits(lockedTrain,"lockedTrain"))
    stop("lockedTrain should be a \"lockedTrain\" object.")
  ## make a histogram without generating a plot
  lTh <- hist(lockedTrain,bw=bw,plot=FALSE)
  Time <- lTh$mids
  nRef <- lTh$nRef
  testFreq <- lTh$testFreq
  bwV <- diff(lTh$breaks)
  Count <- lTh$density*nRef*bwV
  PoissonFit <- gam(Count ~ s(Time,k=k,bs=bs),family=poisson(),...)

  result <- list(gamFit=PoissonFit,
                 Time=Time,
                 nRef=nRef,
                 testFreq=testFreq,
                 bwV=bwV,
                 CCH=lTh$CCH,
                 call=match.call())
  
  class(result) <- "gamlockedTrain"
  result
}

gamObj.gamlockedTrain <- function(object, ...) {
  object$gamFit
}

summary.gamlockedTrain <- function(object, ...) {
  summary(gamObj(object))
}

print.gamlockedTrain <- function(x, ...) {
  print(gamObj(x))
}

plot.gamlockedTrain <- function(x,
                              xlab, ylab, main, xlim, ylim,
                              col,lwd,
                              ...) {

  if (missing(xlab)) 
    xlab <- "Time (s)"
  if (missing(ylab)) 
    ylab <- "PDF of a Test Neuron Spike (1/s)"
  if (missing(main)) {
    nameList <- deparse(x$call[["lockedTrain"]])
    if (x$CCH) main <- paste(nameList, "smoothed cross-intensity")
    else main <- paste(nameList, "smoothed auto-intensity")
  }

  gamFitP <- predict(gamObj(x),type="response",se.fit=TRUE)
  Y <- gamFitP$fit/x$nRef/x$bwV
  Yp <- (gamFitP$fit+1.96*gamFitP$se.fit)/x$nRef/x$bwV
  Ym <- (gamFitP$fit-1.96*gamFitP$se.fit)/x$nRef/x$bwV
  X <- x$Time

  if (missing(xlim)) xlim <- range(X)
  if (missing(ylim)) ylim <- c(min(c(Ym,x$testFreq)),max(c(Yp,x$testFreq)))
  if (missing(col)) col <- 2
  if (missing(lwd)) lwd <- 2
  
  plot(X,Y,type="n",
       xlab=xlab,ylab=ylab,
       xlim=xlim,ylim=ylim,
       main=main,
       ...)
  abline(h=x$testFreq)
  abline(v=0)
  lines(X,Yp,lty=2)
  lines(X,Ym,lty=2)
  lines(X,Y,col=col,lwd=lwd)
  
}
