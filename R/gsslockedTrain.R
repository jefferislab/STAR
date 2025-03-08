.mkDiscreteLockedTrain <- function(lockedTrain,
                                   bw
                                   ) {
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
  result <- data.frame(Count=Count,
                       Time=Time)
  attr(result,"nRef") <- nRef
  attr(result,"testFreq") <- testFreq
  attr(result,"bwV") <- bwV
  attr(result,"CCH") <- lTh$CCH
  return(result)
}




#'Function to Smooth a lockedTrain Object and Related Methods: The Smoothing
#'Spline Approach
#'
#'Smooths a \code{lockedTrain} object using a smoothing spline
#'(\code{\link[gss]{gssanova}} or \code{\link[gss]{gssanova0}}) with the
#'Poisson family after binning the object.
#'
#'\code{gsslockedTrain} calls internally \code{\link[gss]{gssanova}} while
#'\code{gsslockedTrain0} calls \code{\link[gss]{gssanova0}}. See the respective
#'documentations and references therein for an explanation of the differences.
#'\code{gsslockedTrain} and \code{gsslockedTrain0} essentially generate a
#'smooth version of the histogram obtained by \code{\link{hist.lockedTrain}}.
#'The Idea is to build the histogram first with a "too" small bin width before
#'fitting a regression spline to it with a Poisson distribution of the observed
#'counts.
#'
#'@aliases gsslockedTrain gsslockedTrain0 print.gsslockedTrain
#'print.gsslockedTrain0 summary.gsslockedTrain summary.gsslockedTrain0
#'plot.gsslockedTrain plot.gsslockedTrain0
#'@param lockedTrain a \code{\link{lockedTrain}} object.
#'@param bw the bin width (in s) used to generate the observations on which the
#'gss fit will be performed. See details below.
#'@param x an \code{gsslockedTrain} or a \code{gsslockedTrain0} object.
#'@param object an \code{gsslockedTrain} or a \code{gsslockedTrain0} object.
#'@param xlim a numeric (default value supplied). See \code{\link{plot}}.
#'@param ylim a numeric (default value supplied). See \code{\link{plot}}.
#'@param xlab a character (default value supplied). See \code{\link{plot}}.
#'@param ylab a character (default value supplied). See \code{\link{plot}}.
#'@param main a character (default value supplied). See \code{\link{plot}}.
#'@param lwd line width used to plot the estimated density. See
#'\code{\link{plot}}.
#'@param col color used to plot the estimated density. See \code{\link{plot}}.
#'@param \dots in \code{gsslockedTrain}, respectively \code{gsslockedTrain0},
#'the \dots{} are passed to the internally called \code{\link[gss]{gssanova}},
#'repectively \code{\link[gss]{gssanova0}}. Not used in
#'\code{print.gsslockedTrain} and \code{summary.gsslockedTrain} and their
#'counterparts for \code{gsslockedTrain0} objects. Passed to \code{\link{plot}}
#'in \code{plot.gsslockedTrain} and \code{plot.gsslockedTrain0}.
#'@return A list of class \code{gsslockedTrain}, respectively
#'\code{gsslockedTrain0}, is returned by \code{gsslockedTrain}, respectively
#'\code{gsslockedTrain0}. These lists have the following components:
#'
#'\code{print.gsslockedTrain} returns the result of \code{\link[gss]{print}}
#'applied to the \code{gssanova} object generated by \code{gsslockedTrain} and
#'stored in the the component \code{gssFit} of its argument. The same goes for
#'\code{print.gsslockedTrain0}.
#'
#'\code{summary.gsslockedTrain} returns the result of
#'\code{\link[gss]{summary.gssanova}} applied to the \code{gssanova} object
#'generated by \code{gsspsth} and stored in the component \code{gssFit} of its
#'argument. The same goes for \code{summary.gsslockedTrain0}.
#'@returnItem gssFit the \code{gss} object generated by
#'\code{\link[gss]{gssanova}} or \code{\link[gss]{gssanova0}}.
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
#'\code{\link[gss]{gssanova}}, \code{\link[gss]{gssanova0}}
#'@references Gu C. (2002) \emph{Smoothing Spline ANOVA Models}. Springer.
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
#'slt1.3 <- gsslockedTrain(lt1.3)
#'plot(slt1.3)
#'## do some check on the gss fit
#'summary(slt1.3)
#'
#'## do the same with gsslockedTrain0
#'slt1.3 <- gsslockedTrain0(lt1.3)
#'plot(slt1.3)
#'## do some check on the gss fit
#'summary(slt1.3)
#'}
gsslockedTrain <- function(lockedTrain,
                           bw=0.001,
                           ...)
### smooths a lockedTrain object
### uses the gss model with the Poisson family after
### binning the object.
{
  df <- .mkDiscreteLockedTrain(lockedTrain,bw)
  PoissonFit <- gssanova(Count ~ Time, family="poisson", data=df, ...)

  result <- list(gssFit=PoissonFit,
                 Time=df$Time,
                 nRef=attr(df,"nRef"),
                 testFreq=attr(df,"testFreq"),
                 bwV=attr(df,"bwV"),
                 CCH=attr(df,"CCH"),
                 call=match.call())
  
  class(result) <- "gsslockedTrain"
  result
}

gsslockedTrain0 <- function(lockedTrain,
                            bw=0.001,
                            ...)
### smooths a lockedTrain object
### uses the gss model with the Poisson family after
### binning the object.
{
  df <- .mkDiscreteLockedTrain(lockedTrain,bw)
  PoissonFit <- gssanova0(Count ~ Time, family="poisson", data=df, ...)

  result <- list(gssFit=PoissonFit,
                 Time=df$Time,
                 nRef=attr(df,"nRef"),
                 testFreq=attr(df,"testFreq"),
                 bwV=attr(df,"bwV"),
                 CCH=attr(df,"CCH"),
                 call=match.call())
  
  class(result) <- "gsslockedTrain0"
  result
}

gssObj.gsslockedTrain <- function(object, ...) {
  object$gssFit
}

gssObj.gsslockedTrain0 <- gssObj.gsslockedTrain

summary.gsslockedTrain <- function(object, ...) {
  summary(gssObj(object))
}

summary.gsslockedTrain0 <- summary.gsslockedTrain

print.gsslockedTrain <- function(x, ...) {
  print(gssObj(x))
}

print.gsslockedTrain0 <- print.gsslockedTrain

plot.gsslockedTrain <- function(x,
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

  gssFitP <- predict(gssObj(x),newdata=data.frame(Time=x$Time),se.fit=TRUE)
  Y <- exp(gssFitP$fit)/x$nRef/x$bwV
  Yp <- exp((gssFitP$fit+1.96*gssFitP$se.fit))/x$nRef/x$bwV
  Ym <- exp((gssFitP$fit-1.96*gssFitP$se.fit))/x$nRef/x$bwV
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

plot.gsslockedTrain0 <- plot.gsslockedTrain
