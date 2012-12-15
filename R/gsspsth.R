.mkDiscreteTrain <- function(repeatedTrain,
                             binSize
                             ) {
  
  if (!is.repeatedTrain(repeatedTrain)) 
    repeatedTrain <- as.repeatedTrain(repeatedTrain)
  bigTrain <- sort(unlist(repeatedTrain))
  minSpikeTime <- floor(min(bigTrain))
  maxSpikeTime <- ceiling(max(bigTrain))
  Time <- seq(minSpikeTime, maxSpikeTime, by = binSize)
  if (max(bigTrain) > max(Time)) 
    Time <- c(Time, max(Time) + binSize)
  H <- hist(bigTrain, breaks = Time, plot = FALSE)
  result <- data.frame(Count=H$counts,
                       Time=H$mids)
  attr(result,"minSpikeTime") <- minSpikeTime
  attr(result,"maxSpikeTime") <- maxSpikeTime
  return(result)
}

.postFitFormat <- function(gfit,
                           repeatedTrain,
                           binSize,
                           minSpikeTime,
                           maxSpikeTime
                           ) {

  nbTrials <- length(repeatedTrain)
  Time <- gfit$mf$Time
  Count <- gfit$mf$Count
  Y <- predict(gfit,
               newdata=data.frame(Time=Time),
               se.fit = TRUE)
  lambdaFct <- function(t) {
    exp(as.numeric(predict(gfit, newdata = data.frame(Time = t))))/nbTrials/binSize
  }
  LambdaFct <- function(t) {
    t <- sort(t)
    t <- t[minSpikeTime <= t & t <= maxSpikeTime]
    tL <- length(t)
    result <- numeric(tL)
    result[1] <- integrate(lambdaFct, minSpikeTime, t[1])$value
    if (tL > 1) 
      for (i in 2:tL) result[i] <- integrate(lambdaFct, 
                                             t[i - 1], t[i])$value
    cumsum(result)
  }
  result <- list(freq = exp(Y$fit)/nbTrials/binSize,
                 ciUp = exp(Y$fit + 1.96 * Y$se.fit)/nbTrials/binSize, 
                 ciLow = exp(Y$fit - 1.96 * Y$se.fit)/nbTrials/binSize,
                 breaks = c(minSpikeTime, 
                   maxSpikeTime),
                 mids = Time,
                 counts = Count,
                 nbTrials = nbTrials, 
                 lambdaFct = lambdaFct,
                 LambdaFct = LambdaFct,
                 call = match.call())
  return(result)
  
}

gsspsth0 <- function(repeatedTrain,
                     binSize = 0.025, 
                     plot = FALSE,
                     ...) {
  
  df <- .mkDiscreteTrain(repeatedTrain,
                         binSize)
  PoissonF <- gssanova0(Count ~ Time, family = "poisson", data=df, ...)
  result <- .postFitFormat(PoissonF,repeatedTrain,binSize,
                           attr(df,"minSpikeTime"),
                           attr(df,"maxSpikeTime"))
  class(result) <- "gsspsth0"
  if (plot) {
    plot(result, colCI=2)
  }
  else {
    return(result)
  }
}




#'Smooth Peri Stimulus Time Histogram Related Functions and Methods: The
#'Smoothing Spline Approach
#'
#'Function \code{gsspsth} and \code{gsspsth0} compute a smooth psth, while
#'method \code{print.gsspsth} and \code{print.gsspsth0} print and
#'\code{summary.gsspsth} or \code{summary.gsspsth0} summarize the
#'\code{gssanova} / \code{gssanova0} objects contained in the returned
#'\code{gsspsth} or \code{gsspsth0} objects, \code{plot.gsspsth} or
#'\code{plot.gsspsth0} plot them and \code{simulate.gsspsth} or
#'\code{simulate.gsspsth0} simulate data from fitted objects.
#'
#'\code{gsspsth} calls internally \code{\link[gss]{gssanova}} while
#'\code{gsspsth0} calls \code{\link[gss]{gssanova0}}. See the respective
#'documentations and references therein for an explanation of the differences.
#'For both \code{gsspsth} and \code{gsspsth0}, the raw data contained in
#'\code{repeatedTrain} are pre-processed with \code{\link{hist}} using a bin
#'size given by argument \code{binSize}. This \code{binSize} should be small
#'"enough". That is, the rate of the aggregated train created by collapsing the
#'spike times of the different trials onto a single "pseudo" spike train,
#'should not change too much on the scale of \code{binSize} (see Ventura et al
#'(2002) Sec. 4.2 p8 for more details). Argument \code{nbasis} of
#'\code{\link[gss]{gssanova}} called internally by \code{gsspsth} is set to the
#'number of bins of the histogram resulting from the preprocessing stage.
#'
#'\code{simulate.gsspsth} and \code{simulate.gsspsth0} perform exact simuations
#'of inhomogenous Poisson processes whose (time dependent) rate/intensity
#'function is accessible through the componenent \code{lambdaFct} of objects of
#'class \code{gsspsth} and \code{gsspsth0}. The inhomogenous Poisson processes
#'are simulated with the thinning method (Devroye, 1986, pp 253-256).
#'
#'@aliases gsspsth gsspsth0 print.gsspsth print.gsspsth0 summary.gsspsth
#'summary.gsspsth0 plot.gsspsth plot.gsspsth0 simulate.gsspsth0
#'simulate.gsspsth
#'@param repeatedTrain a \code{repeatedTrain} object or a list which can be
#'coerced to such an object.
#'@param binSize the bin size (in s) used to generate the observations on which
#'the gss fit will be performed. See details below.
#'@param plot corresponding argument of \code{\link{hist}}. Should a plot be
#'generated or not?
#'@param object a \code{gsspsth} or a \code{gsspsth0} object.
#'@param x a \code{gsspsth} or a \code{gsspsth0} object.
#'@param stimTimeCourse \code{NULL} (default) or a two elements vector
#'specifying the time boundaries (in s) of a stimulus presentation.
#'@param colStim the background color used for the stimulus.
#'@param colCI if not \code{NULL} (default) a confidence band is plotted with
#'the specified color; two dashed lines are plotted otherwise.
#'@param xlim a numeric (default value supplied). See \code{\link{plot}}.
#'@param ylim a numeric (default value supplied). See \code{\link{plot}}.
#'@param xlab a character (default value supplied). See \code{\link{plot}}.
#'@param ylab a character (default value supplied). See \code{\link{plot}}.
#'@param main a character (default value supplied). See \code{\link{plot}}.
#'@param lwd line width used to plot the estimated density. See
#'\code{\link{plot}}.
#'@param col color used to plot the estimated density. See \code{\link{plot}}.
#'@param nsim number of \code{repeatedTrain} objects to simulate. Defaults to
#'1.
#'@param seed see \code{\link{simulate}}.
#'@param \dots in \code{gsspsth}, respectively \code{gsspsth0}, the \dots{} are
#'passed to the internally called \code{\link[gss]{gssanova}}, repectively
#'\code{\link[gss]{gssanova0}}. In \code{plot.gsspsth} and \code{plot.gsspsth0}
#'they are passed to \code{\link{plot}} which is called internally. They are
#'not used otherwise.
#'@return When \code{plot} is set to \code{FALSE} in \code{gsspsth},
#'repectively \code{gsspsth0}, a list of class \code{gsspsth}, respectively
#'\code{gsspsth0}, is returned and no plot is generated. These list have the
#'following components:
#'
#'When \code{plot} is set to \code{TRUE} nothing is returned and a plot is
#'generated as a side effect. Of course the same occurs upon calling
#'\code{plot.gsspsth} with a \code{gsspsth} object argument or
#'\code{plot.gsspsth0} with a \code{gsspsth0}.
#'
#'\code{print.gsspsth} returns the result of \code{\link[gss]{print}} applied
#'to the \code{gssanova} object generated by \code{gsspsth} and stored in the
#'\code{\link{environment}} of both \code{lambdaFct} and \code{LambdaFct}. The
#'same goes for \code{print.gsspsth0}.
#'
#'\code{summary.gsspsth} returns the result of
#'\code{\link[gss]{summary.gssanova}} applied to the \code{gssanova} object
#'generated by \code{gsspsth} and stored in the \code{\link{environment}} of
#'both \code{lambdaFct} and \code{LambdaFct}. The same goes for
#'\code{summary.gsspsth0}.
#'
#'\code{simulate.gsspsth} and \code{simulate.gsspsth0} return a
#'\code{repeatedTrain} object if argument \code{nsim} is set to one and a list
#'of such objects if it is greater than one.
#'@returnItem freq a vector containing the instantaneous firing rate in the
#'middle of the "thin" bins used for preprocessing.
#'@returnItem ciUp a vector with the upper limit of a pointwise 95\% confidence
#'interval.  Check \code{\link[gss]{predict.ssanova}} for details.
#'@returnItem ciLow a vector with the lower limit of a pointwise 95\%
#'confidence interval.
#'@returnItem breaks a vector with 2 elements the ealiest and the latest spike
#'in \code{repeatedTrain}.
#'@returnItem mids a numeric vector with the mid points of the bins.
#'@returnItem counts a vector with the actual number of spikes in each bin.
#'@returnItem nbTrials the number of trials in \code{repeatedTrain}.
#'@returnItem lambdaFct a function of a single time argument returning the
#'estimated intensity (or instantaneous rate) at its argument.
#'@returnItem LambdaFct a function of a single time argument returning the
#'integrale of estimated intensity (or instantaneous rate) at its argument.
#'That is, the integrated intensity. \code{\link{integrate}} is used by this
#'function.
#'@returnItem call the matched call.
#'@note Most of the components of the list returned by \code{gsspsth} and
#'\code{gsspsth0} are not of direct interest for the user but they are used by,
#'for instance, \code{\link{reportHTML.repeatedTrain}}.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{psth}}, \code{\link{plot.psth}},
#'\code{\link[gss]{gssanova}}, \code{\link[gss]{gssanova0}},
#'\code{\link[gss]{summary.gssanova}}, \code{\link[gss]{summary.gssanova0}},
#'\code{\link{reportHTML.repeatedTrain}}, \code{\link{simulate}}
#'@references Gu C. (2002) \emph{Smoothing Spline ANOVA Models}. Springer.
#'
#'Ventura, V., Carta, R., Kass, R. E., Gettner, S. N. and Olson, C. R. (2002)
#'Statistical analysis of temporal evolution in single-neuron firing rates.
#'\emph{Biostatistics} \bold{3}: 1--20.
#'
#'Kass, R. E., Ventura, V. and Cai, C. (2003) Statistical smoothing of neuronal
#'data. \emph{Network: Computation in Neural Systems} \bold{14}: 5--15.
#'
#'Devroye Luc (1986) \emph{Non-Uniform Random Variate Generation}. Springer.
#'Book available in pdf format at:
#'\url{http://cg.scs.carleton.ca/~luc/rnbookindex.html}.
#'@keywords models smooth regression
#'@examples
#'
#'\dontrun{
#'## Get the e070528citronellal data set into workspace
#'data(e070528citronellal)
#'## Compute gsspsth without a plot for neuron 1
#'## using a smmothing spline with gssanova0, and default bin size of 25 ms.
#'n1CitrGSSPSTH0 <- gsspsth0(e070528citronellal[[1]])
#'## plot the result
#'plot(n1CitrGSSPSTH0,stim=c(6.14,6.64),colCI=2)
#'## get a summary of the gss fit
#'summary(n1CitrGSSPSTH0)
#'## Now take a look at the observation on which the gss
#'## was actually performed
#'plot(n1CitrGSSPSTH0$mids,n1CitrGSSPSTH0$counts,type="l")
#'## Add the estimated smooth psth after proper scaling
#'theBS <- diff(n1CitrGSSPSTH0[["mids"]])[1]
#'Y <- n1CitrGSSPSTH0$lambdaFct(n1CitrGSSPSTH0$mids)*theBS*n1CitrGSSPSTH0$nbTrials
#'lines(n1CitrGSSPSTH0$mids,Y,col=4,lwd=2)
#'
#'## check the (absence of) effect of the pre-binning by using a smaller
#'## and larger one, say 10 and 75 ms
#'n1CitrGSSPSTH0_10 <- gsspsth0(e070528citronellal[[1]],binSize=0.01)
#'n1CitrGSSPSTH0_75 <- gsspsth0(e070528citronellal[[1]],binSize=0.075)
#'## plot the "high resolution" smoothed-psth
#'plot(n1CitrGSSPSTH0_10,colCI="grey50")
#'## add to it the estimate obtained with the "low resolution" one
#'Y_75 <- n1CitrGSSPSTH0_75$lambdaFct(n1CitrGSSPSTH0_10$mids)
#'lines(n1CitrGSSPSTH0_10$mids,Y_75,col=2,lwd=2)
#'
#'## simulate data from the first fitted model
#'s1 <- simulate(n1CitrGSSPSTH0)
#'## look at the result
#'s1
#'
#'## Do the same thing with gsspsth
#'n1CitrGSSPSTH <- gsspsth(e070528citronellal[[1]])
#'plot(n1CitrGSSPSTH,stim=c(6.14,6.64),colCI=2)
#'summary(n1CitrGSSPSTH)
#'plot(n1CitrGSSPSTH$mids,n1CitrGSSPSTH$counts,type="l")
#'theBS <- diff(n1CitrGSSPSTH[["mids"]])[1]
#'Y <- n1CitrGSSPSTH$lambdaFct(n1CitrGSSPSTH$mids)*theBS*n1CitrGSSPSTH$nbTrials
#'lines(n1CitrGSSPSTH$mids,Y,col=4,lwd=2)
#'## check the (absence of) effect of the pre-binning by using a smaller
#'## and larger one, say 10 and 75 ms
#'n1CitrGSSPSTH_10 <- gsspsth(e070528citronellal[[1]],binSize=0.01)
#'n1CitrGSSPSTH_75 <- gsspsth(e070528citronellal[[1]],binSize=0.075)
#'## plot the "high resolution" smoothed-psth
#'plot(n1CitrGSSPSTH_10,colCI="grey50")
#'## add to it the estimate obtained with the "low resolution" one
#'Y_75 <- n1CitrGSSPSTH_75$lambdaFct(n1CitrGSSPSTH_10$mids)
#'lines(n1CitrGSSPSTH_10$mids,Y_75,col=2,lwd=2)
#'## simulate data from the first fitted model
#'s1 <- simulate(n1CitrGSSPSTH)
#'## look at the result
#'s1
#'}
#'
gsspsth <- function(repeatedTrain,
                    binSize = 0.025, 
                    plot = FALSE,
                    ...) {
  
  df <- .mkDiscreteTrain(repeatedTrain,
                         binSize)
  PoissonF <- gssanova(Count ~ Time, family = "poisson", data = df,
                       nbasis=dim(df)[1],...)
  result <- .postFitFormat(PoissonF,repeatedTrain,binSize,
                           attr(df,"minSpikeTime"),
                           attr(df,"maxSpikeTime"))
  class(result) <- "gsspsth"
  if (plot) {
    plot(result, colCI=2)
  }
  else {
    return(result)
  }
}

.plot_spsth <- function(x,
                        stimTimeCourse,
                        colStim,
                        colCI, 
                        xlab, ylab, main, xlim, ylim,
                        lwd, col,
                        ...) {
  
  plot(xlim, ylim, type = "n",
       xlab = xlab, ylab = ylab, main = main, 
       xlim = xlim, ylim = ylim,
       ...)
  
  if (!is.null(stimTimeCourse)) {
    rect(stimTimeCourse[1], 0, stimTimeCourse[2], ylim[2], 
         col = colStim, lty = 0)
  }
  if (is.null(colCI)) {
    lines(x$mids, x$ciLow, lty = 2)
    lines(x$mids, x$ciUp, lty = 2)
  } else {
    polygon(c(x$mids, rev(x$mids)), c(x$ciLow, rev(x$ciUp)), 
            col = colCI, border = NA)
  }
  lines(x$mids, x$freq, col = col, lwd = lwd)
  abline(h = 0)
  
}

plot.gsspsth <- function(x,
                         stimTimeCourse = NULL,
                         colStim = "grey80",
                         colCI = NULL, 
                         xlab, ylab, main, xlim, ylim,
                         lwd = 2, col = 1,
                         ...) {
  
  if (!is.null(stimTimeCourse)) {
    if (length(stimTimeCourse) != 2) 
      stop(paste(deparse(substitute(stimTimeCourse)), "should be a vector with 2 elements."))
  }
  if (missing(xlab)) xlab <- "Time (s)"
  if (missing(ylab)) ylab <- "Freq (Hz)"
  if (missing(main)) {
    nameList <- deparse(x$call[["repeatedTrain"]])
    main <- paste(nameList, "PSTH")
  }
  if (missing(xlim)) 
    xlim <- x$breaks
  if (missing(ylim)) 
    ylim <- c(0, max(x$ciUp[!is.na(x$ciUp)]))

  .plot_spsth(x,stimTimeCourse,colStim,colCI, 
              xlab, ylab, main, xlim, ylim, lwd, col,
              ...)
}

plot.gsspsth0 <- plot.gsspsth




#'Generic Function and Methods for Extracting a gss object
#'
#'Some functions of \code{STAR}, like \code{gsspsth}, \code{gsspsth0} and
#'\code{gsslockedTrain}, \code{gsslockedTrain0} perform gss fits internally and
#'keep as a list component or within the environment of a returned function the
#'result of this fit. Method \code{gssObj} extracts this \code{gss} object.
#'
#'
#'@aliases gssObj gssObj.gsspsth gssObj.gsspsth0 gssObj.gsslockedTrain
#'gssObj.gsslockedTrain0
#'@param object an object containing a \code{gssanova} or a \code{gssanova0}
#'object. Currently the result of a call to \code{\link{gsspsth}},
#'\code{\link{gsspsth0}} or to \code{\link{gsslockedTrain}},
#'\code{\link{gsslockedTrain0}}.
#'@param \dots not used for now
#'@return A \code{\link[gss]{gssanova}} or a \code{\link[gss]{gssanova0}}
#'object.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link[gss]{gssanova}}, \code{\link[gss]{gssanova0}},
#'\code{\link{gsspsth}}, \code{\link{gsspsth0}}, \code{\link{gsslockedTrain}},
#'\code{\link{gsslockedTrain0}}
#'@keywords models smooth regression
#'@examples
#'
#'##
#'
gssObj <- function(object,
                   ...) {

  UseMethod("gssObj")
}

gssObj.gsspsth <- function(object,...) {

  evalq(gfit, envir=environment(object$lambdaFct))

}


summary.gsspsth <- function(object, ...) {
  summary(gssObj(object))
}

print.gsspsth <- function(x, ...) {
  print(gssObj(x))
}

gssObj.gsspsth0 <- gssObj.gsspsth
summary.gsspsth0 <- summary.gsspsth
print.gsspsth0 <- print.gsspsth


simulate.gsspsth0 <- function(object,
                              nsim=1,
                              seed=NULL,
                              ...) {
  
  if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)                       ## initialize the RNG if necessary
  if(is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }

  ## find out the largest frequency
  idx.max <- which.max(object$freq)
  if (idx.max == 1) {
    lwr <- 1
    upr <- 2
  } else {
    freq.length <- length(object$freq)
    if (idx.max == freq.length) {
      lwr <- freq.length-1
      upr <- freq.length
    } else {
      lwr <- idx.max-1
      upr <- idx.max+1
    }
  }
  lwr <- object$mids[lwr]
  upr <- object$mids[upr]
  f.max <- optimize(object$lambdaFct,lower=lwr,upper=upr,maximum=TRUE)$objective

  time.range <- object$breaks
  
  mk.train <- function(f.max,time.range) {
    first.guess <- diff(time.range)*f.max*2
    result <- cumsum(rexp(first.guess,f.max))+time.range[1]
    result.max <- result[first.guess]
    while (result.max < time.range[2]) {
      result <- c(result,cumsum(rexp(first.guess,f.max))+result.max)
      first.guess <- length(result)
      result.max <- result[first.guess]
    }
    result <- result[result < time.range[2]]
    ratio <- object$lambdaFct(result)/f.max
    result <- result[runif(length(result))<=ratio]
    class(result) <- "spikeTrain"
    result
  }
  
  result <- lapply(1:nsim,
                   function(sIdx) {
                     ans <- lapply(1:object$nbTrials,
                                   function(tIdx) mk.train(f.max,time.range)
                                   )
                     class(ans) <- "repeatedTrain"
                     ans
                   }
                   )
  attr(result, "seed") <- RNGstate

  if (nsim==1) result <- result[[1]]

  result
}

simulate.gsspsth <- simulate.gsspsth0
