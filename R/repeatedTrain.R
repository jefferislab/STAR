#'Coerce and Test repeatedTrain Objects
#'
#'\code{as.repeatedTrain} attempts to coerce a list with numeric vector
#'elements to a \code{repeatedTrain} object while \code{is.repeatedTrain} tests
#'if its argument is such an object.
#'
#'A \code{repeatedTrain} object is list of \code{spikeTrain} objects. It is
#'used to store the responses of a given neuron to repeated stimulations.
#'
#'@aliases as.repeatedTrain is.repeatedTrain
#'@param x An object to be coerced to or to test against a \code{repeatedTrain}
#'object.
#'@return \code{as.repeatedTrain} returns a \code{repeatedTrain} object or an
#'error.
#'
#'\code{is.repeatedTrain} returns \code{TRUE} if its argument is a
#'\code{repeatedTrain} object and \code{FALSE} otherwise.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{plot.repeatedTrain}}, \code{\link{print.repeatedTrain}},
#'\code{\link{summary.repeatedTrain}}, \code{\link{psth}},
#'\code{\link{raster}}, \code{\link{as.spikeTrain}},
#'\code{\link{is.spikeTrain}}
#'@keywords ts survival
#'@examples
#'
#'## load CAL1V data
#'data(CAL1V)
#'## convert them to repeatedTrain objects
#'CAL1V <- lapply(CAL1V, as.repeatedTrain)
#'## did the conversion work?
#'sapply(CAL1V, is.repeatedTrain)
#'## look at the raster of the 1st neuron
#'CAL1V[["neuron 1"]]
#'
as.repeatedTrain <- function(x) {
  if (!is.list(x)) 
    x <- list(x)
  if (!is.null(names(x))) {
    xN <- names(x)
  } else {
    xN <- paste("stim.",1:length(x))
  }
  
  x <- lapply(x, as.spikeTrain)
  names(x) <- xN
  class(x) <- "repeatedTrain"
  x
}

is.repeatedTrain <- function(x) {
  if (!("repeatedTrain" %in% class(x))) return(FALSE)
  if (!all(sapply(x,is.spikeTrain))) return(FALSE)
  return(TRUE)
}



#'Generate a Raster Plot
#'
#'Given a list of spike trains (or a \code{repeatedTrain} object) where each
#'train was acquired during, say, one presentation of a given stimulus, a
#'raster plot is generated. If stimulus time properties are specified, the
#'stimulus application time also appears on the plot.
#'
#'Basic raster plot stuff.
#'
#'@aliases raster plot.repeatedTrain
#'@param x a \code{repeatedTrain} object or a list which can be coerced to such
#'an object.
#'@param stimTimeCourse \code{NULL} (default) or a two elements vector
#'specifying the time boundaries (in s) of a stimulus presentation.
#'@param colStim the background color used for the stimulus.
#'@param xlim a numeric (default value supplied). See \code{\link{plot}}.
#'@param pch data symbol used for the spikes. See \code{\link{plot}}.
#'@param xlab a character (default value supplied). See \code{\link{plot}}.
#'@param ylab a character (default value supplied). See \code{\link{plot}}.
#'@param main a character (default value supplied). See \code{\link{plot}}.
#'@param \dots see \code{\link{plot}}.
#'@return Nothing is returned \code{raster} is used for its side effect, a plot
#'is generated on the current graphical device.
#'@note Brillinger (1992) calls these plots "rastor" instead of raster...
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{as.repeatedTrain}}, \code{\link{is.repeatedTrain}},
#'\code{\link{print.repeatedTrain}}, \code{\link{summary.repeatedTrain}},
#'\code{\link{psth}}
#'@references Brillinger, David R. (1992) Nerve Cell Spike Train Data Analysis:
#'A Progression of Technique. \emph{JASA} \bold{87}: 260--271.
#'@keywords ts survival
#'@examples
#'
#'## Load Vanillin responses data (first cockroach data set)
#'data(CAL1V)
#'## convert them into repeatedTrain objects
#'## The stimulus command is on between 4.49 s and 4.99s
#'CAL1V <- lapply(CAL1V,as.repeatedTrain)
#'## look at the individual raster plots
#'raster(CAL1V[["neuron 1"]],stimTimeCourse=c(4.49,4.99),main="N1")
#'plot(CAL1V[["neuron 2"]],stimTimeCourse=c(4.49,4.99),main="N2")
#'plot(CAL1V[["neuron 3"]],stimTimeCourse=c(4.49,4.99),main="N3")
#'plot(CAL1V[["neuron 4"]],stimTimeCourse=c(4.49,4.99),main="N4")
#'
raster <- function(x,
                   stimTimeCourse = NULL,
                   colStim = "grey80", 
                   xlim, pch, xlab, ylab, main,
                   ...) {

  do.call(plot.repeatedTrain,as.list(match.call()[-1]))

}

plot.repeatedTrain <- function (x,
                                stimTimeCourse = NULL,
                                colStim = "grey80", 
                                xlim, pch, xlab, ylab, main,
                                ...) {
  if (!is.repeatedTrain(x)) x <- as.repeatedTrain(x)
  nbTrains <- length(x)
  if (missing(xlim)) 
    xlim <- c(0, ceiling(max(sapply(x, max))))
  if (missing(xlab)) 
    xlab <- "Time (s)"
  if (missing(ylab)) 
    ylab <- "trial"
  if (missing(main)) 
    main <- paste(deparse(substitute(x)), "raster")
  if (missing(pch)) 
    pch <- ifelse(nbTrains <= 20, "|", ".")
  oldpar <- par(mar = c(5, 4, 2, 1))
  on.exit(par(oldpar))
  acquisitionDuration <- max(xlim)
  plot(c(0, acquisitionDuration), c(0, nbTrains + 1), type = "n", 
       xlab = xlab, ylab = ylab, xlim = xlim, ylim = c(1, nbTrains + 
                                                1), bty = "n", main = main, ...)
  if (!is.null(stimTimeCourse)) {
    rect(stimTimeCourse[1], 0.1, stimTimeCourse[2], nbTrains + 
         0.9, col = colStim, lty = 0)
  }
  invisible(sapply(1:nbTrains, function(idx) points(x[[idx]], 
                                                    numeric(length(x[[idx]])) + idx, pch = pch)))
  axis(2, at = 1:nbTrains)
  axis(1)

}

plot.psth <- function(x,
                      stimTimeCourse = NULL,
                      colStim = "grey80", 
                      colCI = NULL,
                      xlab, ylab, main, xlim, 
                      ylim, lwd = 2, col = 1,
                      ...) {
  if (!is.null(stimTimeCourse)) {
    if (length(stimTimeCourse) != 2) 
      stop(paste(deparse(substitute(stimTimeCourse)),
                 "should be a vector with 2 elements.")
           )
  }
  if (missing(xlab)) 
    xlab <- "Time (s)"
  if (missing(ylab)) 
    ylab <- "Freq (Hz)"
  if (missing(main)) {
    nameList <- deparse(x$call[["repeatedTrain"]])
    main <- paste(nameList, "PSTH")
  }

  withCI <- !is.null(x$ciUp)

  originalBreaksArg <- x$call[["breaks"]]
  if (!is.null(originalBreaksArg) && (length(eval(originalBreaksArg))==2))
    smoothOne <- TRUE
  else
    smoothOne <- FALSE
  if (missing(xlim)) {
    if (!smoothOne) xlim <- range(x$breaks)
    else xlim <- range(x$mids) + c(-0.5,0.5)*x$breaks["step"]
  }  
  if (missing(ylim))
    ylim <- c(0, ifelse(withCI, max(x$ciUp), max(x$freq))) ## Bug fix thanks to Gregory Jefferis

  

  plot(xlim,
       ylim,
       type = "n", 
       xlab = xlab,
       ylab = ylab,
       main = main,
       xlim = xlim, 
       ylim = ylim,
       ...)
  
  if (!is.null(stimTimeCourse)) {
    rect(stimTimeCourse[1], 0, stimTimeCourse[2], ylim[2], 
         col = colStim, lty = 0)
  }
  if (withCI) {
    if (is.null(colCI)) {
      if (!smoothOne) {
        sapply(1:length(x$mids),
               function(idx)
               segments(x$breaks[idx], x$ciLow[idx],
                        x$breaks[idx + 1], x$ciLow[idx], 
                        lty = 2)
               )
        sapply(1:length(x$mids),
               function(idx)
               segments(x$breaks[idx], x$ciUp[idx],
                        x$breaks[idx + 1], x$ciUp[idx], 
                        lty = 2)
               )
      } else {
        lines(x$mids, x$ciLow, lty = 2)
        lines(x$mids, x$ciUp, lty = 2)
      }
    } else {
      if (!smoothOne) {
        sapply(1:length(x$mids),
               function(idx)
               rect(x$breaks[idx], x$ciLow[idx],
                    x$breaks[idx + 1], x$ciUp[idx], 
                    lty = 0, col = colCI))
      } else {
        polygon(c(x$mids, rev(x$mids)), c(x$ciLow, rev(x$ciUp)), 
                col = colCI, border = NA)
      }
    }
  }
  if (!smoothOne) {
    invisible(sapply(1:length(x$mids),
                     function(idx) {
                       segments(x$breaks[idx], x$freq[idx],
                                x$breaks[idx + 1], x$freq[idx],
                                lwd = lwd, col = col)
                       if (idx == 1)
                         segments(x$breaks[idx], 0, x$breaks[idx], 
                                  x$freq[idx], lwd = lwd, col = col)
                       if (idx != length(x$mids)) {
                         segments(x$breaks[idx + 1],
                                  x$freq[idx],
                                  x$breaks[idx + 1],
                                  x$freq[idx + 1],
                                  lwd = lwd,
                                  col = col)
                       } else {
                         segments(x$breaks[idx + 1],
                                  x$freq[idx],
                                  x$breaks[idx + 1],
                                  0, lwd = lwd, col = col)
                       }
                     }
                     )
              )
  } else {
    lines(x$mids, x$freq, col = col, lwd = lwd)
  }
  abline(h = 0)
}



#'Compute and Plot Peri-Stimulus Time Histogram
#'
#'\code{psth} computes and \code{plot.psth} plots a peri-stimulus time
#'histogram (called PST, post-stimulus time histogram by Gerstein and Kiang
#'(1960)) from repeated presentations of a stimulation. Confidence bands can be
#'obtained using the Poisson approximation.
#'
#'When confidence bands are requested they are obtained from the qunatiles of
#'the \code{\link{Poisson}} distribution.
#'
#'When a 2 elements vector is used as \code{breaks} argument it is interpreted
#'as specifying a bin width (first element if elements are unnamed, \code{"bw"}
#'element otherwise) and a step (second element if elements are unnamed,
#'\code{"step"} element otherwise). The idea is then to obtain a smoother
#'looking PSTH by counting spikes within overlapping bins. That is if the
#'center of the ith bin is xi the one of the (i+1)th bin will be xi + step.
#'
#'@aliases psth plot.psth
#'@param repeatedTrain a \code{repeatedTrain} object or a list which can be
#'coerced to such an object.
#'@param x a \code{psth} object.
#'@param stimTimeCourse \code{NULL} (default) or a two elements vector
#'specifying the time boundaries (in s) of a stimulus presentation.
#'@param colStim the background color used for the stimulus.
#'@param breaks a numeric. A single number is interpreted has the number of
#'bins; a vector of length 2 is interpreted as the bin width and the step to
#'use (see details); otherwise interpreted as the position of the "breaks"
#'between bins.
#'@param include.lowest corresponding argument of \code{\link{hist}}.
#'@param right corresponding argument of \code{\link{hist}}.
#'@param plot corresponding argument of \code{\link{hist}}.
#'@param CI The coverage probability of the confidence intervals.
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
#'@param \dots see \code{\link{plot}}.
#'@return When \code{plot} is set to \code{FALSE} in \code{psth}, a list of
#'class \code{psth} is returned and no plot is generated. This list has the
#'following components:
#'
#'When \code{plot} is set to \code{TRUE} nothing is returned and a plot is
#'generated as a side effect. Of course the same occurs upon calling
#'\code{plot.psth} with a \code{psth} object argument.
#'@returnItem freq a vector containing the instantaneous firing rate.
#'@returnItem ciUp a vector with the upper limit of the confidence band.
#'@returnItem ciLow a vector with the lower limit of the confidence band.
#'@returnItem breaks a numeric vector with the breaks in between which spikes
#'were counted. Similar to the component of the same name returned by
#'\code{\link{hist}}.
#'@returnItem mids a numeric vector with the mid points of \code{breaks}.
#'Similar to the component of the same name returned by \code{\link{hist}}.
#'@returnItem counts a matrix with as many rows as components in
#'\code{repeatedTrain} and as many columns as bins. Each element of the matrix
#'contains the number of spikes falling in a given trial in a given bin.
#'@returnItem nbTrials the number of stimulations.
#'@returnItem call the matched call.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{as.repeatedTrain}}, \code{\link{is.repeatedTrain}},
#'\code{\link{print.repeatedTrain}}, \code{\link{summary.repeatedTrain}},
#'\code{\link{raster}}
#'@references Gerstein, George L. and Kiang, Nelson Y.-S. (1960) An Approach to
#'the Quantitative Analysis of Electrophysiological Data from Single Neurons.
#'\emph{Biophysical Journal} \bold{1}: 15--28.
#'\url{http://www.pubmedcentral.nih.gov/articlerender.fcgi?tool=pubmed&pubmedid=13704760}
#'
#'Kalbfleisch, J. G. (1985) \emph{Probability and Statistical Inference. Volume
#'2: Statistical Inference}. Springer-Verlag.
#'@keywords ts survival
#'@examples
#'
#'## Load Vanillin responses data (first cockroach data set)
#'data(CAL1V)
#'## convert them into repeatedTrain objects
#'## The stimulus command is on between 4.49 s and 4.99s
#'CAL1V <- lapply(CAL1V,as.repeatedTrain)
#'## look at the individual raster plots
#'plot(CAL1V[["neuron 1"]],stimTimeCourse=c(4.49,4.99),main="N1")
#'## Create a simple black and white PSTH for neuron 1
#'psth(CAL1V[["neuron 1"]],stimTimeCourse=c(4.49,4.99),breaks=20)
#'## Rebuilt the same PSTH but with red confidence bands
#'psth(CAL1V[["neuron 1"]],stimTimeCourse=c(4.49,4.99),breaks=20,colCI=2)
#'## Make the PSTH smoother
#'psth(CAL1V[["neuron 1"]],stimTimeCourse=c(4.49,4.99),breaks=c(bw=0.5,step=0.05),colCI=2)
#'## Make a plot with PSTHs from 4 neurons superposed
#'## First get lists containing PSTHs from each neuron
#'psth1 <- psth(CAL1V[["neuron 1"]],breaks=c(bw=0.5,step=0.05),plot=FALSE)
#'psth2 <- psth(CAL1V[["neuron 2"]],breaks=c(bw=1,step=0.1),plot=FALSE)
#'psth3 <- psth(CAL1V[["neuron 3"]],breaks=c(bw=0.5,step=0.05),plot=FALSE)
#'psth4 <- psth(CAL1V[["neuron 4"]],breaks=c(bw=2,step=0.2),plot=FALSE)
#'## Get the maximal frequency to display
#'maxFreq <- max(max(psth1$ciUp),max(psth2$ciUp),max(psth3$ciUp),max(psth4$ciUp))
#'## Build plot
#'plot(c(0,10),c(0,75),type="n",
#'     xaxs="i",yaxs="i",xlab="Time (s)",
#'     ylab="Freq. (Hz)",
#'     main="PSTHs from 4 simultaneously recorded neurons",
#'     sub="20 stimulations with vanillin were used.")
#'## Add rectangle corresponding to stimulation command
#'rect(4.49,0,4.99,75,col="grey80",lty=0)
#'## Add the neurons PSTHs as confidence bands
#'polygon(c(psth1$mids,rev(psth1$mids)),c(psth1$ciLow,rev(psth1$ciUp)),col=1,border=NA)
#'polygon(c(psth2$mids,rev(psth2$mids)),c(psth2$ciLow,rev(psth2$ciUp)),col=2,border=NA)
#'polygon(c(psth3$mids,rev(psth3$mids)),c(psth3$ciLow,rev(psth3$ciUp)),col=3,border=NA)
#'polygon(c(psth4$mids,rev(psth4$mids)),c(psth4$ciLow,rev(psth4$ciUp)),col=4,border=NA)
#'legend(0.1,maxFreq,legend=paste("neuron",1:4),lty=1,col=1:4,bty="n")
#'
psth <- function (repeatedTrain,
                  breaks=20,
                  include.lowest = TRUE,
                  right = TRUE, 
                  plot = TRUE,
                  CI = 0.95,
                  ...
                  ) {
  
  if (!is.repeatedTrain(repeatedTrain))
    repeatedTrain <- as.repeatedTrain(repeatedTrain)

  if (!inherits(breaks,"numeric"))
    stop("breaks should be a numeric.")
 
  nbTrials <- length(repeatedTrain)
  breaksName <- deparse(substitute(breaks))
  l <- floor(min(sapply(repeatedTrain,function(l) l[1])))
  r <- ceiling(max(sapply(repeatedTrain,function(l) l[length(l)])))
      
  if (length(breaks) != 2) {
    if (length(breaks)==1) breaks <- seq(l,r,length.out=breaks+1)
    counts <- t(sapply(repeatedTrain,
                       function(train)
                       hist(x = unclass(train), breaks = breaks,
                            include.lowest = include.lowest, 
                            right = right, plot = FALSE)$counts
                       )
                )
    h <- list(breaks = breaks,
              counts = counts,
              mids = breaks[-length(breaks)]+diff(breaks)/2)
    f <- colSums(counts)/(nbTrials * diff(h$breaks))
  } else {
    if (is.null(names(breaks))) {
      bw <- breaks[1]
      step <- breaks[2]
    } else {
      if (!all(names(breaks) %in% c("bw", "step"))) 
        stop(paste(breaksName, "should have named elements: bw and step"))
      bw <- as.numeric(breaks["bw"])
      step <- as.numeric(breaks["step"])
    } ## End of conditional on is.null(names(breaks)
    bwh <- bw/2
    breaks <- c(bw = bw, step = step)
    mids <- seq(bwh, r - bwh, by = step)
    counts <- t(sapply(repeatedTrain,
                       function(train)
                       sapply(mids,
                              function(m)
                              ifelse(right,
                                     sum(m - bwh < train & train <= m + bwh),
                                     sum(m - bwh <= train & train < m + bwh))
                              )
                       )
                )
    h <- list(breaks = breaks,
              counts = counts,
              mids = mids)
    f <- colSums(counts)/(nbTrials * bw)
  } ## End of conditional on length(breaks) != 2
  if (!is.null(CI)) {
    expectedCount <- colSums(counts)
    ciUp <- sapply(1:length(h$mids),
                   function(idx) qpois(1-(1-CI)/2,expectedCount[idx])
                   )
    if (length(breaks) > 2) ciUp <- ciUp / (nbTrials * diff(h$breaks))
    else ciUp <- ciUp / (nbTrials * bw)
    ciLow <- sapply(1:length(h$mids),
                    function(idx) qpois((1-CI)/2,expectedCount[idx])
                    )
    if (length(breaks) > 2) ciLow <- ciLow / (nbTrials * diff(h$breaks))
    else ciLow <- ciLow / (nbTrials * bw)
  }

  if (!is.null(CI)) {
    result <- list(freq = f, ciUp = ciUp, ciLow = ciLow, 
                   breaks = h$breaks, mids = h$mids, counts = h$counts, 
                   nbTrials = nbTrials, call = match.call())
    class(result) <- "psth"
  } else {
    result <- list(freq = f, ciUp = NULL, ciLow = NULL, 
                   breaks = h$breaks, mids = h$mids, counts = h$counts, 
                   nbTrials = nbTrials, call = match.call())
    class(result) <- "psth"
  }
  
  if (plot) {
    plot(result,...)
  } else {
    return(result)
  }
  
}



#'Print and Summary Methods for repeatedTrain Objects
#'
#'Print and summary \code{\link{methods}} for \code{repeatedTrain} objects.
#'
#'\code{print.repeatedTrain} calls \code{\link{plot.repeatedTrain}}
#'
#'@aliases print.repeatedTrain summary.repeatedTrain
#'print.summary.repeatedTrain
#'@param x a \code{repeatedTrain} or a \code{summary.repeatedTrain} object.
#'@param object a \code{repeatedTrain} object
#'@param responseWindow a 2 elements vector specifying the begining and the end
#'of the neuron response.
#'@param acquisitionWindow a 2 elements vector specifying the begining and the
#'end of the acquisition. If \code{missing} values are obtained using the
#'\code{\link{floor}} of the smallest spike time and the \code{\link{ceiling}}
#'of the largest one.
#'@param \dots additional arguments passed to function \code{\link{chisq.test}}
#'or \code{\link{print}}.
#'@return \code{summary.repeatedTrain} returns a LIST of class
#'\code{summary.repeatedTrain} with the following components:
#'@returnItem nbRepeates The number of repetitions.
#'@returnItem acquisitionWindow The acquisition window.
#'@returnItem stats A matrix with as many rows as repetitions. The first column
#'contains the total number of spikes generated by the neuron during a given
#'repeat (this column appears under the heading "nb" when the object is
#'printed). The second column contains the corresponding average discharge rate
#'(this column appears under the heading "nu" when the object is printed). If a
#'\code{responseWindow} was specified, the third column contains the number of
#'spikes generated by the neuron during the response period and the fourth
#'column contains the corresponding rate (these column appear under the
#'headings "nbR" and "nuR", respectively when the object is printed).
#'@returnItem globalPval The p value of the chi square test for homogeneity of
#'the total number of spikes generated accross repetitions. Thats a rough
#'stationarity test.
#'@returnItem responsePval If a \code{responseWindow} was specified, the p
#'value of the chi square test for homogeneity of the number of spikes
#'generated within the "response window" accross repetitions.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{as.repeatedTrain}}, \code{\link{is.repeatedTrain}},
#'\code{\link{plot.repeatedTrain}}, \code{\link{raster}}, \code{\link{psth}}
#'@keywords ts survival
#'@examples
#'
#'## Load the Vanillin responses of the first
#'## cockroach data set
#'data(CAL1V)
#'## convert them into repeatedTrain objects
#'## The stimulus command is on between 4.49 s and 4.99s
#'CAL1V <- lapply(CAL1V,as.repeatedTrain)
#'## Generate raster plot for the neurons
#'raster(CAL1V[["neuron 1"]],c(4.49,4.99))
#'plot(CAL1V[["neuron 2"]],c(4.49,4.99))
#'plot(CAL1V[["neuron 3"]],c(4.49,4.99))
#'## Basic summary of neuron 1
#'summary(CAL1V[["neuron 1"]])
#'## Enhanced summary giving a response window between 5 and 5.5s
#'summary(CAL1V[["neuron 1"]],c(5,5.5))
#'
print.repeatedTrain <- function(x,...) {
  main <- "Raster plot"
  plot(x,main=main,...)
}

summary.repeatedTrain <- function(object,
                                  responseWindow,
                                  acquisitionWindow,
                                  ...) {
  repeatedTrain <- object
  rm(object)
  nbRepeates <- length(repeatedTrain)
  if (missing(acquisitionWindow))
    acquisitionWindow <- c(min(sapply(repeatedTrain,function(l) floor(l[1]))),
                           max(sapply(repeatedTrain,function(l) ceiling(l[length(l)])))
                           )
  
  acquisitionDuration <- diff(acquisitionWindow)
  if (!missing(responseWindow)) {
    if (length(responseWindow) != 2) stop("Wrong responseWindow.")
    if ((responseWindow[1] < acquisitionWindow[1]) ||
        (responseWindow[2] > acquisitionWindow[2]) )
      warning("responseWindow does not make much sense.")
    stimDuration <- diff(responseWindow)
  } else {
    responseWindow <- NULL
    stimDuration <- NULL
  }
  theStats <- sapply(repeatedTrain,
                     function(l) {
                       nb <- length(l)
                       nu <- nb/acquisitionDuration
                       result <- c(nb=nb,nu=nu)
                       if (!is.null(stimDuration)) {
                         goodOnes <- responseWindow[1] < l & l < responseWindow[2]
                         nbR <- sum(goodOnes)
                         nuR <- nbR/stimDuration
                         result <- c(result,c(nbR=nbR,nuR=nuR))
                       }
                       result
                     }
                     )
  stationarityTest <- chisq.test(theStats["nb",],...)$p.value
  if (!is.null(stimDuration)) {
    stationarityTestR <- chisq.test(theStats["nbR",],...)$p.value
  } else {
    stationarityTestR <- NULL
  }

  result <- list(nbRepeates=nbRepeates,
                 acquisitionWindow=acquisitionWindow,
                 responseWindow=responseWindow,
                 stats=theStats,
                 globalPval=stationarityTest,
                 responsePval=stationarityTestR)
  class(result) <- "summary.repeatedTrain"
  return(result)
}

print.summary.repeatedTrain <- function(x,...) {

  obj <- x
  rm(x)
  cat(paste(obj$nbRepeates, " repeats in the object.\n"))
  print(t(obj$stats))
  cat(paste("The p value of the chi 2 test for the stationarity accross repeats is:\n",
            obj$globalPval,
            ".\n",sep="")
      )
  if (!is.null(obj$responseWindow))
    cat(paste("The p value of the chi 2 test for the stationarity accross repeats during the stim. is:\n",
            obj$responsePval,
            ".\n",sep="")
      )
  
}




#'Generates a Data Frame from a repeatedTrain Object After Time Binning
#'
#'Generates a \code{\link{data.frame}} object out of a \code{repeatedTrain}
#'object after time binning in order to study trials stationarity with a
#'\code{\link{glm}} fit.
#'
#'The bins are placed between the \code{\link{floor}} of the smallest spike
#'time and the \code{\link{ceiling}} of the largest one when \code{breaks} is a
#'scalar. After time binning the number of spikes of each trial falling in each
#'bin is counted (in the same way as the \code{counts} component of a
#'\code{\link{psth}} list is obtained). This matrix of count is then formatted
#'as a data frame.
#'
#'@param repeatedTrain a \code{repeatedTrain} object or a list which can be
#'coerced to such an object.
#'@param breaks a numeric. A single number is interpreted has the number of
#'bins; a vector is interpreted as the position of the "breaks" between bins.
#'@return A \code{\link{data.frame}} with the following variables:
#'@returnItem Count a count (number of spikes in a given bin at a given trial).
#'@returnItem Bin the bin index (a \code{\link{factor}}.
#'@returnItem Trial the trial index (a \code{\link{factor}}.
#'@returnItem Rate the count divided by the length of the corresponding bin.
#'@returnItem Time the time of the midpoints of the bins.
#'@note When a \code{\link{glm}} of the poisson family is used for subsequent
#'analysis the important implicit hypothesis of an inhomogenous Poisson train
#'is of course made.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{as.repeatedTrain}}, \code{\link{psth}}
#'@keywords ts
#'@examples
#'
#'\dontrun{
#'## Load the Vanillin responses of the first
#'## cockroach data set
#'data(CAL1V)
#'## convert them into repeatedTrain objects
#'## The stimulus command is on between 4.49 s and 4.99s
#'CAL1V <- lapply(CAL1V,as.repeatedTrain)
#'## Generate raster plot for neuron 1
#'raster(CAL1V[["neuron 1"]],c(4.49,4.99))
#'## make a smooth PSTH of these data
#'psth(CAL1V[["neuron 1"]],stimTimeCourse=c(4.49,4.99),breaks=c(bw=0.5,step=0.05),colCI=2,xlim=c(0,10))
#'## add a grid to the plot
#'grid()
#'## The response starts after 4.5 s and is mostly over after 6 s: create
#'## breaks accordingly
#'myBreaks <- c(0,2.25,4.5,seq(4.75,6.25,0.25),seq(6.5,11,0.5))
#'## get a count data frame
#'CAL1Vn1DF <- df4counts(CAL1V[["neuron 1"]],myBreaks)
#'## use a box plot to look at the result
#'boxplot(Rate ~ Time, data=CAL1Vn1DF)
#'## watch out here the time scale is distorted because of our
#'## choice of unequal bins
#'## Fit a glm of the Poisson family taking both Bin and Trial effects
#'CAL1Vn1DFglm <- glm(Count ~ Bin + Trial,family=poisson,data=CAL1Vn1DF)
#'## use an anova to see that both the Bin effect and the trial effect are
#'## highly significant
#'anova(CAL1Vn1DFglm, test="Chisq")
#'}
df4counts <- function(repeatedTrain,
                      breaks=length(repeatedTrain)
                      ) {
  
  if (!is.repeatedTrain(repeatedTrain))
    repeatedTrain <- as.repeatedTrain(repeatedTrain)

  nbTrials <- length(repeatedTrain)
  breaksName <- deparse(substitute(breaks))
  trainNames <- names(repeatedTrain)
  if (is.null(trainNames))
    trainNames <- paste("trial",1:nbTrials)
  
  l <- floor(min(sapply(repeatedTrain,function(l) l[1])))
  r <- ceiling(max(sapply(repeatedTrain,function(l) l[length(l)])))

  if (length(breaks) == 1) {
    ## breaks is interpreted as a number of bins
    breaks <- seq(l,r,length.out=breaks+1)
  } ## End of conditional on length(breaks) == 1

  counts <- sapply(repeatedTrain,
                   function(train)
                   hist(x=train[l <= train & train <= r],
                        breaks=breaks, plot=FALSE)$counts
                   )

  mids <- breaks[-length(breaks)]+diff(breaks)/2
  nb <- length(mids)
  data.frame(Count=as.numeric(counts),
             Bin=factor(rep(1:nb,nbTrials)),
             Trial=factor(rep(trainNames,each=nb)),
             Rate=as.numeric(counts)/diff(breaks),
             Time=rep(mids,nbTrials)
             )

}
