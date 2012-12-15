#'Coerce, Test and Extract from spikeTrain Objects
#'
#'\code{as.spikeTrain} attempts to coerce a numeric vector to a
#'\code{spikeTrain} object while \code{is.spikeTrain} tests if its argument is
#'such an object. \code{[.spikeTrain}, extracts a subset of a \code{spikeTrain}
#'object.
#'
#'A \code{spikeTrain} object is a \code{numeric} vector whose elements are
#'strictly increasing (that is, something which can be interpreted as a
#'sequence of times of successive events with no two events occurring at the
#'same time). The extractor method, \code{[} requires that the extracted
#'elements are without gaps, an error is returned otherwise.
#'
#'@aliases as.spikeTrain is.spikeTrain [.spikeTrain
#'@param x An object to be coerced to or to test against a \code{spikeTrain}
#'object or a \code{spikeTrain} object for \code{[}.
#'@param i indices specifying elements to extract. \emph{No gaps are allowed}.
#'@return \code{as.spikeTrain} returns a \code{spikeTrain} object or an error.
#'
#'\code{is.spikeTrain} returns \code{TRUE} if its argument is a
#'\code{spikeTrain} object and \code{FALSE} otherwise.
#'
#'\code{[} returns a \code{spikeTrain} object or an error.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{plot.spikeTrain}}, \code{\link{print.spikeTrain}},
#'\code{\link{summary.spikeTrain}}
#'@references Perkel D. H., Gerstein, G. L. and Moore G. P. (1967) Neural Spike
#'Trains and Stochastic Point Processes. I. The Single Spike Train.
#'\emph{Biophys. J.}, \bold{7}: 391-418.
#'\url{http://www.pubmedcentral.nih.gov/articlerender.fcgi?tool=pubmed&pubmedid=4292791}
#'@keywords ts survival
#'@examples
#'
#'## load CAL1S data
#'data(CAL1S)
#'## convert the data into spikeTrain objects
#'CAL1S <- lapply(CAL1S,as.spikeTrain)
#'## Are the list eleemnts now spikeTrain objects?
#'sapply(CAL1S, is.spikeTrain)
#'## look at the train of the 1st neuron
#'CAL1S[["neuron 1"]]
#'## look at the window 10-40 using the extractor function
#'CAL1S[["neuron 1"]][10 < CAL1S[["neuron 1"]] & CAL1S[["neuron 1"]] < 40] 
#'
#'
as.spikeTrain <- function(x) {
  if (inherits(x,"CountingProcessSamplePath")) x <- spikeTimes(x)
  if (!is.numeric(x)) 
    x <- as.numeric(x)
  if (!identical(length(unique(x)), length(x))) 
    stop(paste("The elements of", deparse(substitute(x)), 
               "should be all different."))
  isi <- diff(x)
  if (any(isi <= 0)) 
    stop(paste(deparse(substitute(x)), "should have strictly incresing elements."))

  class(x) <- "spikeTrain"
  x
}

is.spikeTrain <- function(x) {
  if (!("spikeTrain" %in% class(x))) return(FALSE) 
  if (!is.numeric(x)) 
    x <- as.numeric(x)
  if (!identical(length(unique(x)), length(x))) return(FALSE)
  isi <- diff(x)
  if (any(isi <= 0)) return(FALSE)
  return(TRUE)
}



#'Display Counting Process Associated with Single Spike Train
#'
#'Adds a counting process display to the classical raster plot of single spike
#'trains.
#'
#'The counting process is obtained by a call to \code{\link{stepfun}}.  When
#'\code{xlab}, \code{ylab}, \code{main}, \code{xlim} or \code{ylim} is (are)
#'missing, default values are used.
#'
#'@param x a \code{spikeTrain} object or a vector which can be coerced to such
#'an object.
#'@param xlab a character. The x label.
#'@param ylab a character. The y label.
#'@param main a character. The title.
#'@param xlim a numeric. See \code{\link{plot}}.
#'@param ylim a numeric. See \code{\link{plot}}.
#'@param do.points see \code{\link{plot.stepfun}}.
#'@param addMeanRate should the expected counting process for a Poisson process
#'with the same rate be added to the plot?
#'@param addRug should a rug representation be added at teh bottom of the plot?
#'See \code{\link{rug}}.
#'@param \dots additional arguments passed to \code{plot}, see
#'\code{\link{plot}} and \code{\link{plot.stepfun}}.
#'@return Nothing is returned, \code{plot.spikeTrain} is used for its side
#'effect, a plot is generated on the current graphic device.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{as.spikeTrain}}, \code{\link{is.spikeTrain}},
#'\code{\link{print.spikeTrain}}, \code{\link{summary.spikeTrain}},
#'\code{\link{renewalTestPlot}}, \code{\link{varianceTime}},
#'\code{\link{stepfun}}, \code{\link{plot.stepfun}}, \code{\link{rug}}
#'@references D. R. Cox and P. A. W. Lewis (1966) \emph{The Statistical
#'Analysis of Series of Events}. John Wiley and Sons.
#'
#'Brillinger, D. R. (1988) Maximum likelihood analysis of spike trains of
#'interacting nerve cells. \emph{Biol. Cybern.} \bold{59}: 189--200.
#'
#'Johnson, D.H. (1996) Point process models of single-neuron discharges.
#'\emph{J. Computational Neuroscience} \bold{3}: 275--299.
#'@keywords ts survival
#'@examples
#'
#'\dontrun{
#'data(ShallowShocks)
#'plot(as.spikeTrain(ShallowShocks$Date),
#'     xlab="Time (days)",
#'     main="Shallow Shocks Counting Process of Ogata 1988")
#'}
#'
plot.spikeTrain <- function(x,
                            xlab="Time (s)",
                            ylab="Cumulative Number of Events",
                            main=paste("Counting Process of", deparse(substitute(x))),
                            xlim=c(floor(x[1]), ceiling(x[length(x)])),
                            ylim=c(0, length(x) + 1),
                            do.points=ifelse(length(x)<100,TRUE,FALSE),
                            addMeanRate=TRUE,
                            addRug=TRUE,
                            ...) {
  
  if (!is.spikeTrain(x)) x <- as.spikeTrain(x)
  
  cp <- stepfun(x, 0:length(x))

  if (!addMeanRate) {
    plot(cp, xlab = xlab, ylab = ylab, main = main, do.points=do.points, ...)
  }
  else {
    plot(xlim, ylim, type = "n", xlab = xlab, ylab = ylab, 
         main = main, ...)
    abline(a = -xlim[1]*length(x)/diff(xlim), b = length(x)/diff(xlim), 
           col = 2)
    plot(cp,add=TRUE,do.points=do.points,...)
  }
  if (addRug) rug(x)
}



#'Print and Summary Methods for spikeTrain Objects
#'
#'Print and summary \code{\link{methods}} for \code{spikeTrain} objects.
#'
#'\code{print.spikeTrain} does in fact call the \code{plot} method for
#'\code{spikeTrain} objects.
#'
#'@aliases print.spikeTrain summary.spikeTrain
#'@param x,object A \code{spikeTrain} object.
#'@param timeUnit The unit with which the occurrence times were measured.
#'@param digits The number of digits used to print the summary (see
#'\code{\link{round}}).
#'@param \dots see \code{\link{print}} and \code{\link{summary}}.
#'@return \code{print.spikeTrain} generates a plot as a side effect.
#'
#'\code{summary.spikeTrain} returns the number of spikes, the times of the
#'first and last spikes, the mean inter-spike interval (ISI) and its sd as well
#'as the mean and sd of the log(ISI) together with the shortest and longest
#'ISIs.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{as.spikeTrain}}, \code{\link{is.spikeTrain}},
#'\code{\link{renewalTestPlot}}, \code{\link{varianceTime}},
#'\code{\link{stepfun}}
#'@keywords ts survival
#'@examples
#'
#'## load spontaneous data of 4 putative projection neurons
#'## simultaneously recorded from the cockroach (Periplaneta
#'## americana) antennal lobe
#'data(CAL1S)
#'## convert data into spikeTrain objects
#'CAL1S <- lapply(CAL1S,as.spikeTrain)
#'## look at the individual trains
#'## first the "raw" data
#'CAL1S[["neuron 1"]]
#'## next some summary information
#'summary(CAL1S[["neuron 1"]])
#'
print.spikeTrain <- function(x,...) plot.spikeTrain(x,main="Spike Train Counting Process")

summary.spikeTrain <- function(object,timeUnit="s",digits=3,...) {
  spikeTrain <- object
  rm(object)
  if (!is.spikeTrain(spikeTrain)) stop("Not a proper spike train.")
  stRange <- range(spikeTrain)
  stNb <- length(spikeTrain)
  isi <- diff(spikeTrain)
  stStat1 <- c(mean(isi),sd(isi))
  stStat2 <- c(mean(log(isi)),sd(log(isi)))
  cat(paste("A spike train with ",
            stNb,
            " events, starting at: ",
            round(stRange[1],digits=digits),
            " and ending at: ",
            round(stRange[2],digits=digits),
            " (",timeUnit,").\n",
            "The mean ISI is: ",
            round(stStat1[1],digits=digits),
            " and its SD is: ",
            round(stStat1[2],digits=digits),
            " (",timeUnit,").\n",
            "The mean log(ISI) is: ",
            round(stStat2[1],digits=digits),
            " and its SD is: ",
            round(stStat2[2],digits=digits),
            "\n",
            "The shortest interval is: ",
            round(min(isi),digits=digits),
            "\n and the longest is: ",
            round(max(isi),digits=digits),
            " (",timeUnit,").\n",
            sep=""
            )
      )
            
}



#'diff method for spikeTrain objects
#'
#'\code{diff} method for \code{spikeTrain} objects.
#'
#'
#'@param x a \code{spikeTrain} object.
#'@param \dots see \code{\link{diff}}
#'@return a \code{numeric}
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{diff}}, \code{\link{as.spikeTrain}},
#'\code{\link{is.spikeTrain}}
#'@keywords ts
#'@examples
#'
#'data(CAL1S)
#'## convert data into spikeTrain objects
#'CAL1S <- lapply(CAL1S,as.spikeTrain)
#'## look at the individual trains
#'## first the "raw" data
#'CAL1S[["neuron 1"]]
#'## get the isi of neuron 1
#'n1.isi <- diff(CAL1S[["neuron 1"]])
#'
diff.spikeTrain <- function(x,...) {
  class(x) <- NULL
  diff(x,...)
}

"[.spikeTrain" <- function(x,i) {
  xClass <- class(x)
  x <- unclass(x)
  newI <- seq(x)[i]
  if (all(diff(newI)==1)) {
    x <- x[newI]
    class(x) <- xClass
    return(x)
  }
  stop("Wrong i specification.")
}
