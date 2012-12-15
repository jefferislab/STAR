#'Evaluates and Plots a Spike Train Firing Rate's Evolution
#'
#'\code{rateEvolution} evaluates and \code{plot.rateEvolution} plots the firing
#'rate evolution of a \code{spikeTrain} object. The evaluation is done by
#'convolving the spike train with a kernel like in \code{density} estimation.
#'
#'\code{rateEvolution} is mainly a wrapper for \code{\link{density}} which also
#'adjusts the result of the latter such that the y component of the returned
#'list is an instantaneous firing rate. If the length of \code{x} is smaller or
#'equal to 1 and if \code{from} or \code{to} is (are) \code{missing} the
#'returned object has then each of its components set to \code{NA} except
#'\code{data.name} (see below). If the length of \code{x} is smaller or equal
#'to 1 and if both \code{from} and \code{to} are specified a \code{missing}
#'\code{bw} is then set to 3 times the spacing between the points of the
#'regular grid on which the density is evaluated.
#'
#'\code{plot.rateEvolution} is also a wrapper for \code{\link{plot.density}}
#'which only adjust the default value of some arguments.
#'
#'@aliases rateEvolution plot.rateEvolution
#'@param x a \code{spikeTrain} object or an object which can be coerced to it
#'for \code{rateEvolution} or a \code{rateEvolution} object for
#'\code{plot.rateEvolution}.
#'@param bw the kernel bin width in seconds. If \code{missing} it is set to 10
#'times the median inter-spike interval of \code{x}.
#'@param kernel see \code{\link{density}}.
#'@param n see \code{\link{density}}.
#'@param from see \code{\link{density}}.
#'@param to see \code{\link{density}}.
#'@param na.rm see \code{\link{density}}.
#'@param main see \code{plot.density}.
#'@param xlab see \code{plot.density}.
#'@param ylab see \code{plot.density}.
#'@param type see \code{plot.density}.
#'@param zero.line see \code{plot.density}.
#'@param \dots see \code{\link{density}} and \code{plot.density}.
#'@return \code{rateEvolution} returns a LIST of class \code{rateEvolution}
#'which inherits from class \code{density}.
#'
#'\code{plot.rateEvolution} is called for its side effect: a plot is generated.
#'@returnItem x the \code{n} coordinates of the points where the density is
#'estimated. See \code{\link{density}}.
#'@returnItem y the estimated rate (in 1/s).  These will be non-negative, but
#'can be zero.
#'@returnItem bw the bandwidth used.
#'@returnItem n the sample size after elimination of missing values.
#'@returnItem call the call which produced the result.
#'@returnItem data.name the deparsed name of the \code{x} argument.
#'@returnItem has.na logical, for compatibility (always \code{FALSE}).
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{as.spikeTrain}}, \code{\link{density}},
#'\code{\link{plot.density}}, \code{\link{mkREdf}}
#'@keywords ts
#'@examples
#'
#'## load Purkinje cell data recorded in cell-attached mode
#'data(sPK)
#'## coerce sPK to a spikeTrain object
#'sPK <- lapply(sPK, as.spikeTrain)
#'## get the rate evolution in ctl condition
#'sPKreCTL <- rateEvolution(sPK[["ctl"]])
#'## plot the result
#'plot(sPKreCTL)
#'## check the bin width which was actually used
#'sPKreCTL$bw
#'## look at the effect of a 10 times larger bw
#'plot(rateEvolution(sPK[["ctl"]],bw=10*sPKreCTL$bw))
#'## look at the effect of a 10 times smaller one
#'plot(rateEvolution(sPK[["ctl"]],bw=sPKreCTL$bw/10))
#'## get the rate evolution in bicuculline conditions
#'sPKreBICU <- rateEvolution(sPK[["bicu"]])
#'## plot results
#'plot(sPKreBICU,col=2)
#'## add the ctl rate evolution
#'lines(sPKreCTL)
#'
rateEvolution <- function(x,
                          bw,
                          kernel = c("gaussian", 
                            "epanechnikov", "rectangular", "triangular", "biweight", 
                            "cosine", "optcosine"), 
                          n = 512,
                          from, to,
                          na.rm = FALSE, 
                          ...) {

  data.name <- deparse(substitute(x))
  ## check if x is a spikeTrain object
  if (!is.spikeTrain(x)) x <- as.spikeTrain(x)

  if (length(x) > 1) {
    if (missing(bw)) bw <- 10*median(diff(x))

    class(x) <- NULL
    result <- density(x,bw,kernel=kernel,n=n,from=from,to=to,na.rm=na.rm,...)
    result$y <- result$y*result$n
    result$data.name <- data.name
  } else {
    if (!missing(from) && !missing(to)) {
      xx <- seq(from,to,length.out=n)
      if (missing(bw)) bw <- 3*diff(xx)[1]
      class(x) <- NULL
      result <- density(x,bw,kernel=kernel,n=n,from=from,to=to,na.rm=na.rm,...)
      result$y <- result$y*result$n
      result$data.name <- data.name
    } else {
      result <- list(x=NA,y=NA,bw=NA,n=NA,call=NA,data.name,has.na=NA)
    }
  }
  
  class(result) <- c("rateEvolution","density")
  result
}

plot.rateEvolution <- function(x,
                               main = NULL,
                               xlab = NULL,
                               ylab = "Rate (Hz)",
                               type = "l", 
                               zero.line = TRUE,
                               ...) {
  
  if (is.null(main)) main <- paste(x$data.name,"Rate Evolution")
  if (is.null(xlab)) xlab <- "Time (s)"
  NextMethod("plot",,main=main,xlab=xlab,ylab=ylab,type=type,zero.line=zero.line,...)
}




#'Evaluates RateEvolutions for spikeTrain Lists and Returns Data Frame
#'
#'Given a list of \code{spikeTrain} or \code{repeatedTrain} objects
#'\code{mkREdf} evaluates the rate evolution of each train and returns a data
#'frame suitable for use with \code{coplot}, \code{xyplot} and \code{qplot}.
#'
#'\code{mkREdf} calls \code{\link{rateEvolution}} on every \code{spikeTrain} in
#'\code{x}. If \code{from} and \code{to} are missing, they are internally set
#'to the \code{floor} of the global minimal spike time contained in \code{x}
#'and to the \code{ceiling} of the global maximal time.
#'
#'@param x a \emph{named} list of \code{spikeTrain} or \code{repeatedTrain}
#'objects.
#'@param longitudinal a \code{character} vector with the names of the different
#'"conditions" applied to each neuron like "ctl", "bicu" or "stim. 1", "stim.
#'2", ..., "stim. 20". Default provided.
#'@param across a \code{character} vector with the names of the different
#'neurons. Default provided.
#'@param bw see \code{\link{rateEvolution}}. This can be a vector.
#'@param kernel see \code{\link{rateEvolution}}.
#'@param n see \code{\link{rateEvolution}}.
#'@param from see \code{\link{rateEvolution}}.
#'@param to see \code{\link{rateEvolution}}.
#'@param na.rm see \code{\link{rateEvolution}}.
#'@param minusMean should the mean of the rate evolution along the across
#'"dimension" be subtracted from each individual rate evolution along this
#'dimension?
#'@return A data frame with the following variables:
#'@returnItem time The time (in s) at which the rate was evaluated.
#'@returnItem rate The rate (in 1/s).
#'@returnItem longitudinal A factor corresponding to the argument with the same
#'name.
#'@returnItem across A factor corresponding to the argument with the same name.
#'@note argument \code{minusMean} is now here as an "experimental" feature. The
#'idea is that it could be used to detect non-stationarities of the reponses
#'(in a repeated stimulation context) which would be correlated across
#'different neurons. I'm not sure yet if this will be useful or not.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{as.spikeTrain}}, \code{\link{as.repeatedTrain}},
#'\code{\link{data.frame}}, \code{\link{factor}}, \code{\link{rateEvolution}},
#'@keywords ts
#'@examples
#'
#'## load Purkinje cell data recorded in cell-attached mode
#'data(sPK)
#'## coerce sPK to a spikeTrain object
#'sPK <- lapply(sPK, as.spikeTrain)
#'## get a rate evolution data frame
#'sPKreDF <- mkREdf(sPK)
#'## display result using coplot
#'coplot(rate ~ time | longitudinal,data=sPKreDF,panel=lines,show.given=FALSE)
#'\dontrun{
#'## make it prettier with with xyplot of package lattice
#'library(lattice)
#'xyplot(rate ~ time | longitudinal, data=sPKreDF,panel=panel.lines)
#'## if ggplot2 is installed, try it out
#'library(ggplot2)
#'qplot(time,rate,data=sPKreDF,geom="line",colour=longitudinal)
#'}
#'
#'## load Purkinje cell data recorded with the NeuroNexus probes
#'data(mPK)
#'mPK <- lapply(mPK, as.repeatedTrain)
#'## get a rate evolution data frame
#'mPKreDF <- mkREdf(mPK)
#'## use coplot to display result
#'coplot(rate ~ time | longitudinal * across,data = mPKreDF,panel=lines)
#'\dontrun{
#'## make it prettier with with xyplot of package lattice
#'library(lattice)
#'xyplot(rate ~ time | across,data = mPKreDF,groups=longitudinal,panel=panel.lines)
#'xyplot(rate ~ time | across * longitudinal,data = mPKreDF, panel=panel.lines)
#'## if ggplot2 is installed, try it out
#'library(ggplot2)
#'qplot(time,rate,data=mPKreDF,geom="line",colour=longitudinal,facets=across ~ .)
#'}
#'
#'## another example with the CAL1V data set
#'data(CAL1V)
#'CAL1V <- lapply(CAL1V,as.repeatedTrain)
#'## generate the data frame specifying the longitudinal argument
#'## to end up with a clearer display
#'CAL1VreDF <- mkREdf(CAL1V,longitudinal=paste(1:20))
#'coplot(rate ~ time | across * longitudinal,data=CAL1VreDF,panel=lines,show.given=FALSE)
#'\dontrun{
#'## if ggplot2 is installed, try it out
#'library(ggplot2)
#'qplot(time,rate,data=CAL1VreDF,geom="line",facets=longitudinal ~ across)
#'}
#'
#'## another example with the CAL2C data set
#'data(CAL2C)
#'CAL2C <- lapply(CAL2C,as.repeatedTrain)
#'## generate the data frame specifying the longitudinal argument
#'## to end up with a clearer display
#'CAL2CreDF <- mkREdf(CAL2C,longitudinal=paste(1:20))
#'coplot(rate ~ time | across * longitudinal,data=CAL2CreDF,panel=lines,show.given=FALSE)
#'\dontrun{
#'## if ggplot2 is installed, try it out
#'library(ggplot2)
#'qplot(time,rate,data=CAL2CreDF,geom="line",facets=longitudinal ~ across)
#'}
#'
#'
mkREdf <- function(x,
                   longitudinal,
                   across,
                   bw,
                   kernel=c("gaussian", 
                     "epanechnikov", "rectangular", "triangular", "biweight", 
                     "cosine", "optcosine"), 
                   n=512,
                   from, to,
                   na.rm=FALSE,
                   minusMean=FALSE
                   ) {

  ## check argument x
  if (!inherits(x,"list")) x <- list(x=as.spikeTrain(x))
  if (!is.spikeTrain(x[[1]]) && !is.repeatedTrain(x[[1]]))
    stop("Wrong x.")

  xST <- is.spikeTrain(x[[1]])
  if (missing(longitudinal)) {
    if (xST) longitudinal <- names(x)
    else longitudinal <- names(x[[1]])
  }
  if (missing(across)) {
    if (xST) across <- deparse(substitute(x))
    else across <- names(x)
  } 
  allDF <- expand.grid(longitudinal,across)
  nbT <- dim(allDF)[1]

  theList <- list(kernel=kernel,n=n,na.rm=na.rm)
  if (!missing(bw)) {
    bwGiven <- TRUE
    bw <- rep(bw,length.out=nbT)
  } else {
    bwGiven <- FALSE
  }
  if (!missing(from)) {
    theList$from <- from
  } else {
    if (xST) theList$from <- floor(min(sapply(x,min)))
    else theList$from <- floor(min(sapply(x,function(l) min(sapply(l,min)))))
  }
  if (!missing(to)) {
    theList$to <- to
  } else {
    if (xST) theList$to <- ceiling(max(sapply(x,max)))
    else theList$to <- ceiling(max(sapply(x,function(l) max(sapply(l,max)))))
  }
  
  myRE <- function(idx) {
    if (xST) list4call <- c(list(x=x[[idx]]),theList)
    else list4call <- c(list(x=x[[allDF[idx,2]]][[allDF[idx,1]]]),
                        theList)
    if (bwGiven) list4call$bw <- bw[idx]
    re <- do.call(rateEvolution,list4call)[c("x","y")]
    re$longitudinal <- rep(allDF[idx,1],length(re$x))
    re$across <- rep(allDF[idx,2],length(re$x))
    re
  }
                        
  result <- lapply(1:nbT, myRE)
  result <- data.frame(time=unlist(lapply(result,function(l) l$x)),
                       rate=unlist(lapply(result,function(l) l$y)),
                       longitudinal=unlist(lapply(result,function(l) l$longitudinal)),
                       across=unlist(lapply(result,function(l) l$across))
                       )

  if (minusMean) {
    meanR <- sapply(across,
                    function(n)
                    apply(sapply(longitudinal,
                                 function(l)
                                 result$rate[result$across==n & result$longitudinal==l]
                                 ),1,mean)
                    )
    colnames(meanR) <- across
    result$rate <- as.vector(sapply(across,
                                    function(n)
                                    sapply(longitudinal,
                                           function(l)
                                           result$rate[result$across==n & result$longitudinal==l]-
                                           meanR[,n]
                                           )
                                    )
                             )
  }
  result
  
}
