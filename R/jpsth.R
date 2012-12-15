jsd <- function(xRT,
                yRT,
                acquisitionWindow,
                xlab,
                ylab,
                main,
                pch=".",
                ...
                )
### computes a joint peristimulus scatter diagram
### xRT and yRT should be repeatedTrain objects.
### The spike times of xRT are used to get the x coordinate
### The spike times of yRT are used to get the y coordinate
### acquisitionWindow the acquisition window used or a subset
### of it.
{
  xRTn <- deparse(substitute(xRT))
  yRTn <- deparse(substitute(yRT))
  
  ## check that xRT and yRT are repeatedTrain objects of the
  ## same length
  if (!is.repeatedTrain(xRT))
    stop(paste(xRTn,"should be a \"repeatedTrain\" object."))
  if (!is.repeatedTrain(yRT))
    stop(paste(yRTn,"should be a \"repeatedTrain\" object."))
  nbTrials <- length(xRT)
  if (nbTrials != length(yRT))
    stop(paste(xRTn,"and", yRTn,
               "should be \"repeatedTrain\" objects of the same length.")
         )

  if (missing(acquisitionWindow)) {
    acquisitionWindow <- range(c(unlist(xRT),unlist(yRT)))
    acquisitionWindow <- c(floor(acquisitionWindow[1]),
                           ceiling(acquisitionWindow[2])
                           )
  }

  if (missing(xlab))
    xlab <- paste(xRTn, "spike times (s)")

  if (missing(ylab))
    ylab <- paste(yRTn, "spike times (s)")

  if (missing(main))
    main <- paste("Joint peri-stimulus scatter diagram of:",
                  xRTn,"and",yRTn)

  plot(acquisitionWindow,
       acquisitionWindow,
       type="n",
       xlab=xlab,
       ylab=ylab,
       main=main,
       ...)

  invisible(sapply(1:nbTrials,
                   function(tIdx) {
                     y <- as.numeric(yRT[[tIdx]])
                     X <- as.numeric(xRT[[tIdx]])
                     sapply(X,
                            function(x) {
                              x <- rep(x,length(y))
                              points(x,y,pch=pch)
                            }
                            )
                   }
                   )
            )
  
}



#'Related Functions and Methods for Joint-PSTHs and Joint Scatter Diagrams
#'
#'Some mainly graphical tools to probe interactions between 2 neurons recorded
#'in the presence of a repeated stimulation.
#'
#'The joint scatter diagram was introduced by Gerstein and Perkel (1972). The
#'joint peristimulus time histogram is a binned version of it (Aertsen et al,
#'1989). \code{jpsth2df} allows the reformating of a \code{jpsth} object in
#'order to compute a smooth version of it with \code{\link[gss]{gssanova}},
#'\code{\link[gss]{gssanova0}} or \code{\link[mgcv]{gam}}.
#'
#'@aliases jpsth jsd contour.jpsth image.jpsth persp.jpsth jpsth2df
#'@param xRT a \code{repeatedTrain} object whose spike times will appear on the
#'abscissa of the plots.
#'@param yRT a \code{repeatedTrain} object whose spike times will appear on the
#'ordinate of the plots. It must have the same length as \code{xRT}.
#'@param x,object \code{jpsth} objects.
#'@param xBreaks,yBreaks A single number (the bin width) or a vector defining
#'bins boundaries on the X and Y axis. If missing a default is provided.
#'@param acquisitionWindow a 2 elements vector specifying the begining and the
#'end of the acquisition. If missing values are obtained using the
#'\code{\link{floor}} of the smallest spike time and the \code{\link{ceiling}}
#'of the largest one.
#'@param nbEvtPerBin If both \code{xBreaks} and \code{xBreaks} are missing a
#'bin width, \code{bw}, is computed such that the expected value of the count
#'per cell (2 dimensional bin) would be \code{nbEvtPerBin} assuming a
#'stationary Poisson discharge for both neurons.
#'@param xlab a character (default value supplied). See \code{\link{plot}}.
#'@param ylab a character (default value supplied). See \code{\link{plot}}.
#'@param main a character (default value supplied). See \code{\link{plot}}.
#'@param pch the type of "points" displayed by \code{jsd}. See
#'\code{\link{plot}}.
#'@param \dots additional arguments passed to \code{\link{plot}} by \code{jsd}
#'and to respective generic methods by \code{contour.jpsth}, \code{image.jpsth}
#'and \code{persp.jpsth}.
#'@return \code{jsd} is used for its side effect, a plot is generated and
#'nothing is returned.
#'
#'\code{jpsth2df} returns a \code{\link{data.frame}} with the following
#'variables: \code{Count}, the counts per cell; \code{X}, the position of the
#'cell on the X axis; \code{Y}, the position of the cell on the Y axis; and
#'attributes: \code{xBreaks}, \code{yBreaks}, \code{xTotal}, \code{yTotal},
#'\code{nbTrials}, \code{acquisitionWindow} corresponding to the components of
#'its argument with the same name and \code{originalCall} corresponding to
#'component \code{call}.
#'
#'\code{jpsth} returns a list of class \code{jpsth} with the following
#'components:
#'@returnItem counts a matrix storing the counts per cell.
#'@returnItem density a matrix storing the density in each cell.
#'@returnItem xMids a vector containing the X positions of the cells.
#'@returnItem yMids a vector containing the Y positions of the cells.
#'@returnItem xBreaks a vector containing the bin boundaries of the cells along
#'the X axis.
#'@returnItem yBreaks a vector containing the bin boundaries of the cells along
#'the X axis.
#'@returnItem xTotal the total number of spikes of the "X" neuron.
#'@returnItem yTotal the total number of spikes of the "Y" neuron.
#'@returnItem xFreq the mean freqency of the "X" neuron.
#'@returnItem yFreq the mean freqency of the "Y" neuron.
#'@returnItem nbTrials the number of trials of \code{xRT} (and \code{yRT}).
#'@returnItem acquisitionWindow the boundaries of the acquisition window.
#'@returnItem call the matched call.
#'@note I use "joint scatter diagram" for what Gerstein and Perkel (1972) more
#'properly call a "joint peristimulus time scatter diagram".
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{lockedTrain}}, \code{\link{plot.lockedTrain}},
#'\code{\link{hist.lockedTrain}}, \code{\link{gsslockedTrain}},
#'\code{\link{plot.gsslockedTrain}}, \code{\link{gsslockedTrain0}},
#'\code{\link{plot.gsslockedTrain0}}, \code{\link{gamlockedTrain}},
#'\code{\link{plot.gamlockedTrain}}, \code{\link{contour}},
#'\code{\link{image}}, \code{\link{persp}}, \code{\link{attr}},
#'\code{\link{attributes}}
#'@references Gerstein, G. L. and Perkel, D. H. (1972) Mutual temporal
#'relationships among neuronal spike trains. Statistical techniques for display
#'and analysis. \emph{Biophys J} \bold{12}: 453--473.
#'\url{http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=1484144}
#'
#'Aertsen, A. M., Gerstein, G. L., Habib, M. K., Palm, G. (1989) Dynamics of
#'neuronal firing correlation: modulation of "effective connectivity". \emph{J
#'Neurophysiol} \bold{61}: 900--917.
#'\url{http://jn.physiology.org/cgi/content/abstract/61/5/900}
#'@keywords models
#'@examples
#'
#'## load e070528citronellal data
#'data(e070528citronellal)
#'## plot a jsd with neuron 1 on X and neuron 2 on Y
#'jsd(e070528citronellal[[1]],e070528citronellal[[2]])
#'## now make the jpsth
#'j1.2 <- jpsth(e070528citronellal[[1]],e070528citronellal[[2]])
#'## make a contour plot
#'contour(j1.2)
#'## make an image plot
#'image(j1.2)
#'## make a persp plot
#'persp(j1.2)
#'\dontrun{
#'## fit a gss model with interactions
#'## use a larger bin width for the jpsth
#'j1.2 <- jpsth(e070528citronellal[[1]],e070528citronellal[[2]],0.2,0.2)
#'## get a data frame
#'j1.2DF <- jpsth2df(j1.2)
#'## To save computation time start analyzing
#'## just before the stimulation time
#'j1.2DF <- j1.2DF[j1.2DF$X > 6 & j1.2DF$Y>6,]
#'gf <- gssanova(Count ~ X*Y, family="poisson", data=j1.2DF,seed=20061001)
#'## Use the project function of gss to check the significance
#'## of the interaction term
#'project(gf2,inc=c("X","Y"))
#'}
#'\dontrun{
#'## fit a gam model assuming no interaction
#'## get a data frame
#'j1.2DF <- jpsth2df(j1.2)
#'fitNoI <- gam(Count ~ s(X,k=100,bs="cr") + s(Y,k=100,bs="cr"),data=j1.2DF,family=poisson())
#'}
#'
jpsth <- function(xRT,
                  yRT,
                  xBreaks,
                  yBreaks,
                  acquisitionWindow,
                  nbEvtPerBin=50
                  )
### computes a joint peristimulus time-histogram
### xRT and yRT should be repeatedTrain objects.
### The spike times of xRT are used to get the x coordinate
### The spike times of yRT are used to get the y coordinate
### xBreaks can be a single number specifying the number
### of equally larged bins to use to split acquisitionWindow
### or a vector of boundaries between bins for the x axis.
### yBreaks can be a single number specifying the number
### of equally larged bins to use to split acquisitionWindow
### or a vector of boundaries between bins for the y axis.
### acquisitionWindow the acquisition window used or a subset
### of it.
{
  ## check that xRT and yRT are repeatedTrain objects of the
  ## same length
  if (!is.repeatedTrain(xRT))
    stop("xRT should be a \"repeatedTrain\" object.")
  if (!is.repeatedTrain(yRT))
    stop("yRT should be a \"repeatedTrain\" object.")
  nbTrials <- length(xRT)
  if (nbTrials != length(yRT))
    stop("xRT and yRT should be \"repeatedTrain\" objects of the same length.")

  if (missing(acquisitionWindow)) 
    acquisitionWindow <- range(c(unlist(xRT),unlist(yRT)))
  
  l <- floor(acquisitionWindow[1])
  r <- ceiling(acquisitionWindow[2])

  if (missing(xBreaks) && missing(yBreaks)) {
    xRT <- lapply(xRT, function(st) as.numeric(st[l < st & st < r]))
    yRT <- lapply(yRT, function(st) as.numeric(st[l < st & st < r]))

    xFreq <- sum(sapply(xRT, length))/(r-l)
    yFreq <- sum(sapply(yRT, length))/(r-l)

    bw <- sqrt(nbEvtPerBin/xFreq/xFreq)
    xBreaks <- seq(l,r,bw)
    yBreaks <- xBreaks
    r <- min(r,xBreaks[length(xBreaks)])
  } ## End of conditional on missing(xBreaks) && missing(yBreaks)

  if (length(xBreaks) == 1) {
    ## xBreaks is interpreted as a bin width
    bw <- xBreaks
    xBreaks <- seq(l,r,bw)
  } ## End of conditional on length(xBreaks) == 1
  
  if (length(yBreaks) == 1) {
    ## yBreaks is interpreted as a bin width
    bw <- yBreaks
    yBreaks <- seq(l,r,bw)
  } ## End of conditional on length(yBreaks) == 1

  if (missing(xBreaks) && !missing(yBreaks)) xBreaks <- yBreaks
  if (missing(yBreaks) && !missing(xBreaks)) yBreaks <- xBreaks

  l <- min(c(xBreaks,yBreaks))
  r <- max(c(xBreaks,yBreaks))
  xRT <- lapply(xRT, function(st) as.numeric(st[l < st & st < r]))
  yRT <- lapply(yRT, function(st) as.numeric(st[l < st & st < r]))

  xTotal <- sum(sapply(xRT, length))
  xFreq <- xTotal/(r-l)
  yTotal <- sum(sapply(yRT, length))
  yFreq <- yTotal/(r-l)

  xRT <- lapply(xRT, function(x) findInterval(x,xBreaks))
  yRT <- lapply(yRT, function(y) findInterval(y,yBreaks))

  counts <- matrix(as.integer(0),
                   nrow=length(xBreaks)-1,
                   ncol=length(yBreaks)-1)

  for (tIdx in 1:nbTrials) {
    x <- xRT[[tIdx]]
    y <- yRT[[tIdx]]
    for (i in seq(x)) counts[x[i],y] <- counts[x[i],y] + 1
  } ## End of for loop on tIdx

  density <- counts/(diff(xBreaks) %o% diff(yBreaks))/xTotal/yTotal

  xMids <- xBreaks[-length(xBreaks)]+diff(xBreaks)/2
  yMids <- yBreaks[-length(yBreaks)]+diff(yBreaks)/2

  result <- list(counts=counts,
                 density=density,
                 xMids=xMids,
                 yMids=yMids,
                 xBreaks=xBreaks,
                 yBreaks=yBreaks,
                 xTotal=xTotal,
                 yTotal=yTotal,
                 xFreq=xFreq,
                 yFreq=yFreq,
                 nbTrials=nbTrials,
                 acquisitionWindow=c(l,r),
                 call=match.call()
                 )

  class(result) <- "jpsth"
  return(result)
                 
}


contour.jpsth <- function(x,
                          xlab,
                          ylab,
                          main,
                          ...)
### contour method or jpsth objects
{
  xN <- deparse(substitute(x))
  ## check x
  if (!inherits(x,"jpsth"))
    stop(paste(xN,"should be a jpsth object."))

  if (missing(xlab))
    xlab <- paste(deparse(x$call[["xRT"]]),"spike times (s)")
  if (missing(ylab))
    ylab <- paste(deparse(x$call[["yRT"]]),"spike times (s)")
  if (missing(main))
    main <- paste("JPSTH of",
                  deparse(x$call[["xRT"]]),
                  "and",
                  deparse(x$call[["yRT"]])
                  )
  
  contour(x=x$xMids,
          y=x$yMids,
          z=x$density,
          xlab=xlab,
          ylab=ylab,
          main=main,
          ...)

}

image.jpsth <- function(x,
                        xlab,
                        ylab,
                        main,
                        ...)
### image method or jpsth objects
{
  xN <- deparse(substitute(x))
  ## check x
  if (!inherits(x,"jpsth"))
    stop(paste(xN,"should be a jpsth object."))

  if (missing(xlab))
    xlab <- paste(deparse(x$call[["xRT"]]),"spike times (s)")
  if (missing(ylab))
    ylab <- paste(deparse(x$call[["yRT"]]),"spike times (s)")
  if (missing(main))
    main <- paste("JPSTH of",
                  deparse(x$call[["xRT"]]),
                  "and",
                  deparse(x$call[["yRT"]])
                  )
  
  image(x=x$xMids,
        y=x$yMids,
        z=x$density,
        xlab=xlab,
        ylab=ylab,
        main=main,
        ...)

}

persp.jpsth <- function(x,
                        xlab,
                        ylab,
                        main,
                        ...)
### persp method or jpsth objects
{
  xN <- deparse(substitute(x))
  ## check x
  if (!inherits(x,"jpsth"))
    stop(paste(xN,"should be a jpsth object."))

  if (missing(xlab))
    xlab <- paste(deparse(x$call[["xRT"]]),"spike times (s)")
  if (missing(ylab))
    ylab <- paste(deparse(x$call[["yRT"]]),"spike times (s)")
  if (missing(main))
    main <- paste("JPSTH of",
                  deparse(x$call[["xRT"]]),
                  "and",
                  deparse(x$call[["yRT"]])
                  )
  
  persp(x=x$xMids,
        y=x$yMids,
        z=x$density,
        xlab=xlab,
        ylab=ylab,
        main=main,
        ...)

}


jpsth2df <- function(object)
### function formating a jpsth object
### into a data frame suitable for use
### in glm or gam for instance
{
  objectN <- deparse(substitute(object))
  ## check object
  if (!inherits(object,"jpsth"))
    stop(paste(objectN,"should be a jpsth object."))

  Count <- object$counts
  CountD <- dim(Count)
  dim(Count) <- NULL
  X <- rep(object$xMids,CountD[2])
  Y <- rep(object$yMids,each=CountD[1])
  result <- data.frame(Count=Count,
                       X=X,
                       Y=Y)
  
  attr(result,"xBreaks") <- object$xBreaks
  attr(result,"yBreaks") <- object$yBreaks
  attr(result,"xTotal") <- object$xTotal
  attr(result,"yTotal") <- object$yTotal
  attr(result,"nbTrials") <- object$nbTrials
  attr(result,"acquisitionWindow") <- object$acquisitionWindow
  attr(result,"originalCall") <- object$call

  return(result)
  
}
