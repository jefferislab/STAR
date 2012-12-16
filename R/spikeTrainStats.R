#'Variance-Time Analysis for Spike Trains
#'
#'Performs Variance-Time Analysis for a Spike Train (or any univariate time
#'series) assuming a Poisson Process with the same Rate as the Spike Train.
#'
#'See Fig. 5 of Ogata (1988) for details. The confidence intervals are obtained
#'with a Normal approximation of the Poisson distribution.
#'
#'@aliases varianceTime plot.varianceTime is.varianceTime
#'@param spikeTrain a \code{spikeTrain} object or a vector which can be coerced
#'to such an object.
#'@param CI a numeric vector with at most two elements. The coverage
#'probability of the confidence intervals.
#'@param windowSizes a numeric increasing vector of positive numbers. The
#'window sizes used to split the spike train.
#'@return \code{varianceTime} returns a list of class \code{varianceTime} with
#'the following elements:
#'
#'\code{plot.varianceTime} is used for its side effect: a graph is produced.
#'
#'\code{is.varianceTime} returns \code{TRUE} if its argument is a
#'\code{varianceTime} object and \code{FALSE} otherwise.
#'@returnItem s2 numeric vector of empirical variance.
#'@returnItem sigma2 numeric vector of expected variance under the Poisson
#'hypothesis.
#'@returnItem ciUp a numeric vector or a 2 rows matrix with the upper limits of
#'the confidence interval(s).
#'@returnItem ciLow a numeric vector or a 2 rows matrix with the lower limits
#'of the confidence interval(s).
#'@returnItem windowSizes numeric vector of window sizes actually used.
#'@returnItem CI a numeric vector, the coverage probabilities of the confidence
#'intervals.
#'@returnItem call the matched call
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com} and Chong Gu
#'\email{chong@@stat.purdue.edu} for a correction on the sampling variance of
#'the variance of a normal distribution.
#'@seealso \code{\link{acf.spikeTrain}}, \code{\link{renewalTestPlot}}
#'@references Ogata, Yosihiko (1988) Statistical Models for Earthquake
#'Occurrences and Residual Analysis for Point Processes. \emph{Journal of the
#'American Statistical Association} \bold{83}: 9-27.
#'@keywords ts survival
#'@examples
#'
#'## Replicate (almost) Fig. 5 of Ogata 1988
#'data(ShallowShocks)
#'vtShallow <- varianceTime(ShallowShocks$Date,,c(5,10,20,40,60,80,seq(100,500,by = 25))*10)
#'is.varianceTime(vtShallow)
#'plot(vtShallow, style="Ogata")
#'
varianceTime <- function (spikeTrain,
                          CI = c(0.95, 0.99),
                          windowSizes) {

  if (!is.spikeTrain(spikeTrain)) spikeTrain <- as.spikeTrain(spikeTrain)
  
  train.l <- length(spikeTrain)
  lastTime <- spikeTrain[train.l]
  mean.rate <- train.l/(lastTime - spikeTrain[1])
  if (!missing(windowSizes)) {
    windowSizes <- sort(windowSizes)
    if (any(windowSizes <= 0)) 
      stop(paste("Elements of", deparse(substitute(windowSizes)), 
                 "should be non-negative"))
    if (max(windowSizes) >= lastTime/5) 
      warning(paste("Some elements of", deparse(substitute(windowSizes)), 
                    "are too large and will be ignored"))
    windowSizes <- windowSizes[windowSizes < lastTime/5]
  } else {
    minSize <- 5/mean.rate
    maxSize <- lastTime/10
    windowSizes <- seq(minSize, maxSize, minSize)
    windowSizes <- windowSizes[windowSizes <= maxSize]
  }
  if (length(CI) > 2) 
    CI <- CI[1:2]
  if (any(CI >= 1 | CI <= 0)) 
    stop(paste(deparse(substitute(CI)), "components should be in (0,1)"))
  dataVariance <- sapply(windowSizes,
                         function(ws) {
                           nBins <- lastTime%/%ws
                           binSeq <- (0:nBins) * ws
                           counts <- hist(spikeTrain[spikeTrain > min(binSeq) & spikeTrain < max(binSeq)],
                                          breaks = binSeq, plot = FALSE)$counts
                           s2 <- var(counts)
                           sigma2 <- mean.rate * ws
                           ciUp <- sapply(CI, function(ci)
                                          qnorm(1 - (1 - ci)/2,sigma2, sqrt((2*sigma2^2+sigma2)/nBins)))
                           ciLow <- sapply(CI, function(ci)
                                           max(qnorm((1 - ci)/2,sigma2, sqrt((2*sigma2^2+sigma2)/nBins)), 0))
                           c(s2, sigma2, ciUp, ciLow)
                         }
                         )
  if (length(CI) == 1) { 
    resutl <- list(s2 = dataVariance[1,],
                   sigma2 = dataVariance[2,],
                   ciUp = dataVariance[3,],
                   ciLow = dataVariance[4,],
                   windowSizes = windowSizes,
                   CI = CI,
                   call = match.call()
                   )
  } else {
    result <- list(s2 = dataVariance[1, ],
                   sigma2 = dataVariance[2,],
                   ciUp = dataVariance[c(3, 4), ],
                   ciLow = dataVariance[c(5,6), ],
                   windowSizes = windowSizes,
                   CI = CI,
                   call = match.call()
                   )
  }
  class(result) <- "varianceTime"
  return(result)

}

#'@param obj a object to test against a \code{varianceTime} object.
#'@return logical
#'@export
#'@rdname varianceTime
is.varianceTime <- function(obj) {

  if(!("varianceTime" %in% class(obj))) return(FALSE)
  
  componentsNames <- c("s2",
                       "sigma2",
                       "ciUp",
                       "ciLow",
                       "call", 
                       "windowSizes",
                       "CI")
  
  if (!all(names(obj) %in% componentsNames)) return(FALSE)
  TRUE
}


#'@param x a \code{varianceTime} object.
#'@param style a character. The style of the plot, \code{"default"} or
#'\code{"Ogata"}.
#'@param unit a character. The unit in which the spike times are expressed.
#'@param xlab a character. The x label.
#'@param ylab a character. The y label.
#'@param main a character. The title.
#'@param sub a character. The subtitle.
#'@param xlim a numeric. See \code{\link{plot}}.
#'@param ylim a numeric. See \code{\link{plot}}.
#'@param \dots see \code{\link{plot}}.
#'@rdname varianceTime
plot.varianceTime <- function (x,
                               style = c("default", "Ogata"), 
                               unit = "s", xlab, ylab, main,
                               sub, xlim, ylim,
                               ...) {

  varianceTimeObj <- x
  if (!is.varianceTime(varianceTimeObj))
    stop(paste(deparse(substitute(varianceTimeObj)), "is not a varianceTime object."))
  
  if (missing(xlab)) 
    xlab <- paste("Time (", unit, ")", sep = "")
  if (missing(ylab)) 
    ylab <- "Variance"
  if (missing(main)) 
    main <- "Estimated Variance-Time Curve and Theoretical Poisson"
  if (missing(xlim)) 
    xlim <- c(0, max(varianceTimeObj$windowSizes))
  if (missing(ylim)) {
    ylim <- c(0, max(c(varianceTimeObj$s2, varianceTimeObj$ciUp)) * 1.01)
  }
  if (missing(sub)) {
    if (is.null(varianceTimeObj$call[["CI"]])) {
      CItxt <- paste(eval(formals(varianceTime)$CI) * 100,collapse = " and ")
    } else {
      CItxt <- paste(eval(varianceTimeObj$call[["CI"]]) * 100, collapse = " and ")
    }
    sub <- paste("Conf. Interval at: ", CItxt, sep = "")
  }
  X <- varianceTimeObj$windowSizes
  if (style[1] == "Ogata") {
    plot(X, varianceTimeObj$s2, type = "n", xlab = xlab, 
         ylab = ylab, main = main, xlim = xlim, ylim = ylim, 
         sub = sub, xaxs = "i", yaxs = "i")
    if (is.null(dim(varianceTimeObj$ciLow))) {
      lines(X, varianceTimeObj$ciLow, lty = 2)
      lines(X, varianceTimeObj$ciUp, lty = 2)
    } else {
      apply(varianceTimeObj$ciLow, 1, function(Y) lines(X,Y, lty = 2))
      apply(varianceTimeObj$ciUp, 1, function(Y) lines(X,Y, lty = 2))
    }
    lines(X, varianceTimeObj$sigma2)
    points(X, varianceTimeObj$s2, pch = 3, ...)
  } else {
    plot(X, varianceTimeObj$s2, type = "n", xlab = xlab, 
         ylab = ylab, main = main, xlim = xlim, ylim = ylim, 
         xaxs = "i", yaxs = "i", sub = sub)
    if (is.null(dim(varianceTimeObj$ciUp))) {
      polygon(c(X, rev(X)),
              c(varianceTimeObj$ciUp,rev(varianceTimeObj$ciLow)),
              lty = 0,
              col = "grey80")
    } else {
      polygon(c(X, rev(X)),
              c(varianceTimeObj$ciUp[2,],rev(varianceTimeObj$ciLow[2, ])),
              lty = 0, 
              col = "grey30")
      polygon(c(X, rev(X)),
              c(varianceTimeObj$ciUp[1,],rev(varianceTimeObj$ciLow[1, ])),
              lty = 0, 
              col = "grey80")
    }
    lines(X, varianceTimeObj$sigma2, lwd = 2)
    lines(X, varianceTimeObj$s2, col = 2, lwd = 2, ...)
  }
  
}




#'Auto- Covariance and -Correlation Function Estimation for Spike Train ISIs
#'
#'The function \code{acf.spikeTrain} computes (and by default plots) estimates
#'of the autocovariance or autocorrelation function of the inter-spike
#'intervals of a spike train.
#'
#'Just a wrapper for \code{\link{acf}} function. The first argument,
#'\code{spikeTrain}, is processed first to extract the inter-spike intervals.
#'\code{acf.spikeTrain} is mainly used to plot what Perkel et al (1967) call
#'the \emph{serial correlation coefficient} (Eq. 8) or \emph{serial covariance
#'coefficient} (Eq. 7), p 400.
#'
#'@param spikeTrain a \code{spikeTrain} object or a vector which can be coerced
#'to such an object.
#'@param lag.max maximum lag at which to calculate the acf.  Default is
#'\eqn{10\log_{10}(N)}{10*log10(N)} where \eqn{N} is the number of ISIs.  Will
#'be automatically limited to one less than the number of ISIs in the spike
#'train.
#'@param type character string giving the type of acf to be computed.  Allowed
#'values are \code{"correlation"} (the default), \code{"covariance"} or
#'\code{"partial"}.
#'@param plot logical. If \code{TRUE} (the default) the acf is plotted.
#'@param na.action function to be called to handle missing values.
#'\code{na.pass} can be used.
#'@param demean logical.  Should the covariances be about the sample means?
#'@param xlab x axis label.
#'@param ylab y axis label.
#'@param main title for the plot.
#'@param \dots further arguments to be passed to \code{plot.acf}.
#'@return An object of class \code{"acf"}, which is a list with the following
#'elements:
#'
#'The lag \code{k} value returned by \code{ccf(x,y)} estimates the correlation
#'between \code{x[t+k]} and \code{y[t]}.
#'
#'The result is returned invisibly if \code{plot} is \code{TRUE}.
#'@returnItem lag A three dimensional array containing the lags at which the
#'acf is estimated.
#'@returnItem acf An array with the same dimensions as \code{lag} containing
#'the estimated acf.
#'@returnItem type The type of correlation (same as the \code{type} argument).
#'@returnItem n.used The number of observations in the time series.
#'@returnItem series The name of the series \code{x}.
#'@returnItem snames The series names for a multivariate time series.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{acf}}, \code{\link{varianceTime}},
#'\code{\link{renewalTestPlot}}
#'@references Perkel D. H., Gerstein, G. L. and Moore G. P. (1967) Neural Spike
#'Trains and Stochastic Point Processes. I. The Single Spike Train.
#'\emph{Biophys. J.}, \bold{7}: 391-418.
#'\url{http://www.pubmedcentral.nih.gov/articlerender.fcgi?tool=pubmed&pubmedid=4292791}
#'@keywords ts
#'@examples
#'
#'## Simulate a log normal train
#'train1 <- c(cumsum(rlnorm(301,log(0.01),0.25)))
#'train1 <- as.spikeTrain(train1)
#'## Get its isi acf
#'acf.spikeTrain(train1,lag.max=100)
#'
acf.spikeTrain <- function (spikeTrain,
                            lag.max = NULL,
                            type = c("correlation", "covariance", "partial"),
                            plot = TRUE,
                            na.action = na.fail, 
                            demean = TRUE,
                            xlab = "Lag (in isi #)",
                            ylab = "ISI acf", 
                            main,
                            ...) {

  if (!is.spikeTrain(spikeTrain)) spikeTrain <- as.spikeTrain(spikeTrain)
  isi <- diff(spikeTrain)
  if (missing(main)) 
    main <- paste("Train", deparse(substitute(spikeTrain)), "ISI acf")
  acf(isi, lag.max = lag.max, type = type, plot = plot, na.action = na.action, 
      demean = demean, xlab = xlab, main = main, ylab = ylab, 
      ...)
  
}



#'Non-Parametric Tests for Renewal Processes
#'
#'Performs and displays rank based tests checking if a spike train is a renewal
#'process
#'
#'\code{renewalTestPlot} generates a 4 panel plot. The 2 graphs making the top
#'row are qualitative and display the rank of inter-spike interval (ISI) k+1
#'versus the rank of ISI k (left graph) and the rank of ISI k+2 versus the one
#'of ISI k (right graph).  The bottom left graph displays the autocorrelation
#'function of the ISIs and is generated by a call to
#'\code{\link{acf.spikeTrain}}. The bottom right graph display the result of a
#'Chi square test performed on the ranks at different lags. More precisely, for
#'each considered lag j (from 1 to \code{lag.max}) the square within which the
#'rank of ISI k+1 vs the one of ISI k is found is splited in \eqn{d^2}{d^2}
#'cells. This decomposition into cells is shown on the two graphs of the top
#'row. Under the renewal process hypothesis the points should be uniformly
#'distributed with a density \eqn{\frac{N}{d^2}}{N/d^2}, where N is the number
#'of ISIs. The sum other rows and other columns is moreover exactly
#'\eqn{\frac{N}{d}}{N/d}. The upper graphs are therefore graphical displays of
#'two-dimensional contingency tables. A chi square test for two-dimensional
#'contingency tables (function \code{\link{chisq.test}}) is performed on the
#'table generated at each lag j. The resulting Chi 2 value is displayed vs the
#'lag. The 95\% confidence region appears as a clear grey rectangle, the value
#'falling within this region appear as black dots and the ones falling out
#'appear as dark grey triangles.
#'
#'@param spikeTrain a \code{spikeTrain} object or a vector which can be coerced
#'to such an object.
#'@param lag.max argument passed to \code{\link{acf.spikeTrain}}.
#'@param d an integer >= 2, the number of divisions used for the Chi 2 test.
#'The default value is such that under the null hypothesis at least 25 events
#'should fall in each division.
#'@param orderPlotPch \code{pch} argument for the order plots.
#'@param \dots additional arguments passed to function
#'\code{\link{chisq.test}}.
#'@return Nothing is returned, the function is used for its side effect: a plot
#'is generated.
#'@note You should not use a too large value for \code{d} otherwise the Chi 2
#'values will be too approximative and warnings will be printed.  If your
#'process is a renewal process you should have on average 5\% of the points on
#'the bottom right graph appearing as dark triangles.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{acf}}, \code{\link{varianceTime}},
#'\code{\link{acf.spikeTrain}}
#'@keywords ts survival
#'@examples
#'
#'\dontrun{
#'## Apply the test of Ogata (1988) shallow shock data
#'data(ShallowShocks)
#'renewalTestPlot(ShallowShocks$Date,d=3)
#'
#'## Apply the test to the second and third neurons of the cockroachAlSpont
#'## data set
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
#'## next the renewal tests
#'renewalTestPlot(CAL1S[["neuron 1"]])
#'
#'## Simulate a renewal log normal train with 500 isi
#'isi.nb <- 500
#'train1 <- c(cumsum(rlnorm(isi.nb+1,log(0.01),0.25)))
#'## make the test
#'renewalTestPlot(train1)
#'
#'## Simulate a (non renewal) 2 states train
#'myTransition <- matrix(c(0.9,0.1,0.1,0.9),2,2,byrow=TRUE)
#'states2 <- numeric(isi.nb+1) + 1
#'for (i in 1:isi.nb) states2[i+1] <- rbinom(1,1,prob=1-myTransition[states2[i],])+1
#'myLnormPara2 <- matrix(c(log(0.01),0.25,log(0.05),0.25),2,2,byrow=TRUE)
#'train2 <-
#'cumsum(rlnorm(isi.nb+1,myLnormPara2[states2,1],myLnormPara2[states2,2]))
#'## make the test
#'renewalTestPlot(train2)}
#'
renewalTestPlot <- function (spikeTrain,
                             lag.max=NULL,
                             d=max(c(2,sqrt(length(spikeTrain)) %/% 5)),
                             orderPlotPch=ifelse(length(spikeTrain)<=600,1,"."), 
                             ...) {

  if (!is.spikeTrain(spikeTrain)) spikeTrain <- as.spikeTrain(spikeTrain)
  if (length(spikeTrain) < 50) 
    stop(paste(deparse(substitute(spikeTrain)), "contains less than 50 events."))
  m <- matrix(c(1:4), nrow = 2, ncol = 2, byrow = TRUE)
  oldpar <- par(mar = c(5, 4, 2, 2))
  layout(m)
  on.exit(par(oldpar))
  isi <- diff(spikeTrain)
  isi.o <- rank(isi)/length(isi)
  isi.l <- length(isi)
  if (is.null(lag.max)) 
    lag.max <- round(10 * log10(isi.l))
  lag.max <- min(isi.l - 1, lag.max)
  grid <- seq(0, 1, length.out = d + 1)
  getChi2 <- function(lag) {
    isi.x <- isi.o[1:(isi.l - lag.max)]
    isi.y <- isi.o[(1 + lag):(isi.l - lag.max + lag)]
    isi.x <- as.integer(cut(isi.x, breaks = grid))
    isi.y <- as.integer(cut(isi.y, breaks = grid))
    counts <- matrix(0, nrow = d, ncol = d)
    for (i in seq(along.with = isi.x))
      counts[isi.x[i], isi.y[i]] <- counts[isi.x[i], isi.y[i]] + 1
    chisq.test(counts, ...)
  }
  chi2seq <- lapply(1:lag.max, getChi2)
  minChi2 <- qchisq(0.025, df = chi2seq[[1]]$parameter)
  maxChi2 <- qchisq(0.975, df = chi2seq[[1]]$parameter)
  chi2V <- sapply(chi2seq, function(l) l$statistic)
  outOf95 <- chi2V < minChi2 | chi2V > maxChi2

  plot(isi.o[-length(isi.o)] * isi.l, isi.o[-1] * isi.l, xlab = quote(O[k]), 
       ylab = quote(O[k + 1]), main = "Order Statistic Correlation at Lag 1", 
       type = "n")

  sapply(grid[-c(1, d + 1)],
         function(idx) {
           abline(v = idx * isi.l, col = "grey")
           abline(h = idx * isi.l, col = "grey")
         }
         )

  points(isi.o[-length(isi.o)] * isi.l, isi.o[-1] * isi.l, 
         pch = orderPlotPch)
  plot(isi.o[-(0:-1 + length(isi.o))] * isi.l,
       isi.o[-(1:2)] * isi.l,
       xlab = quote(O[k]), ylab = quote(O[k + 2]),
       main = "Order Statistic Correlation at Lag 2", 
       type = "n")
  
  sapply(grid[-c(1, d + 1)],
         function(idx) {
           abline(v = idx * isi.l, col = "grey")
           abline(h = idx * isi.l, col = "grey")
         }
         )
  
  points(isi.o[-(0:-1 + length(isi.o))] * isi.l,
         isi.o[-(1:2)] * isi.l,
         pch = orderPlotPch)
  
  acf.spikeTrain(spikeTrain,
                 lag.max = lag.max, main = "")

  plot(1:lag.max, chi2V, type = "n",
       xlab = "Lag (in isi #)", 
       ylab = quote(chi^2),
       main = expression(paste(chi^2, " Statistics")), 
       xlim = c(0, lag.max),
       ylim = c(min(qchisq(0.01, df = chi2seq[[1]]$parameter), min(chi2V)),
         max(qchisq(0.99, df = d^2 - 2 * (d - 1) - 1), max(chi2V)))
       )
  
  polygon(c(0, 0, lag.max + 1, lag.max + 1),
          c(minChi2, maxChi2, maxChi2, minChi2),
          col = "grey80", lty = 0
          )

  points(1:lag.max, chi2V,
         pch = c(16, 17)[outOf95 + 1],
         col = c("black", "grey30")[outOf95 + 1]
         )

}
