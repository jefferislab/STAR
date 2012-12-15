#'Counting Process Sample Paths
#'
#'Functions to create and explore \code{CountingProcessSamplePath} objects.
#'These objects are complementary to the \code{spikeTrain} objects, the latter
#'being in fact point processes representations.
#'
#'\code{CountingProcessSamplePath} objects are complementary to
#'\code{spikeTrain} objects. They are also used to represente slightly more
#'general properties of these objects and are directed towards model testing.
#'
#'More formaly, if we observe n events at times
#'\eqn{\{t_1,\ldots,t_n\}}{{t1,...,tn}} such that, \eqn{from < t_1 < }{from <
#'t1 < ... < tn <= to}\eqn{ \ldots < t_n \le to}{from < t1 < ... < tn <= to},
#'the \emph{counting process sample path} is the right continuous function
#'defined by: \deqn{N(t) = \sharp \{ t_j \; \mathrm{with} \; from < t_j \leq t
#'\}}{N(t) = # {tj : from < tj <= t}} where \eqn{\sharp}{#} stands for the
#'number of elements of a set.
#'
#'@aliases mkCPSP print.CountingProcessSamplePath
#'plot.CountingProcessSamplePath lines.CountingProcessSamplePath as.CPSP
#'@param st A \code{numeric} vector with \emph{strictly} increasing elements.
#'@param from A \code{numeric}, the time at which the counting process
#'obeservation started.
#'@param to A \code{numeric}, the time at which the counting process
#'obeservation ended.
#'@param x A \code{numeric} or a \code{spikeTrain} object for \code{as.CPSP}, a
#'\code{CountingProcessSamplePath} object for \code{print}, \code{plot} and
#'\code{lines}.
#'@param digits An \code{integer}, the number of digits to be used while
#'printing summaries. See \code{\link{round}}.
#'@param y Not used but required by the \code{plot} method definition.
#'@param col,lwd,xlim,ylim,xlab,ylab,main,xaxs,yaxs See \code{\link{plot}}.
#'@param \dots Not used in \code{print} (but included for compatibility with
#'the method definition) otherwise used like in \code{\link{plot}} and
#'\code{\link{lines}}.
#'@return \code{mkCPSP} returns an object of class
#'\code{CountingProcessSamplePath}. This object is a \code{list} with the
#'following components:
#'
#'Functions \code{plot} and \code{lines} are used for their side effects,
#'function \code{print} returns a short description of the object corresponding
#'to the \code{summary} returned by function \code{\link{summary.spikeTrain}}
#'for \code{spikeTrain} objects. Function \code{as.CPSP} returns a
#'\code{CountingProcessSamplePath}.
#'@returnItem cpspFct a right continuous \code{function} of \code{t} returning
#'the number of events whose occurrence time is strictly larger than
#'\code{from} and smaller of equal than \code{t}. \code{t} can be a vector. If
#'\code{missing} the cumulative number of events at the events occurrence time
#'is returned.
#'@returnItem ppspFct a \code{function} that does not take any argument and
#'that returns the sequence of events times, that is, the "point process sample
#'path".
#'@returnItem spikeTrainFct a \code{function} that does not take any argument
#'and that returns the \code{spikeTrain} object associated with the
#'\code{CountingProcessSamplePath} object.
#'@returnItem from argument \code{from} of \code{mkCPSP}.
#'@returnItem to argument \code{to} of \code{mkCPSP}
#'@returnItem call the matched call.
#'@note This functions are directed towards model testing, don't be surprised
#'if they look redundant with the corresponding functions for \code{spikeTrain}
#'objects. An apparent difference of detail with the latter is that no scale
#'(like seconds) is assumed by default for \code{CountingProcessSamplePath}
#'objects. This is to cope in a natural way with the time transformation /
#'rescaling procedures used to test conditional intensity models.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{summary.CountingProcessSamplePath}},
#'\code{\link{print.CountingProcessSamplePath.summary}},
#'\code{\link{plot.CountingProcessSamplePath.summary}},
#'\code{\link{summary.spikeTrain}}, \code{\link{print.spikeTrain}},
#'\code{\link{plot.spikeTrain}}, \code{\link{as.spikeTrain}}
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
#'## A simple illustration with Ogata's hearthquakes data set
#'data(ShallowShocks)
#'plot(mkCPSP(ShallowShocks$Date),
#'     xlab="Time (days)",
#'     main="Shallow Shocks Counting Process of Ogata 1988")
#'## An illustration with on of STAR's data neuroanl dicharge data set
#'data(e060824spont)
#'## Create the object from a spikeTrain
#'n1spt.cp <- as.CPSP(e060824spont[["neuron 1"]])
#'## print it
#'n1spt.cp
#'## plot it
#'plot(n1spt.cp)
#'
mkCPSP <- function(st,
                   from=floor(min(st)),
                   to=ceiling(max(st))
                   ) {

  ## mkCPSPFct returns a CountingProcessSamplePath object
  ## Arguments:
  ##  st: A spike train object or a vector with strictly increasing elements
  ##      containing the spike times (measured in s)
  ##  from: The time (in s) at which the observations started
  ##  to: The time (in s) at which the observations ended
  ## Return
  ##  a list with the following components:
  ##    cpspFct: a function of a single time argument, t, instance of
  ##    the class "CountingProcessSamplePath". 
  ##    The function is cadlag (right continous and admits a limit on the left) between
  ##    from and to. The function returns NA outside of this interval. For a vector 
  ##    argument, a vector is returned. The function contains the original data in its
  ##    closure. The matched call is contained in attribute "call".
  ##    from: argument from
  ##    to: argument to
  ##    call: the matched call

  ## check that st is suitable
  st <- unclass(st)
  if (any(diff(st) <= 0)) stop("Wrong st.")
  
  cpspFct <- function(t) {
    if (missing(t)) t <- st
    result <- numeric(length(t))
    good <- from <= t & t <= to
    result[!good] <- NA
    result[good] <- sapply(t[good], function(x) sum(st <= x))
    result
  }

  ppspFct <- function() st

  spikeTrainFct <- function() as.spikeTrain(st)
  
  result <- list(cpspFct=cpspFct,
                 ppspFct=ppspFct,
                 spikeTrainFct=spikeTrainFct,
                 from=from,
                 to=to,
                 call=match.call()
                 )
  class(result) <- "CountingProcessSamplePath"
  result
  
}

print.CountingProcessSamplePath <- function(x,
                                            digits=5,
                                            ...) {

  cat(paste("A CountingProcessSamplePath object\n"))
  xx <- x$ppspFct()


#'Get Lagged Inter Spike Intervals (ISIs) From Data Frames Generated by mkGLMdf
#'
#'A utility function to create a vector containing the ith preceding inter
#'spike interval (isi) at a given time.
#'
#'Look at the (short) source file for details.
#'
#'@param dataFrame a \code{\link{data.frame}} typically generated by
#'\code{\link{mkGLMdf}}. Should at least contain an \code{event} and a
#'\code{time} variable.
#'@param lag a strictly positive integer. Set to 1 if the previous isi is
#'required, to 2 is the isi preceding the last one is required, etc...
#'@return A \code{numeric} vector with the value of the \code{lag}th isi
#'preceding the time of the corresponding bin center.
#'@note Before plugging the result into \code{\link[gss]{gssanova}}, do not
#'forget to remove the \code{NA} elements (see the example).
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{mkGLMdf}}, \code{\link[gss]{gssanova}},
#'\code{\link{%tt%}}
#'@keywords ts survival htest
#'@examples
#'
#'\dontrun{
#'## load e060517spont data set
#'data(e060517spont)
#'## make a data frame using a 2 ms bin width
#'e060517spontDF <- mkGLMdf(e060517spont,0.002,0,60)
#'## Keep data relevant to neuron 1
#'e060517spontDFn1 <- e060517spontDF[e060517spontDF$neuron == "1",]
#'## get the isi at lag 1 and 2
#'e060517spontDFn1$isi1 <- isi(e060517spontDFn1,lag=1)
#'e060517spontDFn1$isi2 <- isi(e060517spontDFn1,lag=2)
#'## keep only defined elements
#'e060517spontDFn1 <- e060517spontDFn1[!is.na(e060517spontDFn1$isi2),]
#'## split the data set into an "early" and a "late" part
#'e060517spontDFn1e <- e060517spontDFn1[e060517spontDFn1$time <= 30,]
#'e060517spontDFn1l <- e060517spontDFn1[e060517spontDFn1$time > 30,]
#'## Fit the late part
#'e060517spontDFn1lGF <- gssanova(event ~ lN.1*isi1+isi2, data=e060517spontDFn1l, family="binomial", seed=20061001)
#'## Time transform the early part and perform goodness of fit tests
#'e060517spont.n1e.tt <- e060517spontDFn1lGF %tt% e060517spontDFn1e
#'e060517spont.n1e.tt
#'summary(e060517spont.n1e.tt)
#'plot(summary(e060517spont.n1e.tt))
#'}
#'
  isi <- diff(xx)
  isi.m <- min(isi)
  isi.M <- max(isi)
  isi.mu <- mean(isi)
  isi.sd <- sd(isi)
  n <- length(xx)
  from <- x$from
  to <- x$to
  cat(paste("with",n,"events, from:",from,"to:",to,"\n"))
  cat(paste("Smallest inter event interval:",round(isi.m,digits=digits),
            ", largest inter event interval:",round(isi.M,digits=digits),
            "\n"))
  cat(paste("Mean inter event interval:",round(isi.mu,digits=digits),
      ", sd:", round(isi.sd,digits=digits),"\n"))
  
}

plot.CountingProcessSamplePath <- function(x,y,
                                           col,lwd,
                                           xlim,ylim,
                                           xlab,ylab,
                                           xaxs,yaxs,
                                           main,
                                           ...) {

  ## plot.CountingProcessSamplePath plot method for CountingProcessSamplePath objects
  ## Arguments:
  ##  x: A CountingProcessSamplePath object
  ##  y: Not used but required for plot methods.
  ##  xlim, ylim, xlab, ylab, main, lwd, col: Same as in plot
  ##  ...: additional arguments passed to plot

  xx <- x$ppspFct()
  n <- length(xx)
  if (missing(xlim)) xlim <- c(x$from,
                               x$to
                               )
  xx <- c(xx,xlim[2])
  if (missing(ylim)) ylim <- c(sum(xx<xlim[1]),sum(xx<=xlim[2]))
  if (missing(xlab)) xlab <- ""
  if (missing(ylab)) ylab <- "Cumulative number of evts."
  if (missing(main)) main <- "Counting Process Sample Path"
  if (missing(xaxs)) xaxs <- "i"
  if (missing(yaxs)) yaxs <- "i"
  plot(xlim,ylim,type="n",
       xlab=xlab,
       ylab=ylab,
       main=main,
       xaxs=xaxs,
       yaxs=yaxs,
       ...)
  yy <- 1:n
  if (missing(lwd)) lwd <- 1
  if (missing(col)) col <- 1
  segments(xx[-n],yy,xx[-1],yy,lwd=lwd,col=col)
  
}

lines.CountingProcessSamplePath <- function(x,
                                            ...) {

  ## lines.CountingProcessSamplePath lines method for CountingProcessSamplePath objects
  ## Arguments:
  ##  x: A CountingProcessSamplePath object
  ##  ...: additional arguments passed to segments
  
  xx <- evalq(st,envir=environment(x$cpspFct))
  n <- length(xx)
  xlim <- c(x$from,
            x$to
            )
  xx <- c(xx,xlim[2])
  yy <- 1:n
  segments(xx[-n],yy,xx[-1],yy,...)
  
}


as.CPSP <- function(x) {

  if (is.numeric(x)) return(mkCPSP(x))
  if (is.spikeTrain(x)) return(mkCPSP(unclass(x)))
  
}



#'Create and Explore Counting Process Sample Path Summaries
#'
#'These functions / methods are designed to test a
#'\code{CountingProcessSamplePath} object against a uniform Poisson process
#'with rate 1.
#'
#'If the \code{CountingProcessSamplePath} object \code{x} is a the realization
#'of a homogeneous Poisson process then, \emph{conditioned on the number of
#'events observed}, the location of the events (jumps in \eqn{N(t)}{N(t)}) is
#'uniform on the period of observation. This is a basic property of the
#'homogeneous Poisson process derived in Chap. 2 of Cox and Lewis (1966) and
#'Daley and Vere-Jones (2003). Component \code{UniformGivenN} of a
#'\code{CountingProcessSamplePath.summary} list contains the p.value of the
#'Kolmogorov test of this uniform null hypothesis. The first graph generated by
#'the \code{plot} method displays the Kolgorov test graphically (\emph{i.e.},
#'the empirical cumulative distribution function and the null hyptohesis (the
#'diagonal). The two dotted lines on both sides of the diagonal correspond to
#'95 and 99\% (asymptotic) confidence intervals. This is the graph shown on
#'Fig. 9 (p 19) of Ogata (1988). Notice that the \code{summary} method allows
#'you to compute the exact p.value.
#'
#'If we write \eqn{x_i}{x[i]} the jump locations of the
#'\code{CountingProcessSamplePath} object \code{x} and if the latter is the
#'realization of a homogeneous Poisson process then the intervals: \deqn{y_i =
#'x_{i+1}-x_i}{y[i]=x[i+1]-x[i]} are realizations of iid rvs from an
#'exponential distribution with rate 1 and the: \deqn{u_i = 1-\exp
#'(-y_i}{u[i]=1 - exp(-y[i])} are realizations of iid rvs from a uniform
#'distribution on [0,1). The second graph generated by the \code{plot} method
#'tests this uniform distribution hypotheses with a Kolmogorov Test. This is
#'the graph shown on Fig. 10 (p 19) of Ogata (1988) which was suggested by
#'Berman. This is also the one of the graphs proposed by Brown et al (2002)
#'(the other one is a Q-Q plot for the same quantities). The two dotted lines
#'on both sides of the diagonal correspond to 95 and 99\% (asymptotic)
#'confidence intervals. Component \code{BermanTest} of a
#'\code{CountingProcessSamplePath.summary} list contains the p.value of the
#'Kolmogorov test of this uniform null hypothesis.
#'
#'Following the line of the previous paragraph, if the distribution of the
#'\eqn{y_i}{y[i]} is an exponential distribution with rate 1, then their
#'survivor function is: \eqn{\exp (-y)}{exp(-y)}. This is what's shown on the
#'third graph generated by the \code{plot} method, using a log scale for the
#'ordinate. The point wise CI at 95 and 99\% are also drawn (dotted lines).
#'This is the graph shown on Fig. 12 (p 20) of Ogata (1988)
#'
#'If the \eqn{u_i}{u[i]} of the second paragraph are realizations of iid
#'uniform rvs on [0,1) then a plot of \eqn{u_{i+1}}{u[i+1]} vs \eqn{u_i}{u[i]}
#'should fill uniformly the unit square [0,1) x [0,1). This is the fourth
#'generated graph (the one shown on Fig. 11 (p 20) of Ogata (1988)) by the
#'\code{plot} method while the seventh graph shows the autocorrelation function
#'of the \eqn{u_i}{u[i]}s. Component \code{RenewalTest} of a
#'\code{CountingProcessSamplePath.summary} list contains a slightly more
#'elaborated (and quantitative) version of this test summarizing the fourth
#'graph (bottom right) generated by a call to \code{\link{renewalTestPlot}}.
#'This list element is itself a list with elements: \code{chi2.95} (a
#'\code{logical}), \code{chi2.99} (a \code{logical}) and \code{total} (an
#'\code{integer}). The bounds resulting from repetitively testing a sequence of
#'what are, under the null hypothesis, iid \eqn{\chi^2}{chi2} rvs are obtained
#'using Donsker's Theorem (see bellow). For each lag the number of degrees of
#'freedom of the \eqn{\chi^2}{chi2} distribution is subtracted from each
#'\eqn{\chi^2}{chi2} value. These centered values are then divided by their sd
#'(assuming the null hypothesis is correct). The cumulative sum of the centered
#'and scaled sequence is formed and is divided by the square root of the
#'maximal lag used. This is "plugged-in" the Donsker's Theorem. The eighth
#'graph of the \code{plot} method displays the resulting Wiener process. With
#'the tight confidence regions of Kendall et al (2007), see bellow.
#'
#'If the \eqn{x_i}{x[i]} are realization of a homogeneous Poisson process
#'observed between 0 and T, then the number of events observed on
#'non-overlapping windows of length t should be iid Poisson rv with mean t (and
#'variance t). The observation period is therefore chopped into non-overlapping
#'windows of increasing length and the empirical variance of the event count is
#'graphed versus the empirical mean, together with 95 and 99\% CI (using a
#'normal approximation). This is done by calling internally
#'\code{\link{varianceTime}}. That's what's generated by the fifth graph of the
#'\code{plot} method. This is the graph shown on Fig. 13 (p 20) of Ogata
#'(1988). Component \code{varianceTimeSummary} of a
#'\code{CountingProcessSamplePath.summary} list contains a summary of this
#'test, counting the number of events out of each band.
#'
#'The last graph generated by the \code{plot} method and the companions
#'components, \code{Wiener95} and \code{Wiener99}, of a
#'\code{CountingProcessSamplePath.summary} list represent "new" tests (as far
#'as I know). They are based on the fact that if the \eqn{y_i}{y[i]} above are
#'realizations of iid rvs following an exponential distribution with rate 1,
#'then the \eqn{w_i=y_i-1}{w[i]=y[i]-1} are realizations of iid rvs with mean 0
#'and variance 1. We can then form the partial sums:
#'\deqn{S_n=w_1+\cdots+w_n}{S[n]=w[1]+...+w[n]} and define the random right
#'continuous with a left-hand limit functions on [0,1]:
#'\deqn{\frac{1}{\sqrt{n}}S_{\lfloor nt \rfloor}}{S[floor(n*t)]/sqrt(n)} These
#'functions are realizations of a process which converges (weakly) to a Wiener
#'process on [0,1]. The proof of this statement is a corollary of Donsker's
#'Theorem and can be found on pp 146-147, Theorem 14.1, of Billingsley (1999).
#'I thank Vilmos Prokaj for pointing this reference to me.What is then done is
#'testing if the putative Wiener process is entirely within the tight
#'boundaries defined by Kendall et al (2007) for a true Wiener process, see
#'\code{\link{crossTight}}.
#'
#'@aliases summary.CountingProcessSamplePath
#'print.CountingProcessSamplePath.summary
#'plot.CountingProcessSamplePath.summary
#'@param object A \code{CountingProcessSamplePath} object.
#'@param exact Should an exact Kolmogorov test be used? See
#'\code{\link{ks.test}}.
#'@param lag.max See \code{\link{renewalTestPlot}}.
#'@param d See \code{\link{renewalTestPlot}}.
#'@param x A \code{CountingProcessSamplePath.summary} object.
#'@param digits An \code{integer}, the number of digits to be used while
#'printing summaries. See \code{\link{round}}.
#'@param y Not used but required for compatibility with the \code{\link{plot}}
#'method.
#'@param which If a subset of the test plots is required, specify a subset of
#'the numbers \code{1:6}.
#'@param main Title to appear above the plots, if missing the corresponding
#'element of \code{caption} will be used.
#'@param caption Default caption to appear above the plots or, if \code{main}
#'is given, bellow it
#'@param ask A \code{logical}; if \code{TRUE}, the user is \emph{ask}ed to hit
#'the return key before each plot generation, see \code{\link{par}(ask=.)}.
#'@param \dots Passed to \code{\link{chisq.test}} used internally by
#'\code{summary}, not used in \code{plot} and \code{print}.
#'@return \code{summary.CountingProcessSamplePath} returns a
#'\code{CountingProcessSamplePath.summary} object which is a \code{list} with
#'the following components:
#'@returnItem UniformGivenN A \code{numeric}, the p.value of the Kolmogorov
#'test of uniformity of the events times \emph{given} the number of events.
#'@returnItem Wiener95 A \code{logical}: is the scaled martingale within the
#'tight 95\% confidence band?
#'@returnItem Wiener99 A \code{logical}: is the scaled martingale within the
#'tight 99\% confidence band?
#'@returnItem BermanTest A \code{numeric}, the p.value of the Kolmogorov test
#'of uniformity of the scaled inter events intervals.
#'@returnItem RenewalTest A \code{list} with components: \code{chi2.95},
#'\code{chi2.99} and \code{total}. \code{chi2.95} resp. \code{chi2.99} is a
#'\code{logical} and is \code{TRUE} if the Wiener process obtained as described
#'above is within the "tight" 95\% resp. 99\% confidence band of Kendall et al
#'(2007).  \code{total} gives the total number of lags. See
#'\code{\link{renewalTestPlot}}.
#'@returnItem varianceTime A \code{\link{varianceTime}} object.
#'@returnItem varianceTimeSummary A \code{numeric} vector with components:
#'\code{total}, \code{out95} and \code{out99}. \code{total} gives the total
#'number of window sizes explored. \code{out95} gives the number of windows
#'within which the variance is out of the 95\% confidence band. \code{out99}
#'gives the number of windows within which the variance is out of the 99\%
#'confidence band. See \code{\link{varianceTime}}.
#'@returnItem n An \code{integer}, the number of events.
#'@returnItem call The matched call.
#'@note These functions / methods are designed to replace the
#'\code{summary.transformedTrain} and \code{plot.transformedTrain} ones. The
#'former have a more general design.
#'
#'Of course to be fully usable, these functions must be coupled to functions
#'allowing users to fit conditional intensity models.The support for that in
#'\code{STAR} is not complete yet but is coming soon. See for now the example
#'bellow.
#'
#'The end of the example bellow (not ran by default) shows that the coverage
#'probability of the Wiener Process confidence bands are really good even for
#'small (50) sample sizes.
#'@section Acknowledgments : I thank Vilmos Prokaj for pointing out Donsker's
#'Theorem and for indicating me the proof's location (Patrick Billingsley's
#'book).
#'
#'I also thank Olivier Faugeras and Jonathan Touboul for pointing out Donsker's
#'therom to me.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{mkCPSP}}, \code{\link{as.CPSP}},
#'\code{\link{plot.CountingProcessSamplePath}},
#'\code{\link{print.CountingProcessSamplePath}}, \code{\link{varianceTime}},
#'\code{\link{crossTight}}, \code{\link{renewalTestPlot}},
#'\code{\link{ks.test}}
#'@references Patrick Billingsley (1999) \emph{Convergence of Probability
#'Measures}. Wiley - Interscience.
#'
#'Brillinger, D. R. (1988) Maximum likelihood analysis of spike trains of
#'interacting nerve cells. \emph{Biol. Cybern.} \bold{59}: 189--200.
#'
#'Brown, E. N., Barbieri, R., Ventura, V., Kass, R. E. and Frank, L. M. (2002)
#'The time-rescaling theorem and its application to neural spike train data
#'analysis. \emph{Neural Computation} \bold{14}: 325-346.
#'
#'D. R. Cox and P. A. W. Lewis (1966) \emph{The Statistical Analysis of Series
#'of Events}. John Wiley and Sons.
#'
#'Daley, D. J. and Vere-Jones D. (2003) \emph{An Introduction to the Theory of
#'Point Processes. Vol. 1}. Springer.
#'
#'Ogata, Yosihiko (1988) Statistical Models for Earthquake Occurrences and
#'Residual Analysis for Point Processes. \emph{Journal of the American
#'Statistical Association} \bold{83}: 9-27.
#'
#'Johnson, D.H. (1996) Point process models of single-neuron discharges.
#'\emph{J. Computational Neuroscience} \bold{3}: 275--299.
#'
#'W. S. Kendall, J.- M. Marin and C. P. Robert (2007) Brownian Confidence Bands
#'on Monte Carlo Output. \emph{Statistics and Computing} \bold{17}: 1--10.
#'Preprint available at:
#'\url{http://www.ceremade.dauphine.fr/%7Exian/kmr04.rev.pdf}
#'@keywords distribution ts survival htest
#'@examples
#'
#'\dontrun{
#'## load one spike train data set of STAR
#'data(e060824spont)
#'## Create the CountingProcessSamplePath object
#'n1spt.cp <- as.CPSP(e060824spont[["neuron 1"]])
#'## print it
#'n1spt.cp
#'## plot it
#'plot(n1spt.cp)
#'## get the summary
#'## Notice the warning due to few identical interspike intervals
#'## leading to an inaccurate Berman's test.
#'summary(n1spt.cp)
#'
#'## Simulate data corresponding to a renewal process with
#'## an inverse Gaussian ISI distribution in the spontaneous
#'## regime modulated by a multiplicative stimulus whose time
#'## course is a shifted and scaled chi2 density.
#'## Define the "stimulus" function
#'stimulus <- function(t,
#'                     df=5,
#'                     tonset=5,
#'                     timeFactor=5,
#'                     peakFactor=10) {
#'  dchisq((t-tonset)*timeFactor,df=df)*peakFactor
#'}
#'## Define the conditional intensity / hazard function
#'hFct <- function(t,
#'                 tlast,
#'                 df=5,
#'                 tonset=5,
#'                 timeFactor=5,
#'                 peakFactor=10,
#'                 mu=0.075,
#'                 sigma2=3
#'                 ) {
#'  
#'  hinvgauss(t-tlast,mu=mu,sigma2=sigma2)*exp(stimulus(t,df,tonset,timeFactor,peakFactor))
#'
#'}
#'## define the function simulating the train with the thinning method                   
#'makeTrain <- function(tstop=10,
#'                      peakCI=200,
#'                      preTime=5,
#'                      df=5,
#'                      tonset=5,
#'                      timeFactor=5,
#'                      peakFactor=10,
#'                      mu=0.075,
#'                      sigma2=3) {
#'
#'  result <- numeric(500) - preTime - .Machine$double.eps
#'  result.n <- 500
#'  result[1] <- 0
#'  idx <- 1
#'  currentTime <- result[1]
#'  while (currentTime < tstop+preTime) {
#'    currentTime <- currentTime+rexp(1,peakCI)
#'    p <- hFct(currentTime,
#'              result[idx],
#'              df=df,
#'              tonset=tonset+preTime,
#'              timeFactor=timeFactor,
#'              peakFactor=peakFactor,
#'              mu=mu,
#'              sigma2=sigma2)/peakCI
#'    rthreshold <- runif(1)
#'    if (p>1) stop("Wrong peakCI")
#'    while(p < rthreshold) {
#'      currentTime <- currentTime+rexp(1,peakCI)
#'      p <- hFct(currentTime,
#'                result[idx],
#'                df=df,
#'                tonset=tonset+preTime,
#'                timeFactor=timeFactor,
#'                peakFactor=peakFactor,
#'                mu=mu,
#'                sigma2=sigma2)/peakCI
#'      if (p>1) stop("Wrong peakCI")
#'      rthreshold <- runif(1)
#'    }
#'    idx <- idx+1
#'    if (idx > result.n) {
#'      result <- c(result,numeric(500)) - preTime - .Machine$double.eps
#'      result.n <- result.n + 500
#'    }
#'    result[idx] <- currentTime
#'  }
#'
#'  result[preTime < result & result <= tstop+preTime] - preTime
#'  
#'}
#'## set the seed
#'set.seed(20061001)
#'## "make" the train
#'t1 <- makeTrain()
#'## create the corresponding CountingProcessSamplePath
#'## object
#'cpsp1 <- mkCPSP(t1)
#'## print it
#'cpsp1
#'## test it
#'cpsp1.summary <- summary(cpsp1)
#'cpsp1.summary
#'plot(cpsp1.summary)
#'## Define a function returning the conditional intensity function (cif)
#'ciFct <- function(t,
#'                  tlast,
#'                  df=5,
#'                  tonset=5,
#'                  timeFactor=5,
#'                  peakFactor=10,
#'                  mu=0.075,
#'                  sigma2=3
#'                  ) {
#'
#'  sapply(t, function(x) {
#'    if (x <= tlast[1]) return(1/mu)
#'    y <- x-max(tlast[tlast<x])
#'    hinvgauss(y,mu=mu,sigma2=sigma2)*exp(stimulus(x,df,tonset,timeFactor,peakFactor))
#'  }
#'         )
#'
#'}
#'## Compute the cif of the train
#'tt <- seq(0,10,0.001)
#'lambda.true <- ciFct(tt,cpsp1$ppspFct())
#'## plot it together with the events times
#'## Notice that the representation is somewhat inaccurate, the cif
#'## is in fact a left continuous function
#'plot(tt,lambda.true,type="l",col=2)
#'rug(cpsp1$ppspFct())
#'## plot the integrated intensity function and the counting process
#'plot(tt,cumsum(lambda.true)*0.001,type="l",col=2)
#'lines(cpsp1)
#'## define a function doing the time transformation / rescaling
#'## by integrating the cif and returning another CountingProcessSamplePath
#'transformCPSP <- function(cpsp,
#'                          ciFct,
#'                          CIFct,
#'                          method=c("integrate","discrete"),
#'                          subdivisions=100,
#'                          ...
#'                          ) {
#'
#'  if (!inherits(cpsp,"CountingProcessSamplePath"))
#'    stop("cpsp should be a CountingProcessSamplePath objet")
#'  st <- cpsp$ppspFct()
#'  n <- length(st)
#'  from <- cpsp$from
#'  to <- cpsp$to
#'  if (missing(CIFct)) {
#'    if (method[1] == "integrate") {
#'      lwr <- c(from,st)
#'      upr <- c(st,to)
#'      Lambda <- sapply(1:(n+1),
#'                       function(idx)
#'                       integrate(ciFct,
#'                                 lower=lwr[idx],
#'                                 upper=upr[idx],
#'                                 subdivisions=subdivisions,
#'                                 ...)$value
#'                       )
#'      Lambda <- cumsum(Lambda)
#'      st <- Lambda[1:n]
#'      from <- 0
#'      to <- Lambda[n+1]
#'    } ## End of conditional on method[1] == "integrate"
#'    if (method[1] == "discrete") {
#'      lwr <- c(from,st)
#'      upr <- c(st,to)
#'      xx <- unlist(lapply(1:(n+1),
#'                          function(idx) seq(lwr[idx],
#'                                            upr[idx],
#'                                            length.out=subdivisions)
#'                          )
#'                   )
#'      Lambda <- cumsum(ciFct(xx[-length(xx)])*diff(xx))
#'      Lambda <- Lambda - Lambda[1]
#'      st <- Lambda[(1:n)*subdivisions]
#'      from <- 0
#'      to <- Lambda[length(Lambda)]
#'    } ## End of conditional on method[1] == "discrete"
#'  } else {
#'    result <- CIFct(c(from,st,to))
#'    result <- result-result[1]
#'    from <- result[1]
#'    to <- result[n+2]
#'    st <- result[2:(n+1)]
#'  } ## End of conditional on missing(CIFct)
#'  mkCPSP(st,from,to)
#'}
#'## transform cpsp1
#'cpsp1t <- transformCPSP(cpsp1,function(t) ciFct(t,cpsp1$ppspFct()))
#'## test it
#'cpsp1t.summary <- summary(cpsp1t)
#'cpsp1t.summary
#'plot(cpsp1t.summary)
#'## compare the finite sample performances of the
#'## Kolmogorov test (test the uniformity of the
#'## jump times given the number of events) with the
#'## ones of the new "Wiener process test"
#'empiricalCovProb <- function(myRates=c(10,(1:8)*25,(5:10)*50,(6:10)*100),
#'                             nbRep=1000,
#'                             exact=NULL
#'                             ) {
#'
#'  b95 <- function(t) 0.299944595870772 + 2.34797018726827*sqrt(t)
#'  b99 <- function(t) 0.313071417065285 + 2.88963206734397*sqrt(t)
#'  result <- matrix(numeric(4*length(myRates)),nrow=4)
#'  colnames(result) <- paste(myRates)
#'  rownames(result) <- c("ks95","ks99","wp95","wp99")
#'  for (i in 1:length(myRates)) {
#'    rate <- myRates[i]
#'    partial <- sapply(1:nbRep,
#'                      function(repIdx) {
#'                        st <- cumsum(rexp(5*rate,rate))
#'                        while(max(st) < 1) st <- c(st,max(st)+cumsum(rexp(5*rate,rate)))
#'                        st <- st[st<=1]
#'                        ks <- ks.test(st,punif,exact=exact)$p.value
#'                        w <- (st*rate-seq(st))/sqrt(rate)
#'                        c(ks95=0.95 < ks,
#'                          ks99=0.99 < ks,
#'                          wp95=any(w < -b95(st) | b95(st) < w),
#'                          wp99=any(w < -b99(st) | b99(st) < w)
#'                          )
#'                      }
#'                      )
#'    
#'    result[,i] <- apply(partial,1,sum)
#'  }
#' 
#'  attr(result,"nbRep") <- nbRep
#'  attr(result,"myRates") <- myRates
#'  attr(result,"call") <- match.call()
#'  result/nbRep
#'  
#'}
#'
#'plotCovProb <- function(covprob,ci=0.95) {
#'
#'  nbMax <- max(attr(covprob,"myRates"))
#'  plot(c(0,nbMax),c(0.94,1),
#'       type="n",
#'       xlab="Expected number of Spikes",
#'       ylab="Empirical cov. prob.",xaxs="i",yaxs="i")
#'  nbRep <- attr(covprob,"nbRep")
#'  polygon(c(0,nbMax,nbMax,0),
#'          c(rep(qbinom((1-ci)/2,nbRep,0.95)/nbRep,2),rep(qbinom(1-(1-ci)/2,nbRep,0.95)/nbRep,2)),
#'          col="grey50",border=NA)
#'  polygon(c(0,nbMax,nbMax,0),
#'          c(rep(qbinom((1-ci)/2,nbRep,0.99)/nbRep,2),rep(qbinom(1-(1-ci)/2,nbRep,0.99)/nbRep,2)),
#'          col="grey50",border=NA)
#'  nbS <- attr(covprob,"myRates")
#'  points(nbS,1-covprob[1,],pch=3)
#'  points(nbS,1-covprob[2,],pch=3)
#'  points(nbS,1-covprob[3,],pch=1)
#'  points(nbS,1-covprob[4,],pch=1)
#'
#'}
#'system.time(covprobA <- empiricalCovProb())
#'plotCovProb(covprobA)
#'}
#'
summary.CountingProcessSamplePath <- function(object,
                                              exact = TRUE,
                                              lag.max = NULL,
                                              d = max(c(2, sqrt(length(object$ppspFct()))%/%5)),
                                              ...) {

  x <- object$ppspFct()
  n <- length(x)
  xn <- (x-object$from)/(x[n]-object$from)
  UniformGivenN <- ks.test(xn[-n],punif,exact=exact)$p.value

  b95 <- function(t) 0.299944595870772 + 2.34797018726827*sqrt(t)
  b99 <- function(t) 0.313071417065285 + 2.88963206734397*sqrt(t)

  ## Modification 2008 10 20 by C Pouzat
  ## y.w <- (seq(x) - x)/sqrt(object$to-object$from)
  y.w <- (x-x[1])[-1]
  ny <- length(y.w)
  y.w <- (y.w-1:ny)/sqrt(ny)
  x.w <- (1:ny)/ny

  Wiener95 <- all(-b95(x.w) < y.w & y.w < b95(x.w))
  Wiener99 <- all(-b99(x.w) < y.w & y.w < b99(x.w))

  xs <- sort(diff(x))
  BermanTest <- ks.test(xs,pexp,exact=exact)$p.value

  isi <- diff(x)
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
    for (i in seq(along.with = isi.x)) counts[isi.x[i], isi.y[i]] <- counts[isi.x[i], 
                    isi.y[i]] + 1
    chisq.test(counts, ...)
  }
  chi2seq <- lapply(1:lag.max, getChi2)
  chi2.df <- chi2seq[[1]]$parameter
  ## minChi2 <- qchisq(0.025, df = chi2.df)
  ## maxChi2 <- qchisq(0.975, df = chi2.df)
  chi2V <- sapply(chi2seq, function(l) l$statistic)
  ## outOf95 <- chi2V < minChi2 | chi2V > maxChi2
  chi2.w <- cumsum((chi2V-chi2.df)/sqrt(2*chi2.df*lag.max))
  chi2.95 <- all(abs(chi2.w) < b95((1:lag.max)/lag.max))
  chi2.99 <- all(abs(chi2.w) < b99((1:lag.max)/lag.max))
                 
  vt <- varianceTime(x,CI=c(0.95,0.99))
  if (!is.null(dim(vt$ciLow))) {
    vt.out95 <- vt$s2 < vt$ciLow[1,] | vt$ciUp[1,] < vt$s2
    vt.out99 <- vt$s2 < vt$ciLow[2,] | vt$ciUp[2,] < vt$s2
  } else {
    vt.out95 <- vt$s2 < vt$ciLow[1] | vt$ciUp[1] < vt$s2
    vt.out99 <- vt$s2 < vt$ciLow[2] | vt$ciUp[2] < vt$s2
    warning("Only one window size for varianceTime")
  }
  
  result <- list(UniformGivenN=UniformGivenN,
                 Wiener95=Wiener95,
                 Wiener99=Wiener99,
                 BermanTest=BermanTest,
                 RenewalTest=list(chi2.95=chi2.95,chi2.99=chi2.99,total=lag.max),


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
#'@param obj a object to test against a \code{varianceTime} object.
#'@param x a \code{varianceTime} object.
#'@param CI a numeric vector with at most two elements. The coverage
#'probability of the confidence intervals.
#'@param windowSizes a numeric increasing vector of positive numbers. The
#'window sizes used to split the spike train.
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
                 varianceTime=vt,
                 varianceTimeSummary=c(total=length(vt$s2),out95=sum(vt.out95),out99=sum(vt.out99)),
                 n=n,
                 call=match.call()
                 )

  class(result) <- "CountingProcessSamplePath.summary"
  result
  
}

print.CountingProcessSamplePath.summary <- function(x,
                                                    digits=5,
                                                    ...) {

  cat(paste(" *** Test of uniformity on the time axis \n",
            "    Prob. of the Kolmogorov statistic under H0:",
            round(x$UniformGivenN,digits=digits),"\n"
            )
      )
  cat(paste(" *** Wiener process test \n",
            "    Inside 95% domain:",x$Wiener95,", inside 99% domain:",x$Wiener99,"\n"
            )
      )
  cat(paste(" *** Berman test  \n",
            "    Prob. of the Kolmogorov statistic under H0:",
            round(x$BermanTest,digits=digits),"\n"
            )
      )
  cat(paste(" *** Renewal test \n",
            "    Inside 95% domain:",
            x$RenewalTest[["chi2.95"]],
            ", inside 99% domain:",
            x$RenewalTest[["chi2.99"]],"\n",
            "    Maximum lag:", x$RenewalTest[["total"]],"\n"
            )
      )
  ## if (x$RenewalTest["out"]>1)
  ##   cat(paste("    ",x$RenewalTest["out"],"lags out on a total of:",x$RenewalTest["total"],"\n"))
  ## else
  ##  cat(paste("    ",x$RenewalTest["out"],"lag out on a total of:",x$RenewalTest["total"],"\n"))

  total <- x$varianceTimeSummary["total"]
  out95 <- x$varianceTimeSummary["out95"]
  out99 <- x$varianceTimeSummary["out99"]
  cat(paste(" *** Variance vs \"time\" with", total,"time windows:\n"))
  if (out95>1) cat(paste("    ",out95,"windows out at 95% level\n"))
  else cat(paste("    ",out95,"window out at 95% level\n"))
  if (out99>1) cat(paste("    ",out99,"windows out at 99% level\n"))
  else cat(paste("    ",out99,"window out at 99% level\n"))
  cat(paste(" *** The object contains", x$n,"events.\n"))

}


plot.CountingProcessSamplePath.summary <- function (x, y,
                                                    which = c(1,2,6,8),
                                                    main,
                                                    caption = c(
                                                      expression(paste("Uniform on ", Lambda," Test")), 
                                                      "Berman's Test",
                                                      "Log Survivor Function",
                                                      expression(paste(U[k+1]," vs ", U[k])), 
                                                      "Variance vs Mean Test",
                                                      "Wiener Process Test",
                                                      "Autocorrelation Fct.",
                                                      "Renewal Test"),
                                                    ask = FALSE,
                                                    lag.max = NULL,
                                                    d = max(c(2, sqrt(length(eval(x$call[[2]])$ppspFct()))%/%5)),
                                                    ...) {
  
    if (!inherits(x, "CountingProcessSamplePath.summary")) 
        stop("use only with \"CountingProcessSamplePath.summary\" objects")

    force(d)
    ## Get the spike times
    cpsp <- eval(x$call[[2]])
    xx <- cpsp$ppspFct()
    isi <- diff(xx)
    Y <- seq(xx)
    nbSpikes <- length(xx)
    b95 <- function(t) 0.299944595870772 + 2.34797018726827*sqrt(t)
    b99 <- function(t) 0.313071417065285 + 2.88963206734397*sqrt(t)
      
    show <- logical(8)
    show[which] <- TRUE
    if (ask) {
        op <- par(ask = TRUE)
        on.exit(par(op))
    }
    else {
        if (sum(show) == 2) 
            layout(matrix(1:2, nrow = 2))
        if (2 < sum(show) && sum(show) < 5) 
            layout(matrix(1:4, nrow = 2))
        if (sum(show) >= 5) 
            layout(matrix(1:9, nrow = 3))
    }
    mainGiven <- !missing(main)
    if (show[1]) {
      slopeKS <- length(xx)/max(xx)
      plot(as.numeric(xx), Y, type = "n", xlab = expression(Lambda), 
           ylab = expression(N(Lambda)),
           main = ifelse(mainGiven, main, caption[1]),
           sub = ifelse(mainGiven,caption[1],"")
           )
        abline(a = 0, b = slopeKS)
        abline(a = 1.36 * sqrt(nbSpikes), slopeKS, lty = 2)
        abline(a = -1.36 * sqrt(nbSpikes), slopeKS, lty = 2)
        abline(a = 1.63 * sqrt(nbSpikes), slopeKS, lty = 2)
        abline(a = -1.63 * sqrt(nbSpikes), slopeKS, lty = 2)
        lines(as.numeric(xx), Y, col = 2, lwd = 2)
    }
    lambda <- 1 - exp(-isi)
    if (show[2]) {
        plot(c(0, 1), c(0, 1), type = "n", xlab = expression(U[(k)]), 
            ylab = "ECDF", main = ifelse(mainGiven, 
                main, caption[2]), sub = ifelse(mainGiven, caption[2], 
                ""))
        abline(a = 0, b = 1)
        abline(a = 1.36/sqrt(nbSpikes - 1), 1, lty = 2)
        abline(a = -1.36/sqrt(nbSpikes - 1), 1, lty = 2)
        abline(a = 1.63/sqrt(nbSpikes - 1), 1, lty = 2)
        abline(a = -1.63/sqrt(nbSpikes - 1), 1, lty = 2)
        lines(sort(lambda), (1:(nbSpikes - 1))/(nbSpikes - 1), 
            col = 2, lwd = 2)
    }
    if (show[3]) {
        nI <- length(isi)
        Y <- (nI:1)/nI
        X <- sort(isi)
        Yth <- exp(-X)
        Y95p <- qbinom(0.975, nI, Yth)/nI
        Y95m <- qbinom(0.025, nI, Yth)/nI
        Y99p <- qbinom(0.995, nI, Yth)/nI
        Y99m <- qbinom(0.005, nI, Yth)/nI
        maxId <- max(which(Yth > 0.001))
        plot(c(0, X[maxId]), c(0.001, 1), type = "n", xlab = expression(Y[(k)]), 
            ylab = "Survivor Fct", main = ifelse(mainGiven, main, 
                caption[3]), sub = ifelse(mainGiven, caption[3], 
                ""), log = "y")
        lines(X, Y95p, lty = 2)
        lines(X, Y95m, lty = 2)
        lines(X, Y99p, lty = 2)
        lines(X, Y99m, lty = 2)
        lines(X, Y, col = 2, lwd = 2)
    }
    if (show[4]) {
        plot(lambda[-length(lambda)], lambda[-1], xlab = expression(U[k]), 
            ylab = expression(U[k + 1]), pch = 3, main = ifelse(mainGiven, 
                main, caption[4]), sub = ifelse(mainGiven, caption[4], 
                ""))
    }
    if (show[5]) {
        plot(x$varianceTime, style = "Ogata", xlab = "Window Length", 
            main = ifelse(mainGiven, main, caption[5]), sub = ifelse(mainGiven, 
                caption[5], ""))
    }
    if (show[6]) {
      y.w <- (xx-xx[1])[-1]
      n <- length(y.w)
      y.w <- c(0,(y.w-1:n)/sqrt(n))
      x.w <- (0:n)/n
      tt <- seq(0,1,0.01)
      ylim <- c(min(-b99(tt[length(tt)]),min(y.w)),
                max(b99(tt[length(tt)]),max(y.w))
                )
      plot(x.w, y.w, type = "n", xlab = expression(t),
           ylab = expression(X[t]^n), 
           ylim = ylim,
           main = ifelse(mainGiven, main, caption[6]),
           sub = ifelse(mainGiven, caption[6], ""),
           xaxs="i",yaxs="i"
           )
        abline(h = 0)
        lines(tt, b95(tt), lty = 2)
        lines(tt, -b95(tt), lty = 2)
        lines(tt, b99(tt), lty = 2)
        lines(tt, -b99(tt), lty = 2)
        lines(x.w, y.w, col = 2, lwd = 2)
    }
    if (show[7]) {
      acf(isi,lag.max=lag.max,
          main = ifelse(mainGiven, main, caption[7]),
           sub = ifelse(mainGiven, caption[7], ""),
          ...)
    }
    if (show[8]) {
      isi <- diff(xx)
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
        for (i in seq(along.with = isi.x)) counts[isi.x[i], isi.y[i]] <- counts[isi.x[i], 
                        isi.y[i]] + 1
        chisq.test(counts, ...)
      }
      chi2seq <- lapply(1:lag.max, getChi2)
      chi2.df <- chi2seq[[1]]$parameter
      chi2V <- sapply(chi2seq, function(l) l$statistic)
      y.w <- c(0,cumsum((chi2V-chi2.df)/sqrt(2*chi2.df*lag.max)))
      x.w <- (0:lag.max)/lag.max
      tt <- seq(0,1,0.01)
      ylim <- c(min(-b99(tt[length(tt)]),min(y.w)),
                max(b99(tt[length(tt)]),max(y.w))
                )
      sub <- ifelse(mainGiven,
                    paste(caption[8],", lag.max = ",lag.max,", d =",d,sep=""),
                    paste("lag.max = ",lag.max,", d = ",d,sep="")
                    )
      plot(x.w, y.w, type = "n", xlab = expression(t),
           ylab = expression(X[t]^n), 
           ylim = ylim,
           main = ifelse(mainGiven, main, caption[8]),
           sub = sub,
           xaxs="i",yaxs="i"
           )
        abline(h = 0)
        lines(tt, b95(tt), lty = 2)
        lines(tt, -b95(tt), lty = 2)
        lines(tt, b99(tt), lty = 2)
        lines(tt, -b99(tt), lty = 2)
        lines(x.w, y.w, col = 2, lwd = 2)
      
    }

}
