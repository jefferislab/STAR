#'Plot Diagnostics for an transformedTrain Object
#'
#'Six plots (selectable by \code{which}) are currently available: the first 5
#'of which correspond to Fig. 9 to 13 of Ogata (1988). The sixth one is new (as
#'far as I know) and is still "experimental". They are all testing the first
#'argument of \code{plot.transformedTrain} against the Poisson process
#'hypothesis..
#'
#'If the \code{\link{transformedTrain}} object \code{x} is a the realization of
#'a homogeneous Poisson process then, conditioned on the number of events
#'observed, the location of the events is uniform on the (time transformed)
#'period of observation. This is a basic property of the homogeneous Poisson
#'process derived in Chap. 2 of Cox and Lewis (1966) and Daley and Vere-Jones
#'(2003). This is what the first plot generated (by default) tests with a
#'Kolmogorov-Smirnov Test. The two dotted lines on both sides of the diagonal
#'correspond to 95 and 99\% confidence intervals. This is the plot shown on
#'Fig. 9 (p 19) of Ogata (1988).
#'
#'If we write \eqn{x_i}{x[i]} the elements of the
#'\code{\link{transformedTrain}} object \code{x} and if the latter is the
#'realization of a homogeneous Poisson process then the intervals: \deqn{y_i =
#'x_{i+1}-x_i}{y[i]=x[i+1]-x[i]} are iid rv from an exponential distribution
#'with rate 1 and the: \deqn{u_i = 1-\exp (-y_i}{u[i]=1 - exp(-y[i])} are iid
#'rv from a uniform distribution on [0,1). The second plot generated (by
#'default) tests this uniform distribution hypotheses with a Kolmogorov-Smirnov
#'Test. This is the plot shown on Fig. 10 (p 19) of Ogata (1988) which was
#'suggested by Berman. This is also the plot proposed by Brown et al (2002).
#'The two dotted lines on both sides of the diagonal correspond to 95 and 99\%
#'confidence intervals.
#'
#'Following the line of the previous paragraph, if the distribution of the
#'\eqn{y_i}{y[i]} is an exponential distribution with rate 1, then their
#'survivor function is: \eqn{\exp (-y)}{exp(-y)}. This is what's shown on the
#'third plot generated (by default) using a log scale for the ordinate. The
#'point wise CI at 95 and 99\% are also drawn (dotted lines). This is the plot
#'shown on Fig. 12 (p 20) of Ogata (1988)
#'
#'If the \eqn{u_i}{u[i]} of the second paragraph are iid uniform rv on [0,1)
#'then a plot of \eqn{u_{i+1}}{u[i+1]} vs \eqn{u_i}{u[i]} should fill uniformly
#'the unit square [0,1) x [0,1). This is the fourth generated plot (by
#'default). This is the plot shown on Fig. 11 (p 20) of Ogata (1988)
#'
#'If the \eqn{x_i}{x[i]} are realization of a homogeneous Poisson process
#'observed between 0 and T (on the transformed time scale), then the number of
#'events observed on non-overlapping windows of length t should be iid Poisson
#'rv with mean t (and variance t). The observation period is chopped into
#'non-overlapping windows of increasing length and the empirical variance of
#'the event count is plotted versus the empirical mean, together with 95 and
#'99\% CI (using a normal approximation). This is done by calling internally
#'\code{\link{varianceTime}}. That's what's generated by the fifth plot (by
#'default). This is the plot shown on Fig. 13 (p 20) of Ogata (1988)
#'
#'The last plot is experimental and irrelevant for spike trains transformed
#'after a \code{\link{gam}} or a \code{\link{glm}} fit. It should be useful for
#'parametric models fitted with the maximum likelihood method.
#'
#'@param x a \code{\link{transformedTrain}} object.
#'@param which if a subset of the plots is required, specify a subset of the
#'numbers \code{1:6}.
#'@param main title to appear above the plots, if missing the corresponding
#'element of \code{caption} will be used.
#'@param caption Default caption to appear above the plots or, if \code{main}
#'is given, bellow it
#'@param ask logical; if \code{TRUE}, the user is \emph{ask}ed before each
#'plot, see \code{\link{par}(ask=.)}.
#'@param \dots not used only there for compatibility with \code{\link{plot}}
#'generic method.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{transformedTrain}},
#'\code{\link{summary.transformedTrain}}, \code{\link{mkGLMdf}}
#'@references Cox, D. R. and Lewis, P. A. W. (1966) \emph{The Statistical
#'Analysis of Series of Events}. John Wiley and Sons.
#'
#'Daley, D. J. and Vere-Jones D. (2003) \emph{An Introduction to the Theory of
#'Point Processes. Vol. 1}. Springer.
#'
#'Ogata, Yosihiko (1988) Statistical Models for Earthquake Occurrences and
#'Residual Analysis for Point Processes. \emph{Journal of the American
#'Statistical Association} \bold{83}: 9-27.
#'
#'Brown, E. N., Barbieri, R., Ventura, V., Kass, R. E. and Frank, L. M. (2002)
#'The time-rescaling theorem and its application to neural spike train data
#'analysis. \emph{Neural Computation} \bold{14}: 325-346.
#'@keywords models smooth regression hplot
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
plot.transformedTrain <- function(x,
                                  which=1:5,
                                  main,
                                  caption=c("Uniform on Trans. Obs. Time Test",
                                    "Berman's Test",
                                    "Log Survivor Function",
                                    "Lag 1 Transformed Intervals",
                                    "Variance vs Mean",
                                    "Martingale vs Trans. Time"),
                                  ask=TRUE,
                                  ...
                                  ) {

  ## Check that x is a transformedTrain
  if (!inherits(x, "transformedTrain"))
    stop("use only with \"transformedTrain\" objects")

  Y <- seq(x)
  nbSpikes <- length(x)
  
  show <- logical(6)
  show[which] <- TRUE

  if (ask) {
    op <- par(ask = TRUE)
    on.exit(par(op))
  } else {
    if (sum(show)==2) layout(matrix(1:2,nrow=2))
    if (2<sum(show) && sum(show)<5) layout(matrix(1:4,nrow=2))
    if (sum(show)>=5) layout(matrix(1:6,nrow=3))
  }

  mainGiven <- !missing(main)
  if (show[1]) {
    slopeKS <- length(x)/max(x)
    plot(as.numeric(x),Y,type="n",
         xlab="Transformed time",
         ylab="Cumulative number of events",
         main=ifelse(mainGiven,main,caption[1]),
         sub=ifelse(mainGiven,caption[1],"")
         )
    abline(a=0,b=slopeKS)
    abline(a=1.36*sqrt(nbSpikes),slopeKS,lty=2)
    abline(a=-1.36*sqrt(nbSpikes),slopeKS,lty=2)
    abline(a=1.63*sqrt(nbSpikes),slopeKS,lty=2)
    abline(a=-1.63*sqrt(nbSpikes),slopeKS,lty=2)
    lines(as.numeric(x),Y,col=2,lwd=2)
  }

  lambda <- 1-exp(-diff(x))
  if (show[2]) {
    plot(c(0,1),c(0,1),type="n",
         xlab=expression(U[(k)]),
         ylab="Cumulative distribution",
         main=ifelse(mainGiven,main,caption[2]),
         sub=ifelse(mainGiven,caption[2],"")
         )
    abline(a=0,b=1)
    abline(a=1.36/sqrt(nbSpikes-1),1,lty=2)
    abline(a=-1.36/sqrt(nbSpikes-1),1,lty=2)
    abline(a=1.63/sqrt(nbSpikes-1),1,lty=2)
    abline(a=-1.63/sqrt(nbSpikes-1),1,lty=2)
    lines(sort(lambda),(1:(nbSpikes-1))/(nbSpikes-1),col=2,lwd=2)
  }

  if (show[3]) {
    isi <- diff(x)
    nI <- length(isi)
    Y <- (nI:1)/nI
    X <- sort(isi)
    Yth <- exp(-X)
    Y95p <- qbinom(0.975,nI,Yth)/nI
    Y95m <- qbinom(0.025,nI,Yth)/nI
    Y99p <- qbinom(0.995,nI,Yth)/nI
    Y99m <- qbinom(0.005,nI,Yth)/nI
    maxId <- max(which(Yth>0.001))
    plot(c(0,X[maxId]),c(0.001,1),type="n",
         xlab=expression(Y[(k)]),
         ylab="Survivor Fct",
         main=ifelse(mainGiven,main,caption[3]),
         sub=ifelse(mainGiven,caption[3],""),
         log="y"
         )
    lines(X,Y95p,lty=2)
    lines(X,Y95m,lty=2)
    lines(X,Y99p,lty=2)
    lines(X,Y99m,lty=2)
    lines(X,Y,col=2,lwd=2)
  }

  if (show[4]) {
    plot(lambda[-length(lambda)],
         lambda[-1],
         xlab=expression(U[k]),
         ylab=expression(U[k+1]),
         pch=3,
         main=ifelse(mainGiven,main,caption[4]),
         sub=ifelse(mainGiven,caption[4],"")
         )
  }

  if (show[5]) {
    plot(varianceTime(x),
         style="Ogata",
         xlab="Length of Trans. Time Window",
         main=ifelse(mainGiven,main,caption[5]),
         sub=ifelse(mainGiven,caption[5],"")
         )
  }

  if (show[6]) {
    X <- c(0,as.numeric(x))
    M <- 0:(length(X)-1) - X
    mTh <- -1.63*sqrt(length(X)-1)
    MTh <- 1.63*sqrt(length(X)-1)
    plot(X,
         M,
         type="n",
         xlab="Transformed time",
         ylab="Martingale",
         ylim=c(min(c(M,mTh)),max(c(M,MTh))),
         main=ifelse(mainGiven,main,caption[6]),
         sub=ifelse(mainGiven,caption[6],"")
         )
    abline(h=0)
    abline(h=1.36*sqrt(length(X)-1),lty=2)
    abline(h=-1.36*sqrt(length(X)-1),lty=2)
    abline(h=1.63*sqrt(length(X)-1),lty=2)
    abline(h=-1.63*sqrt(length(X)-1),lty=2)
    lines(X,M,col=2,lwd=2)
  }
}
