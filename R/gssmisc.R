#'Time Transformation Using a gssanova Objet
#'
#'Performs time transformation using a \code{gssanova} fit. If the model is
#'correct, the result of the transformation should be a Poisson process with
#'rate 1.
#'
#'The binary operator applies \code{\link[gss]{predict.ssanova}} with the left
#'side as the first argument and the right side as the second argument. The
#'right side (\code{dataFrame}) must therefore contain the variables included
#'in the \code{formula} used in the call giving rise to \code{gssObj}. The
#'result of the \code{predict} method call is then transformed with an inverse
#'logistic function or with an exponential (depending on the \code{family}
#'argument, \code{"binomial"} or \code{"poisson"}, used in the previous
#'\code{\link[gss]{gssanova}} call). The cumulative sum is computed, that is,
#'the integrated conditional intensity, and its value at the events times is
#'returned as a \code{CountingProcessSamplePath} object.
#'
#'@param gssObj a \code{\link[gss]{gssanova}} or a \code{\link[gss]{gssanova0}}
#'object.
#'@param dataFrame a \code{data.frame} with variables corresponding to the ones
#'used in the \code{\link[gss]{gssanova}} call giving rise to \code{gssObj}.
#'@return A \code{CountingProcessSamplePath} object.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link[gss]{gssanova}}, \code{\link[gss]{predict.ssanova}},
#'\code{\link{mkGLMdf}}, \code{\link{mkCPSP}},
#'\code{\link{summary.CountingProcessSamplePath}}
#'@references Gu C. (2002) \emph{Smoothing Spline ANOVA Models}. Springer.
#'
#'Brillinger, D. R. (1988) Maximum likelihood analysis of spike trains of
#'interacting nerve cells. \emph{Biol. Cybern.} \bold{59}: 189--200.
#'
#'Brown, E. N., Barbieri, R., Ventura, V., Kass, R. E. and Frank, L. M. (2002)
#'The time-rescaling theorem and its application to neural spike train data
#'analysis. \emph{Neural Computation} \bold{14}: 325-346.
#'
#'Ogata, Yosihiko (1988) Statistical Models for Earthquake Occurrences and
#'Residual Analysis for Point Processes. \emph{Journal of the American
#'Statistical Association} \bold{83}: 9-27.
#'@keywords ts survival htest
#'@examples
#'
#'\dontrun{
#'## load e060517spont data set
#'data(e060517spont)
#'## make a data frame using a 2 ms bin width
#'e060517spontDF <- mkGLMdf(e060517spont,0.002,0,60)
#'## Keep data relevant to neuron 3
#'e060517spontDFn3 <- e060517spontDF[e060517spontDF$neuron == "3",]
#'## Split data in an "early" and a "late" part
#'e060517spontDFn3e <- e060517spontDFn3[e060517spontDFn3$time <= 30,]
#'e060517spontDFn3l <- e060517spontDFn3[e060517spontDFn3$time > 30,]
#'## fit the late part with a nonparametric renewal model
#'e060517spontDFn3lGF <- gssanova(event ~ lN.3, data=e060517spontDFn3l,family="binomial")
#'## transform the time of the early part
#'e060517spont.n3e.tt <- e060517spontDFn3lGF %tt% e060517spontDFn3e
#'## Test the goodness of fit
#'e060517spont.n3e.tt
#'summary(e060517spont.n3e.tt)
#'plot(summary(e060517spont.n3e.tt))
#'}
#'
"%tt%" <- function(gssObj,
                   dataFrame
                   ) {



  ## check the class of gssObj
  if (!inherits(gssObj,c("gssanova","gssanova0")))
    stop("gssObj should be a gssanova or a gssanova0 object.")

  ## check the dataFrame has the right variable
  if (!all(names(gssObj$mf) %in% names(dataFrame)))
    stop("dataFrame does not contain the appropriate variables.")

  ## make sure the variables of dataFrame are in the right range
  for (vn in names(gssObj$mf)[!(names(gssObj$mf) %in% "event")]) {
    m <- min(gssObj$mf[[vn]])
    M <- max(gssObj$mf[[vn]])
    dataFrame[[vn]] <- sapply(dataFrame[[vn]], function(x) min(x,M))
    dataFrame[[vn]] <- sapply(dataFrame[[vn]], function(x) max(x,m))
  }

  pred <- predict(gssObj,dataFrame)

  tolambda <- switch(gssObj$family,
                     binomial=function(x) exp(x)/(1+exp(x)),
                     poisson=function(x) exp(x)
                     )
  
  mkCPSP(cumsum(tolambda(pred))[dataFrame$event==1])
  
}

isi <- function (dataFrame,
                 lag = 1
                 ){

  if (!all(c("event","time") %in% names(dataFrame))) {
    stop("dataFrame should contain both an event and a time variable")
  }
  lag <- round(lag)
  if (lag < 1) stop("lag should be a strictly positive integer")
  
  event <- dataFrame$event
  time <- dataFrame$time
  if (!("trial" %in% names(dataFrame))) {
    n <- length(event)
    if (length(time) != n) 
      stop("event and time should have the same length.")
    result <- numeric(n)
    nb <- sum(event == 1)
    event.idx <- (1:n)[event == 1]
    isi <- diff(time[event == 1])
    result[1:event.idx[1 + lag]] <- NA
    for (i in (lag + 1):(nb - 1)) {
      result[(event.idx[i] + 1):event.idx[i + 1]] <- isi[i - lag]
    }
    if (n > event.idx[nb]) 
      result[(event.idx[nb] + 1):n] <- isi[nb - lag]
  } else {
    trial <- dataFrame$trial
    result <- lapply(levels(trial),
                     function(tIdx) {
                       evt <- event[trial == tIdx]
                       timetrial <- time[trial == tIdx]
                       n <- length(evt)
                       if (length(timetrial) != n) 
                         stop("eventt and time length mismatch among trials.")
                       partialResult <- numeric(n)
                       nb <- sum(evt == 1)
                       evt.idx <- (1:n)[evt == 1]
                       isi <- diff(timetrial[evt == 1])
                       partialResult[1:evt.idx[1 + lag]] <- NA
                       for (i in (lag + 1):(nb - 1)) {
                         partialResult[(evt.idx[i] + 1):evt.idx[i + 1]] <- isi[i - lag]
                       }
                       if (n > evt.idx[nb]) 
                         partialResult[(evt.idx[nb] + 1):n] <- isi[nb - lag]

                       partialResult
                     }
                     )
    result <- unlist(result)
  }
  result
}
