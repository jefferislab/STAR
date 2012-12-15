##############################################################################
##############################################################################
##############################################################################


#'Makes a Smooth Function Mapping a Data Frame Variable Onto a Variable Uniform
#'on Its Definition Domain
#'
#'The smooth transformation function is a smooth version of the
#'\code{\link{ecdf}}. A smooth density estimate as well as the inverse
#'transformation (the quantile function) are also returned as attributes.
#'
#'The smooth mapping to uniform function returned by \code{mkM2U} is obtained
#'by first selecting a subset of the variable values for which the variable
#'\code{time} of \code{df} is between \code{low} and \code{high}. The values
#'are then binned between the \code{min} and the \code{max} of the (complete)
#'variable values with a bin width \code{delta}. Function
#'\code{\link[gss]{ssden}} is then called on the histogram and the result is
#'stored in object \code{ii.fit} (This object is stored in the \code{closure}
#'of the returned function). The returned function is the result of a call of
#'\code{pssden} on ii.fit and the argument.
#'
#'A function inverting the "mapping to uniform function", that is, a quantile
#'function, is also returned as \link{attributes} \code{qFct}. This inverse
#'function is obtained by numerical inversion, calling \code{\link{uniroot}}
#'internally. Additional arguments can be passed to \code{\link{uniroot}} via
#'the \dots{} argument of the function.
#'
#'A function returning the smooth density estimate is returned as
#'\link{attributes} \code{dFct}.
#'
#'@param df a data frame. This data frame should contain a variable \code{time}
#'like data frames returned by \code{\link[STAR]{mkGLMdf}}.
#'@param vN a character string corresponding to the name of one of the
#'variables of \code{df} or an integer, its index. Variable \code{vN} is the
#'one for which the mapping to uniform is looked for.
#'@param low a numeric, the smallest value of variable \code{time} from which
#'the transformation is looked for. If missing defaults to the smallest time.
#'@param high a numeric, the largest value of variable \code{time} up to which
#'the transformation is looked for. If missing defaults to the largest time.
#'@param delta a numeric, the bin width used to build the variable values
#'histogram. This histogram is subsequently smoothed. Default provided if
#'missing.
#'@param alpha see \code{\link[gss]{ssden}}.
#'@param \dots additional arguments passed to \code{\link[gss]{ssden}} called
#'internally by the function
#'@return A function returning the probability for \code{vN} random variable to
#'have a value smaller or equal to its first argument. The returned function
#'calls internally \code{\link{integrate}}. Additional arguments can be passed
#'to the latter via the \dots{} argument of the returned function.
#'
#'As explained in the \code{details} section, the returned function has the
#'smooth density function, \code{dFct}, as well as the inverse function,
#'\code{qFct}, as attributes. Attribute \code{call} contains the matched call
#'and \code{range} contains the full range of the mapped variable.
#'@note Since the density returned by \code{\link[gss]{dssden}} can sometime
#'integrate to a value slightly different from 1 on its definition domain, the
#'actual integral is evaluated with \code{\link{integrate}} and the returned
#'density is renormalised. A look-up table of 101 regularly spaced quantiles
#'and the corresponding probabilities is also created and stored in the
#'returned function closure. This look-up table is used to speed up the
#'computations performed by the returned function which uses
#'\code{\link{integrate}} and not \code{pssden}. It is also used to speed up
#'the computations of the inverse function (returned as attribute \code{qFct})
#'which uses \code{\link{uniroot}} and not \code{qssden}.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link[gss]{ssden}}, \code{\link[gss]{dssden}},
#'\code{\link{integrate}}, \code{\link{uniroot}}, \code{\link[STAR]{mkGLMdf}}
#'@keywords models regression
#'@examples
#'
#'\dontrun{
#'require(STAR)
#'data(e060824spont)
#'DFA <- subset(mkGLMdf(e060824spont,0.004,0,59),neuron==1)
#'DFA <- within(DFA,i1 <- isi(DFA,lag=1))
#'DFA <- DFA[complete.cases(DFA),]
#'m2u1 <- mkM2U(DFA,"lN.1",0,29)
#'m2ui <- mkM2U(DFA,"i1",0,29,maxiter=200)
#'DFA <- within(DFA,e1t <- m2u1(lN.1))
#'DFA <- within(DFA,i1t <- m2ui(i1))
#'with(DFA,plot(ecdf(e1t[time>29]),pch="."))
#'abline(a=0,b=1,col=2,lty=2)
#'with(DFA,plot(ecdf(i1t[time>29]),pch="."))
#'abline(a=0,b=1,col=2,lty=2)}
#'
mkM2U <- function(df,
                  vN,
                  low,
                  high,
                  delta,
                  alpha=2,
                  ...) {
#######################################################################
### Function mkM2U
### Finds the transformation making the distribution of one variable of 
### a data frame uniform on its definition domain. The transformatin
### is a smooth version of the ecdf computed optionally on a subset of
### of the data frame.
### ----------------------------------------------------------
### Arguments:
###  df: a data frame.
###  vN: a character or an integer specifying the variable to transform.
###  low: a numeric specifying the time from which the ecdf should be
###       be evaluated. If missing it is set to the smallest time.
###  high: a numeric specifying the time up to which the ecdf should be
###        be evaluated. If missing it is set to the largest time.
###  delta: the bin width used to build the histogram of the variable
###         values which is later smoothed.
###  alpha: alpha argument of ssden see this gss function documentation
###  ...: additional arguments passed to ssden
### -------------------------------------------------------------------
### Value: a smooth version of the ecdf.
### The returned function has in addition 4 attributes:
### qFct: the smooth quantile function, that is if mySmooth is what is
###       returned by mkM2U, then attr(mySmooth,"qFct")(mySmooth(x))=x
### dFct: the smooth pdf of the selected variable
### range: a two elements vector with the range of the variable
### call: the matched call
#######################################################################
  
  ii <- df[[vN]]
  if (missing(low)) low <- min(df$time)
  if (missing(high)) high <- max(df$time)  
  if (missing(delta)) {
    iii <- sort(diff(sort(unique(ii))))
    delta <- min(iii[iii>=diff(range(iii))/1000])
    rm(iii)
  } ## End of conditional on missing(delta)
  
  iiD <- c(min(df[[vN]])-.Machine$double.eps,
           max(df[[vN]])+.Machine$double.eps)
  iiRef <- ii[low <= df$time & df$time <= high]
  rm(df,ii,low,high)

  iiB <- seq(iiD[1],iiD[2]+delta,delta)
  iiH <- as.data.frame(hist(iiRef,breaks=iiB,plot=FALSE)[c("mids","counts")])
  names(iiH) <- c("x","counts")
  riiD <- range(iiH$x)+delta*c(-1,1)/2
  ii.fit <- ssden(~x, data=iiH,
                  domain=data.frame(x=riiD),
                  weights=iiH$counts,alpha=alpha,
                  ...
                  )
  rm(iiH,iiB,iiRef,iiD)

  mn <- min(ii.fit$domain)
  mx <- max(ii.fit$domain)
  ## find out the true integral of the density function on its definition
  ## domain
  Z <- integrate(function(x) dssden(ii.fit,x),mn,mx,subdivisions=1000)$value
  ## define dFct returning a density actually summing to 1
  dFct <- function(x) {
    result <- numeric(length(x))
    good <- mn <= x & x <= mx
    if (any(good)) result[good] <- dssden(ii.fit,x[good])/Z
    result
  }

  ## create a lookup table of quantiles and corresponding probabilities
  Q <- seq(mn,mx,len=101)
  P <- numeric(101)
  for (i in 2:101) P[i] <- P[i-1] + integrate(dFct,Q[i-1],Q[i])$value
  
  result <- function(q,...) {
    order.q <- rank(q)
    p <- q <- sort(q)
    q.dup <- duplicated(q)
    p[q <= mn] <- 0
    p[q >= mx] <- 1
    kk <- (1:length(q))[q > mn & q < mx]
    p.dup <- 0
    for (i in kk) {
        if (q.dup[i]) {
            p[i] <- p.dup
            next
        }
        idx <- findInterval(q[i],Q)
        lwr <- Q[idx]
        if (q[i]-lwr < 1e-5) {
          p[i] <- P[idx]+(P[idx+1]-P[idx])/(Q[idx+1]-lwr)*(q[i]-lwr)
        } else {
          p[i] <- integrate(dFct,
                            lower = lwr, upper = q[i],
                            ...)$value + P[idx]
          p.dup <- p[i]
        } ## End of conditional on q[i]-lwr < 1e-5 
    } ## End of loop on i
    p[order.q]
  }

  qFct <- function(p,...) {
    q <- p <- sort(p)
    p.dup <- duplicated(p)
    q[p<=0] <- mn
    q[p>=1] <- mx
    kk <- (1:length(p))[p > 0 & p < 1]
    for (i in kk) {
      if (p.dup[i]) {
        q[i] <- q[i-1]
      } else {
        idx <- findInterval(p[i],P)
        if (identical(p[i],P[idx])) {
          q[i] <- Q[idx]
        } else {
          q[i] <- uniroot(function(x) result(x) - p[i],Q[idx:(idx+1)],...)$root
        }
      } ## End of conditional on p.dup[i]
    } ## End of loop on i
    q
  }
  
  attr(result,"qFct") <- qFct 
  attr(result,"dFct") <- dFct
  attr(result,"range") <- riiD
  attr(result,"call") <- match.call()
  result
}
##############################################################################
##############################################################################
##############################################################################


#'Generate a Data Frame With Variables Suitable For an AR Like Model
#'
#'The variables added to the data frame corresponding to the first argument of
#'the function are the former inter spike intervals. These variables are
#'moreover transformed with \code{mkM2U} so that they have an approximately
#'uniform distribution on their definition domain.
#'
#'When \code{max.order} > 1 the previous inter spike intervals are all
#'transformed using the "map to uniform" function estimated from the inter
#'spike intervals at lag 1.
#'
#'@param df a data frame. This data frame should contain a variable \code{time}
#'like data frames returned by \code{\link[STAR]{mkGLMdf}}.
#'@param low a numeric, the smallest value of variable \code{time} from which
#'the transformation is looked for. If missing defaults to the smallest time.
#'@param high a numeric, the largest value of variable \code{time} up to which
#'the transformation is looked for. If missing defaults to the largest time.
#'@param max.order a postive integer, the maximal order of the AR model. How
#'many previous inter spike intervals should be used in order to predict the
#'duration of the next interval?
#'@param selfName a character string or an integer specifying the variable of
#'\code{df} containing the elapsed time since the last spike of the considered
#'neuron.
#'@param \dots additional arguments passed to \code{\link[STAR]{mkM2U}}
#'@return A data frame is returned. In addition to the variables of df the
#'returned data frame contains a variable \code{est} with the transformed
#'elapsed time since the last spike of the neuron and \code{i1t},
#'\code{i2t},\dots{},\code{i max.order t}, the transformed previous inter spike
#'intervals.  The returned data frame has also four attributes:
#'@returnItem fmla a formula suitable for a first argument of, say,
#'\code{\link[gss]{gssanova}}.
#'@returnItem m2uL the function returned by mkM2U transforming the elasped time
#'since the last spike of the neuron.
#'@returnItem m2uI the function returned by mkM2U transforming the first former
#'inter spike interval.
#'@returnItem call the matched call.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link[STAR]{mkM2U}}, \code{\link[gss]{gssanova}}
#'@keywords models regression ts
#'@examples
#'
#'\dontrun{
#'require(STAR)
#'data(e060824spont)
#'DFA <- subset(mkGLMdf(e060824spont,0.004,0,59),neuron==1)
#'DFA <- mkAR(DFA, 0, 29, 5, maxiter=200)
#'head(DFA)
#'tail(DFA)
#'ar.fit <- gssanova(attr(DFA,"fmla"), data=DFA,family="binomial",seed=20061001)
#'plot(ar.fit %qp% "est")
#'plot(ar.fit %qp% "i1t")
#'plot(ar.fit %qp% "i2t")
#'plot(ar.fit %qp% "i3t")
#'plot(ar.fit %qp% "i4t")
#'plot(ar.fit %qp% "i5t")
#'}
#'
mkAR <- function(df,
                 low,
                 high,
                 max.order,
                 selfName="lN.1",
                 ...
                 ) {
#######################################################################
### Function mkAR
### Generates a data frame with variables suitable for an AR like model.
### These variables are the former inter spike intervals. The variables
### are moreover transformed with mkM2U inorder to have a uniform
### distribution on their definition domain.
### ----------------------------------------------------------
### Arguments:
###  df: a data frame.
###  low: a numeric specifying the time from which the ecdf should be
###       be evaluated. If missing it is set to the smallest time.
###  high: a numeric specifying the time up to which the ecdf should be
###        be evaluated. If missing it is set to the largest time.  
###  max.order: the maximal order of the process.
###  selfName: the name of the variable of df containing the elapsed time
###            since the last spike of the neuron whose discharge is
###            modeled.
###  ... : additional arguments passed to mkM2U
### -------------------------------------------------------------------
### Value: a data frame.
### In addition to the variables of df the returned data frame contains
### a variable "est" with the transformed elapsed time since the last
### spike of the neuron and "i1t", "i2t",...,"i max.order t", the
### transformed previous inter spike intervals.
### The returned data frame has also four attributes:
###  fmla: a formula suitable for a first argument of, say, gssanova
###  tfL: the function returned by mkM2U transforming the elasped time
###       since the last spike of the neuron.
###  tfI: the function returned by mkM2U transforming the first former
###       inter spike interval.
###  call: the matched call.
#######################################################################  
  if (missing(low)) low <- min(df$time) - .Machine$double.eps
  if (missing(high)) high <- max(df$time) + .Machine$double.eps
  vNames <- paste("i",1:max.order,sep="")
  vNames2 <- paste("i",1:max.order,"t",sep="")
  for (i in 1:max.order) df[[vNames[i]]] <- isi(df,lag=i)
  df <- df[complete.cases(df),]
  m2uL <- mkM2U(df,selfName,low,high,...)
  m2uI <- mkM2U(df,"i1",low,high,...)
  df <- within(df, est <- m2uL(lN.1))
  for (i in 1:max.order) df[[vNames2[i]]] <- m2uI(df[[vNames[i]]])
  for (i in 1:max.order) df[[vNames[i]]] <- NULL
  fmla <- as.formula(paste("event ~ est+", paste(vNames2, collapse= "+")))
  attr(df,"fmla") <- fmla
  attr(df,"m2uL") <- m2uL
  attr(df,"m2uI") <- m2uI
  attr(df,"call") <- match.call()
  df
}
##############################################################################
##############################################################################
##############################################################################


#'Change the Scales of a quickPredict Object for an Interaction Term
#'
#'Designed to transform results of \code{\link[STAR]{quickPredict}} obtained on
#'interaction terms from the transformed scale (on which the variables are
#'approximately uniformly distributed) onto the "native", linear scale.
#'
#'
#'@param obj a \code{\link[STAR]{quickPredict}} object.
#'@param xFct a function to be applied on the \code{xx} element of \code{obj}.
#'This function should be the \code{qFct} attribute of the function, returned
#'by \code{\link[STAR]{mkM2U}}, used to transform the variable from the
#'"native" to the "uniform" scale.
#'@param yFct a function to be applied on the \code{yy} element of \code{obj}.
#'This function should be the \code{qFct} attribute of the function, returned
#'by \code{\link[STAR]{mkM2U}}, used to transform the variable from the
#'"native" to the "uniform" scale.
#'@return A \code{\link[STAR]{quickPredict}} object.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link[STAR]{quickPredict}},
#'\code{\link[STAR]{plot.quickPredict}}
#'@keywords models regression
#'@examples
#'
#'\dontrun{
#'data(e060824spont)
#'DFA <- subset(mkGLMdf(e060824spont,0.004,0,59),neuron==1)
#'DFA <- within(DFA,i1 <- isi(DFA,lag=1))
#'DFA <- DFA[complete.cases(DFA),]
#'m2u1 <- mkM2U(DFA,"lN.1",0,29)
#'m2ui <- mkM2U(DFA,"i1",0,29,maxiter=200)
#'DFA <- within(DFA,e1t <- m2u1(lN.1))
#'DFA <- within(DFA,i1t <- m2ui(i1))
#'with(DFA,plot(ecdf(e1t[time>29]),pch="."))
#'abline(a=0,b=1,col=2,lty=2)
#'with(DFA,plot(ecdf(i1t[time>29]),pch="."))
#'abline(a=0,b=1,col=2,lty=2)
#'m1.fit <- gssanova(event~e1t*i1t, data=subset(DFA,time>29), family="binomial", seed=20061001)
#'inter.pred <- m1.fit %qp% "e1t:i1t"
#'contour(inter.pred,what="mean",nlevels=10,col=2,lwd=2)
#'contour(inter.pred,what="sd",nlevels=5,col=1,lwd=1,lty=2,add=TRUE)
#'inter.predN <- changeScale(inter.pred,attr(m2u1,"qFct"),attr(m2ui,"qFct"))
#'contour(inter.predN,what="mean",nlevels=5,col=2,lwd=1)
#'contour(inter.predN,what="sd",nlevels=3,col=1,lwd=1,lty=2,add=TRUE)
#'}
#'
changeScale <- function(obj,
                        xFct,
                        yFct
                        ) {
#######################################################################
### Function changeScale
### Designed to transform results of quickPredict obtained on interaction
### terms from the transformed scale (on which the variables are approxi-
### mately uniformly distributed) onto the "native", linear scale.
### ----------------------------------------------------------
### Arguments:
###  obj: a quickPredict object.
###  xFct: a function to be applied on the "xx" element of obj.
###  yFct: a function to be applied on the "yy" element of obj.
### -------------------------------------------------------------------
### Value: a quickPredict object.
#######################################################################  
  ## Check that obj inherits of "quickPredict"
  if (!inherits(obj,"quickPredict")) stop("obj should be a quickPredict object.")

  if (missing(xFct)) xFct <- function(x) x
  if (missing(yFct)) yFct <- function(y) y

  est.mean <- obj$est.mean
  if (!is.null(obj$est.sd)) {
    est.sd <- obj$est.sd
  } else {
    est.sd <- NULL
  }

  theCall <- match.call()
  theCall[["obj"]] <- obj$call
  xx <- xFct(obj$xx)
  n <- length(xx)
  equal2min <- sum(abs(xx-min(xx)) <= .Machine$double.eps)
  equal2max <- sum(abs(xx-max(xx)) <= .Machine$double.eps)
  goodX <- !logical(n)
  if (equal2min > 1) goodX[1:(equal2min-1)] <- FALSE
  if (equal2max > 1) goodX[(n-equal2max+2):n] <- FALSE
  est.mean <- est.mean[goodX,]
  if (!is.null(est.sd)) est.sd <- est.sd[goodX,]
  xx <- xx[goodX]
  goodX <- !duplicated(xx)
  xx <- xx[goodX]
  est.mean <- est.mean[goodX,]
  if (!is.null(est.sd)) est.sd <- est.sd[goodX,]
  
  yy <- yFct(obj$yy)
  equal2min <- sum(abs(yy-min(yy)) <= .Machine$double.eps)
  equal2max <- sum(abs(yy-max(yy)) <= .Machine$double.eps)
  goodY <- !logical(n)
  if (equal2min > 1) goodY[1:(equal2min-1)] <- FALSE
  if (equal2max > 1) goodY[(n-equal2max+2):n] <- FALSE
  est.mean <- est.mean[,goodY]
  if (!is.null(est.sd)) est.sd <- est.sd[,goodY]
  yy <- yy[goodY]
  goodY <- !duplicated(yy)
  yy <- yy[goodY]
  est.mean <- est.mean[,goodY]
  if (!is.null(est.sd)) est.sd <- est.sd[,goodY]

  result <- list(xx=xx,
                 yy=yy,
                 call=theCall,
                 include=obj$include,
                 est.mean=est.mean,
                 est.sd=est.sd)
  class(result) <- "quickPredict"
  result
  
}
