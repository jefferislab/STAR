##############################################################################
##############################################################################
##############################################################################


#'A Simple Interface to predict method for ssanova and ssanova0 objects
#'
#'Designed to quickly compute the effect of a \emph{single} model term. This
#'term can correspond to a single variable effect or to the interaction of two
#'variables.
#'
#'\code{\%qp\%} is the binary version of \code{quickPredict}.
#'
#'@aliases quickPredict %qp%
#'@param object an object inheriting from ssanova and ssanova0
#'(\code{\link[gss]{gssanova}} and \code{\link[gss]{gssanova0}} objects are
#'therefore suitable).
#'@param include a character string corresponding to a \emph{single} model
#'term. See \code{\link[gss]{predict.ssanova}} and
#'\code{\link[gss]{predict.ssanova}}.
#'@param se.fit logical flag indicating if standard errors are required. See
#'\code{\link[gss]{predict.ssanova}} and \code{\link[gss]{predict.ssanova}}.
#'@param length.out a positive integer, the number of points at which the
#'prediction should be performed. These points are uniformly spread on the
#'definition domain of the variable(s) implicitely specified by argument
#'\code{include}. If \code{missing} a default of 501 for terms involving a
#'single variable and of 101 for interaction terms involving two variables is
#'provided.
#'@param otherTermsFct a function applied to the other variables required for
#'model specification.
#'@return A \code{quickPredict} object. This object is a \code{\link{list}}
#'with the following components:
#'@returnItem xx a numeric vector with the values of the variable specified by
#'the model term selected by argument \code{include}. When an interaction term
#'was selected the values of the first variable are stored here.
#'@returnItem yy a numeric vector with the values of the \emph{second} variable
#'specified by the \emph{interaction} term selected by argument \code{include}.
#'When selected term is not an interaction term, this component is \code{NULL}.
#'@returnItem include the value of the argument with this name.
#'@returnItem call the matched call.
#'@returnItem est.mean a numeric vector or matrix, for intercation terms,
#'containing the estimated mean of the term.
#'@returnItem est.sd a numeric vector or matrix, for intercation terms,
#'containing the estimated SD of the term. Is \code{NULL} is argument
#'\code{se.fit} was \code{FALSE}.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link[gss]{predict.ssanova}},
#'\code{\link[gss]{predict.ssanova}}, \code{\link{plot.quickPredict}},
#'\code{\link{image.quickPredict}}, \code{\link{contour.quickPredict}},
#'\code{\link{persp.quickPredict}}, \code{\link{plot.ssanova}}
#'@keywords models smooth regression
#'@examples
#'
#'## Follow up of ssanova example of gss
#'data(nox)
#'nox.fit <- ssanova(log10(nox)~comp*equi,data=nox)
#'## get prediction for the first term, comp
#'comp.pred <- quickPredict(nox.fit)
#'## plot result with method plot for quickPredict objects
#'plot(comp.pred)
#'## get prediction for the second term, equi using the binary version
#'equi.pred <- nox.fit \%qp\% "equi"
#'plot(equi.pred)
#'## get prediction for the interaction term, comp:equi
#'comp.equi.pred <- nox.fit %qp% "comp:equi"
#'## use image method image
#'image(comp.equi.pred)
#'## use contour method
#'contour(comp.equi.pred,col=2,lwd=2,labcex=1.5)
#'contour(comp.equi.pred,what="sd",lty=3,labcex=1.2,add=TRUE)
#'## use persp method
#'persp(comp.equi.pred,theta=-10,phi=20)
#'
quickPredict <- function(object,
                         include=object$terms$labels[2],
                         se.fit=TRUE,
                         length.out,
                         otherTermsFct=median
                         )
{
######################################################################
### Function quickPredict: an "easy" interface to predict.ssanova and
### predict.ssanova0
### Designed to quickly compute the effect of a _single_ model term.
### This term can correspond to a single variable effect or to the
### interaction of two variables
### -----------------------------------------------------------------
### Arguments:
###  object: an object inheriting from ssanova and ssanova0 (gssanova
###          and gssanova0 objects are therefore suitable).
###  include: a character, the model term for which one wants the
###           the prediction.
###  se.fit: see predict.ssanova and predict.ssanova0
###  length.out: a positive integer, the number of points on the
###              variable / term definition domain at which the
###              "prediction" will be made.
###  otherTermsFct: a function used to set the value of the other
###                 model terms in the data frame fed to
###                 predict.ssanova / predict.ssanova0
### -------------------------------------------------------------------
### Value:
###  A "quickPredict" object. This object is a list with the following
###  components:
###  xx: a numeric vector with the values of the selected term at
###      which prediction was made. When an interaction term was
###      selected the values of the first variable are stored here.
###  yy: values of the second variable (for interaction terms) at
###      which prediction was made. NULL for none interaction terms.
###  include: the value of the argument with this name.
###  call: the matched call.
###  est.mean: a numeric vector, the estimated mean value of the term,
###            or a matrix for interaction terms.
###  est.sd: a numeric vector or NULL, the estimated standard
###          deviation of the term. It is NULL if argument "se.fit"
###          is set to FALSE. It is a matrix for interaction terms.
##########################################################################
  
  ## Make sure that object inherits from ssanova or ssanova0
  if (!inherits(object,c("ssanova","ssanova0")))
    stop("object should be a ssanova or a ssanova0 obbject.")

  allTerms <- object$terms$labels[-1]
  if (length(include) != 1) stop("include should be a character of length 1.")
  if (!(include %in% allTerms)) stop("include should be one of the model terms.")

  ## Find out if the term is an interaction term
  isInteraction <- !(include %in% names(object$mf))

  if (missing(length.out)) {
    if (isInteraction) length.out <- 101
    else length.out <- 501
  } ## End of conditional on missing(length.out)
  
  if (!isInteraction) {
    xx <- seq(from=object$terms[[include]]$rk$env$env$min,
              to=object$terms[[include]]$rk$env$env$max,
              length.out=length.out)
    yy <- NULL
    newdata <- data.frame(xx)
    names(newdata) <- include
    otherTerms <- allTerms[!(allTerms %in% include) &
                           (allTerms %in% names(object$mf))
                           ]
  } else {
    vNames <- strsplit(include,":")[[1]]
    xx <- seq(from=object$terms[[include]]$rk$env$rk[[1]]$env$min,
              to=object$terms[[include]]$rk$env$rk[[1]]$env$max,
              length.out=length.out)
    yy <- seq(from=object$terms[[include]]$rk$env$rk[[2]]$env$min,
              to=object$terms[[include]]$rk$env$rk[[2]]$env$max,
              length.out=length.out)
    newdata <- data.frame(rep(xx,length.out),
                          rep(yy,rep(length.out,length.out))
                          )
    names(newdata) <- vNames
    otherTerms <- allTerms[!(allTerms %in% vNames) &
                           (allTerms %in% names(object$mf))
                           ]
  } ## End of conditional on !isInteraction

  nbOtherTerms <- length(otherTerms)
  if (nbOtherTerms > 0) {
    otherVal <- sapply(otherTerms, function(n) otherTermsFct(object$mf[,n]))
    for (idx in 1:nbOtherTerms) {
      newdata[[otherTerms[idx]]] <- rep(otherVal[idx],dim(newdata)[1])
    } ## End of for loop on idx
  } ## End of conditional on nbOtherTerms > 0 

  est <- predict(object,newdata,se.fit=se.fit,include=include)

  result <- list(xx=xx,
                 yy=yy,
                 call=match.call(),
                 include=include)
  if (se.fit) {
    if (isInteraction) {
      result$est.mean <- matrix(est$fit,length.out,length.out)
      result$est.sd <- matrix(est$se.fit,length.out,length.out)
    } else {
      result$est.mean <- est$fit
      result$est.sd <- est$se.fit
    }                  
  } else {
    if (isInteraction) {
      result$est.mean <- matrix(est,length.out,length.out)
    } else {
      result$est.mean <- est
    }
    result$est.sd <- NULL
  } ## End of conditional on se.fit

  class(result) <- "quickPredict"
  result
}
##############################################################################
##############################################################################
##############################################################################
"%qp%" <- function(object,
                   include
                   )
{
#######################################################################
### Function %qp%
### Binary version of quickPredict
### ----------------------------------------------------------
### Arguments:
###  object: an object inheriting from ssanova and ssanova0 (gssanova
###          and gssanova0 objects are therefore suitable).
###  include: a character, the model term for which one wants the
###           the prediction.
### The other arguments of quickPredict are set to their default values.
### -------------------------------------------------------------------
### Value: see quickPredict.
#######################################################################
  x <- deparse(substitute(object))
  y <- deparse(substitute(include))
  cmd <- paste("quickPredict(",
               x,
               ",",
               y,
               ")",sep="")
  cmd <- parse(text=cmd)
  eval(cmd)

}
