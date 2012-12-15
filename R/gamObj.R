#'Generic Function and Methods for Extracting a gamObject
#'
#'Some functions of \code{STAR}, like \code{gampsth} and \code{gamlockedTrain}
#'perform gam fits internally and keep as a list component or within the
#'environment of a returned function the result of this fit. Method
#'\code{gamObj} extracts this gam object.
#'
#'
#'@aliases gamObj gamObj.gampsth gamObj.gamlockedTrain
#'@param object an object containing a \code{gamObject}. Currently the result
#'of a call to \code{\link{gampsth}} or to \code{\link{gamlockedTrain}}.
#'@param \dots not used for now
#'@return A \code{\link[mgcv]{gamObject}}
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link[mgcv]{gam}}, \code{\link[mgcv]{gamObject}},
#'\code{\link{gampsth}}, \code{\link{gamlockedTrain}}
#'@keywords models smooth regression
#'@examples
#'
#'##
#'
gamObj <- function(object,
                   ...) {

  UseMethod("gamObj")
}

gamObj.gampsth <- function(object,...) {

  evalq(PoissonF, envir=environment(object$lambdaFct))

}

gamObj.gamlockedTrain <- function(object, ...) {
  object$gamFit
}
