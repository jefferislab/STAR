#'Generates a Data Frame of Dummy Variables for Use in gam
#'
#'Using argument \code{by} in \code{\link[mgcv]{s}} or \code{\link[mgcv]{te}}
#'of \code{\link[mgcv]{gam}} requires dummy variables to be set up. This is the
#'job of this function.
#'
#'
#'@param x a \code{factor}.
#'@return A \code{\link{data.frame}} with as many variables as there are
#'\code{levels} in \code{x} and as many rows as elements in \code{x}.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{mkGLMdf}}, \code{\link[mgcv]{gam}},
#'\code{\link[mgcv]{s}}, \code{\link[mgcv]{te}}
#'@keywords models
#'@examples
#'
#'## coming soon
#'
mkDummy <- function(x) {

  xN <- deparse(substitute(x))
  if (!inherits(x,"factor")) stop("x should be a factor")
  lev <- levels(x)

  result <- lapply(lev,function(l) as.integer(x==l))
  names(result) <- paste(xN,".",lev,sep="")
  result

}
