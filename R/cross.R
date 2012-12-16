crossTight <- function(tMax=1,
                       h=0.001,
                       a=0.3,
                       b=2.35,
                       withBounds=TRUE,
                       logScale=FALSE) {

  n <- round(tMax/h)
  tMax <- n*h

  if (!withBounds && !logScale) {
    result <- .C(crossTight_sym,
                 as.double(tMax),
                 as.integer(n),
                 as.double(a),
                 as.double(b),
                 G=double(n)
                 )$G
    result <- list(G=result)
  }
  if (!withBounds && logScale) {
    result <- .C(crossTightlog,
                 as.double(tMax),
                 as.integer(n),
                 as.double(a),
                 as.double(b),
                 G=double(n)
                 )$G
    result <- list(G=result)
  }
  if (withBounds && !logScale) {
    result <- .C(crossTightWithB,
                 as.double(tMax),
                 as.integer(n),
                 as.double(a),
                 as.double(b),
                 G=double(n),
                 Gu=double(n),
                 Gl=double(n)
                 )[c("Gu","G","Gl")]
  }
  if (withBounds && logScale) {
    result <- .C(crossTightWithBlog,
                 as.double(tMax),
                 as.integer(n),
                 as.double(a),
                 as.double(b),
                 G=double(n),
                 Gu=double(n),
                 Gl=double(n)
                 )[c("Gu","G","Gl")]
  }

  result$time <- (1:n)*h
  result$g <- diff(c(0,result$G))/h
  result$mids <- (1:n)*h-0.5*h
  result$h <- h
  result$call <- match.call()
  class(result) <- "FirstPassageTime"
  result
  
}

mkTightBMtargetFct <- function(ci=0.95,
                               tMax=1,
                               h=0.001,
                               logScale=FALSE) {

  ## Check ci
  if (!(0<ci && ci <1)) stop("Expect 0 < ci < 1")
  trueTarget <- (1-ci)/2
  n <- round(tMax/h)
  function(logp) {
    p <- exp(logp)
    (trueTarget - crossTight(tMax,h,p[1],p[2],FALSE,logScale)$G[n])^2
  }
  
}



#'Computations of Boundary Crossing Probabilities for the Wiener Process
#'
#'Computes the distribution of the first passage time through an arbitrary
#'(\code{crossGeneral}) or a "tight" (\code{crossTight}) boundary for a Wiener
#'process. The method of Loader and Deely (1987) is used. A tight boundary is a
#'boundary generating the tighest confidence band for the process (Kendall et
#'al, 2007). Utility function and methods: \code{mkTightBMtargetFct},
#'\code{print}, \code{summary}, \code{plot}, \code{lines}, are also provided to
#'use and explore the results.
#'
#'The Loader and Deely (1987) method to compute the probability
#'\eqn{G(t)}{G(t)} that the first passage of a Wiener process / Brownian motion
#'occurs between 0 and \eqn{t}{t} (argument \code{tMax} of \code{crossGeneral}
#'and \code{crossTight}) through a boundary defined by \eqn{c(t)}{c(t)} is
#'based on the numerical solution of a Volterra integral equation of the first
#'kind satisfied by \eqn{G()}{G()} and defined by their Eq. 2.2: \deqn{F(t) =
#'\int_0^t K(t,u) dG(u)}{F(t) = \int_0^t K(t,u) dG(u)} where, \eqn{F(t)}{F(t)}
#'is defined by: \deqn{F(t)=\Phi(-\frac{c(t)}{\sqrt{t}})+\exp \big( -2 b(t) \,
#'}{F(t) = pnorm(-c(t)/sqrt(t)) +
#'exp(-2*b(t)*(c(t)-t*b(t)))*pnorm((-c(t)+2*t*b(t))/sqrt(t))}\deqn{
#'(c(t)-tb(t))\big) \, \Phi(\frac{-c(t)+2\, t \, b(t)}{\sqrt{t}})}{F(t) =
#'pnorm(-c(t)/sqrt(t)) +
#'exp(-2*b(t)*(c(t)-t*b(t)))*pnorm((-c(t)+2*t*b(t))/sqrt(t))}
#'\eqn{K(t,u)}{K(t,u)} is defined by:
#'\deqn{K(t,u)=\Phi(\frac{c(u)-c(t)}{\sqrt{t-u}})+\exp \big(-2 b(t) \, (c(t)
#'}{K(t,u)=pnorm((c(u)-c(t))/sqrt(t-u)) +
#'exp(-2*b(t)*(c(t)-c(u)-(t-u)*b(t)))*pnorm((c(u)-c(t)+2*(t-u)*b(t))/(sqrt(t-u)))
#'}\deqn{ -c(u) -(t-u) b(t))\big) \, \Phi(\frac{c(u)-c(t)+2\, (t-u) \,
#'}{K(t,u)=pnorm((c(u)-c(t))/sqrt(t-u)) +
#'exp(-2*b(t)*(c(t)-c(u)-(t-u)*b(t)))*pnorm((c(u)-c(t)+2*(t-u)*b(t))/(sqrt(t-u)))
#'}\deqn{ b(t)}{\sqrt{t-u}})}{K(t,u)=pnorm((c(u)-c(t))/sqrt(t-u)) +
#'exp(-2*b(t)*(c(t)-c(u)-(t-u)*b(t)))*pnorm((c(u)-c(t)+2*(t-u)*b(t))/(sqrt(t-u)))
#'} and \eqn{b(t)}{b(t)} is an additional function (that can be uniformly 0)
#'that is chosen to improve convergence speed and to get error bounds. Argument
#'\code{h} is the step size used in the numerical solution of the above
#'Volterra integral equation. The mid-point method (Eq. 3.1 and 3.2 of Loader
#'and Deely (1987)) is implemented. If \code{tMax} is not a multiple of
#'\code{h} it is modified as follows: \code{tMax <- round(tMax/h)*h}.
#'
#'\code{crossGeneral} generates functions \eqn{F()}{F()} and \eqn{K(,)}{K(,)}
#'internally given \eqn{c()}{c()} (argument \code{cFct}) and \eqn{b()}{b()}
#'(argument \code{bFct}). If \code{bFct} is not given (\emph{i.e.},
#'\code{missing(bFct)} returns \code{TRUE}) it is taken as uniformly 0. If a
#'numeric is given for \code{cFct} then \eqn{c()}{c()} is defined as a uniform
#'function returning the first element of the argument (\code{cFct}).
#'
#'Function \code{crossTight} assumes the following functional form for
#'\eqn{c()}{c()}: \eqn{c(t) = a + b \, \sqrt(t)}{a + b * sqrt(t)}.
#'\eqn{b()}{b()} is set to \eqn{c'()}{c'()} (the derivative of \eqn{c()}{c()}).
#'Arguments \code{a} and \code{b} of \code{crossTight} correspond to the 2
#'parameters of \eqn{c()}{c()}.
#'
#'If argument \code{withBounds} is set to \code{TRUE} then bounds on
#'\eqn{G()}{G()} are computed. Function \code{crossTight} uses Eq. 3.6 and 3.7
#'of Loader and Deely (1987) to compute these bounds, \eqn{G_{L}(t)}{Gl(t)} and
#'\eqn{G_U(t)}{Gu(t)}. Function \code{crossGeneral} uses Eq. 3.6 and 3.7 (if
#'argument \code{Lplus} is set to \code{TRUE}) or Eq. 3.10 and 3.11 (if
#'argument \code{Lplus} is set to \code{FALSE}). Here \code{Lplus} stands for
#'the sign of the partial derivative of the kernel \eqn{K(,)}{K(,)} with
#'respect to its second argument. If the sign is not known the user can provide
#'the derivative \eqn{c'()}{c'()} of \eqn{c()}{c()} as argument
#'\code{cprimeFct}. A (slow) numerical check is then performed to decide wether
#'\code{Lplus} should be \code{TRUE} or \code{FALSE} or if it changes sign (in
#'which case bounds cannot be obtained and an error is returned).
#'
#'In function \code{crossTight} argument \code{logScale} controls the way some
#'intermediate computations of the mid-point method are implemented. If set to
#'\code{FALSE} (the default) a literal implementation of Eq. 3.2 of Loader and
#'Deely (1987) is used. If set to \code{TRUE} then additions subtractions are
#'computed on the log scale using functions \code{logspace_add} and
#'\code{logspace_sub} of the \code{R API}. The computation is then slightly
#'slower and it turns out that the gain in numerical precision is not really
#'significant, so you can safely leave this argument to its default value.
#'
#'@aliases crossGeneral crossTight mkTightBMtargetFct print.FirstPassageTime
#'summary.FirstPassageTime plot.FirstPassageTime lines.FirstPassageTime
#'@param tMax A positive \code{numeric}. The "time" during which the Wiener
#'process is followed.
#'@param h A positive \code{numeric}. The integration time step used for the
#'numerical solution of the Volterra integral equation (see \code{details}).
#'@param cFct A \code{function} defining the boundary to be crossed. The first
#'argument of the function should be a "time" argument. If the first argument
#'is a vector, the function should return a vector of the same length.
#'@param cprimeFct A \code{function} defining time derivative of the boundary
#'to be crossed. Needs to be specified only if a check of the sign of the
#'kernel derivative (see \code{details}) is requested. The first argument of
#'the function should be a "time" argument. If the first argument is a vector,
#'the function should return a vector of the same length.
#'@param bFct A \code{function}. The "b" function of Loader and Deely (1987).
#'Does not need to be specified (\emph{i.e.}, can be missing) but can be used
#'to improve convergence. The first argument of the function should be a "time"
#'argument. If the first argument is a vector, the function should return a
#'vector of the same length.
#'@param withBounds A \code{logical}. Should bounds on the distribution be
#'calculated? If yes, set it to \code{TRUE}, leave it to its default value,
#'\code{FALSE}, otherwise.
#'@param Lplus A \code{logical}. If bounds are requested
#'(\code{withBounds=TRUE}) and if the sign of the time derivative of the kernel
#'is known to be positive or null, set to \code{TRUE}, if it is known to be
#'negative, set it to \code{FALSE}. If the sign is unknown, leave \code{Lplus}
#'unspecified and provide a \code{cprimeFct} function.
#'@param logScale A \code{logical}. Should intermediate calculations in
#'\code{crossTight} be carried out on the log scale for numerical precision? If
#'yes, set it to \code{TRUE}, leave it to its default, \code{FALSE}, otherwise.
#'@param a,b \code{numerics}, the two parameters of the "tight" boundary:
#'\code{c(t) = a + b*sqrt(t)}. See details.
#'@param ci A \code{numeric} larger than 0 and smaller than 1. The nominal
#'coverage probability desired for a "tight" confidence band (see
#'\code{details}).
#'@param x,object A \code{FirstPassageTime} object returned by
#'\code{crossGeneral} or \code{crossTight}.
#'@param y Not used but required for a \code{plot} method.
#'@param which A \code{character} string, "\code{Distribution}" or
#'"\code{density}", specifying if a probability distribution or a probability
#'density should be graphed.
#'@param xlab,ylab See \code{\link{plot}}.
#'@param digits A positive \code{integer}. The number of digits to print in
#'\code{summary}. If bounds were computed, the value of \code{digits} is
#'computed internally based on the bounds width.
#'@param \dots Used in \code{plot} and \code{lines} to pass further arguments
#'(see \code{\link{plot}} and \code{\link{lines}}), not used in \code{print}
#'and \code{summary}.
#'@return \code{crossGeneral} and \code{crossTight} return a
#'\code{FirstPassageTime} object which is a \code{list} with the following
#'components:
#'
#'\code{mkTightBMtargetFct} returns a \code{function} which can be used in
#'optim. This function returns the square of the difference between
#'\code{(1-ci)/2} (remember the "symmetry" of the Wiener processes paths, that
#'is, for every path there is a symmetric one with respect to the abscissa
#'having \emph{with the same probability}) and the probability to have the
#'first passage time of the Wiener process smaller or equal to 1 when the
#'boundary is the "tight" boundary defined by: \eqn{a + b \, \sqrt{t}}{a +
#'b*sqrt(t)}. The function takes a single vector argument containing the
#'\emph{log} of the parameters \code{a} (vector's first element) and \code{b}
#'(vector's second element).
#'
#'Methods \code{print.FirstPassageTime} and \code{summary.FirstPassageTime}
#'output the probability to observe the first exit between 0 and \code{tMax}.
#'If bounds were computed, the precision on the probability is also returned
#'(as an attribute for \code{print.FirstPassageTime}).
#'\code{summary.FirstPassageTime} also gives the integration time step,
#'\code{h}, used.
#'@returnItem time A \code{numeric} vector of "times" at which the first
#'passage time probability has been evaluated.
#'@returnItem G A \code{numeric} vector of first passage probability.
#'@returnItem Gu A \code{numeric} vector with the upper bound of first passage
#'probability. Only if \code{withBounds} was set to \code{TRUE}.
#'@returnItem Gl A \code{numeric} vector with the lower bound of first passage
#'probability. Only if \code{withBounds} was set to \code{TRUE}.
#'@returnItem mids A \code{numeric} vector of "times" at which the first
#'passage time probability \emph{density} has been evaluated. Mid points of
#'component \code{time}.
#'@returnItem g A \code{numeric} vector of first passage probability
#'\emph{density}.
#'@returnItem h A \code{numeric}. The value of argument \code{h} of
#'\code{crossGeneral} or \code{crossTight}.
#'@returnItem call The matched call.
#'@note Using \code{logScale = TRUE} in \code{crossTight} seems to be an
#'overkill, that is, it doubles computation's time without bringing significant
#'numerical improvement.
#'
#'\code{crossGeneral} is for now pure \code{R} code. The first passage
#'probability is obtained by solving the lower triangular system (Eq. 3.1 of
#'Loader and Deely (1987)) with \code{\link{forwardsolve}} and is therefore
#'rather fast (but can be memory "hungry"). The bounds are computed by 2 nested
#'loops and can therefore be long to compute.
#'
#'\code{crossTight} is calling a \code{C} code and is fast.
#'
#'Loader and Deely paper also describes a method where \eqn{G(t)}{G(t)} is
#'solution of a Volterra integral equation of the second kind (their Eq. 2.7).
#'This approach is presently not implemented in the above functions.
#'@section Warning : \code{crossGeneral} with \code{withBounds = TRUE} and a
#'negative kernel derivative is presently poorly tested, so be careful and let
#'me know if mistakes show up.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{print}}, \code{\link{summary}}, \code{\link{plot}},
#'\code{\link{lines}}, \code{\link{pinvgauss}}
#'@references C. R. Loader and J. J. Deely (1987) Computations of Boundary
#'Crossing Probabilities for the Wiener Process. \emph{J. Statist. Comput.
#'Simul.} \bold{27}: 95--105.
#'
#'W. S. Kendall, J.- M. Marin and C. P. Robert (2007) Brownian Confidence Bands
#'on Monte Carlo Output. \emph{Statistics and Computing} \bold{17}: 1--10.
#'Preprint available at:
#'\url{http://www.ceremade.dauphine.fr/\%7Exian/kmr04.rev.pdf}
#'@keywords distribution htest
#'@examples
#'
#'\dontrun{
#'## Reproduce Table 1 (p 101) of Loader and Deely (1987)
#'## define a vector of n values
#'nLD <- c(8,16,32,64,128)
#'
#'## Part 1: c(t) = sqrt(1+t) and tMax=1
#'## define cFct
#'cFT1p1 <- function(t) sqrt(1+t)
#'## define the different bFct
#'bFT1p1.ii <- function(t) 0.5/sqrt(1+t)
#'bFT1p1.iii <- function(t) (cFT1p1(t)-cFT1p1(0))/t 
#'bFT1p1.iv <- function(t) 0.5*(bFT1p1.ii(t)+bFT1p1.iii(t)) 
#'bFT1p1.v <- function(t) (2*t-4/5*((1+t)^2.5-1))/t^3+3*cFT1p1(t)/2/t
#'## Do the calculations
#'round(t(sapply(nLD,
#'               function(n) {
#'                 c(n=n,
#'                   i=crossGeneral(tMax=1,h=1/n,cFct=cFT1p1)$G[n],
#'                   ii=crossGeneral(tMax=1,h=1/n,cFct=cFT1p1,bFct=bFT1p1.ii)$G[n],
#'                   iii=crossGeneral(tMax=1,h=1/n,cFct=cFT1p1,bFct=bFT1p1.iii)$G[n],
#'                   iv=crossGeneral(tMax=1,h=1/n,cFct=cFT1p1,bFct=bFT1p1.iv)$G[n],
#'                   v=crossGeneral(tMax=1,h=1/n,cFct=cFT1p1,bFct=bFT1p1.v)$G[n])})),
#'      digits=6)
#'
#'## Part 2: c(t) = exp(-t) and tMax=1
#'## define cFct
#'cFT1p2 <- function(t) exp(-t)
#'## define the different bFct
#'cFT1p2 <- function(t) exp(-t)
#'bFT1p2.ii <- function(t) -exp(-t)
#'bFT1p2.iii <- function(t) (cFT1p2(t)-cFT1p2(0))/t 
#'bFT1p2.iv <- function(t) 0.5*(bFT1p2.ii(t)+bFT1p2.iii(t)) 
#'bFT1p2.v <- function(t) 3*(1-t-exp(-t))/t^3+3*cFT1p2(t)/2/t
#'## Do the calculations
#'round(t(sapply(nLD,
#'               function(n) {
#'                 c(n=n,
#'                   i=crossGeneral(tMax=1,h=1/n,cFct=cFT1p2)$G[n],
#'                   ii=crossGeneral(tMax=1,h=1/n,cFct=cFT1p2,bFct=bFT1p2.ii)$G[n],
#'                   iii=crossGeneral(tMax=1,h=1/n,cFct=cFT1p2,bFct=bFT1p2.iii)$G[n],
#'                   iv=crossGeneral(tMax=1,h=1/n,cFct=cFT1p2,bFct=bFT1p2.iv)$G[n],
#'                   v=crossGeneral(tMax=1,h=1/n,cFct=cFT1p2,bFct=bFT1p2.v)$G[n])})),
#'      digits=6)
#'
#'## Part 3: c(t) = t^2 + 3*t + 1 and tMax=1
#'## define cFct
#'cFT1p3 <- function(t) t^2+3*t+1
#'## define the different bFct
#'bFT1p3.ii <- function(t) 2*t+3
#'bFT1p3.iii <- function(t) (cFT1p3(t)-cFT1p3(0))/t 
#'bFT1p3.v <- function(t) 5*t/4+3
#'bFT1p3.vi <- function(t) rep(3,length(t))
#'round(t(sapply(nLD,
#'               function(n) {
#'                 c(n=n,
#'                   i=crossGeneral(tMax=1,h=1/n,cFct=cFT1p3)$G[n],
#'                   ii=crossGeneral(tMax=1,h=1/n,cFct=cFT1p3,bFct=bFT1p3.ii)$G[n],
#'                   iii=crossGeneral(tMax=1,h=1/n,cFct=cFT1p3,bFct=bFT1p3.iii)$G[n],
#'                   v=crossGeneral(tMax=1,h=1/n,cFct=cFT1p3,bFct=bFT1p3.v)$G[n],
#'                   vi=crossGeneral(tMax=1,h=1/n,cFct=cFT1p3,bFct=bFT1p3.vi)$G[n])})),
#'      digits=6)
#'
#'## Part 3: c(t) = t^2 + 3*t + 1 and tMax=1
#'## define cFct
#'cFT1p4 <- function(t) 1/t
#'## Here only column (i) and (vii) are reproduced.
#'## If one attempts to reproduce (ii) directly with crossGeneral
#'## a NaN appears (when a -Inf is the correct value) in functions
#'## F() and K(,) (see details) for instance when t=0 in F.
#'## Then as crossGeneral is presently written R attempts to
#'## compute t*b(t), where b(t) is c'(t), that is, t*(-1/t^2) which is
#'## NaN (for R) when t=0.
#'bFT1p4.vii <- function(t) rep(-1,length(t))
#'round(t(sapply(nLD,
#'               function(n) {
#'                 c(n=n,
#'                   i=crossGeneral(tMax=1,h=1/n,cFct=cFT1p4)$G[n],
#'                   vii=crossGeneral(tMax=1,h=1/n,cFct=cFT1p4,bFct=bFT1p4.vii)$G[n])})),
#'      digits=6)
#'## The last 3 rows of column (vii) are not the same as in the paper
#'
#'## Reproduce Table 4 (p 104) of Loader and Deely (1987)
#'## As before the probability of first passage between
#'## 0 and 1 is computed. This time only three boundary
#'## functions are considered. The error bounds are
#'## obtained
#'
#'## Part 1: c(t) = sqrt(1+t)
#'## Left columns pair: b(t) = c'(t)
#'round(t(sapply(nLD,
#'               function(n) {
#'                 res <- crossGeneral(tMax=1,h=1/n,cFct=cFT1p1,bFct=bFT1p1.ii,withBounds=TRUE,Lplus=TRUE)
#'                 c(Gl=res$Gl[n],Gu=res$Gu[n])
#'               }
#'               )
#'         ),
#'       digits=5)
#'
#'## Right columns pair: b(t) = 0
#'round(t(sapply(nLD,
#'               function(n) {
#'                 res <- crossGeneral(tMax=1,h=1/n,cFct=cFT1p1,withBounds=TRUE,Lplus=TRUE)
#'                 c(n=n,Gl=res$Gl[n],Gu=res$Gu[n])
#'               }
#'               )
#'         ),
#'       digits=5)
#'
#'## Part 2: c(t) = t^2 + 3*t + 1
#'## Left columns pair: b(t) = 3 - 2*t
#'round(t(sapply(nLD,
#'               function(n) {
#'                 res <- crossGeneral(tMax=1,h=1/n,cFct=cFT1p3,bFct=function(t) 3-2*t,withBounds=TRUE,Lplus=TRUE)
#'                 c(n=n,Gl=res$Gl[n],Gu=res$Gu[n])
#'               }
#'               )
#'         ),
#'       digits=5)
#'
#'## Right columns pair: b(t) = 3 - t
#'round(t(sapply(nLD,
#'               function(n) {
#'                 res <- crossGeneral(tMax=1,h=1/n,cFct=cFT1p3,bFct=function(t) 3-2*t,withBounds=TRUE,Lplus=TRUE)
#'                 c(n=n,Gl=res$Gl[n],Gu=res$Gu[n])
#'               }
#'               )
#'         ),
#'       digits=5)
#'
#'## Part 3: c(t) = 1 + sin(t)
#'## Left columns pair: b(t) = c'(t)
#'round(t(sapply(nLD,
#'               function(n) {
#'                 res <- crossGeneral(tMax=1,h=1/n,cFct=function(t) 1+sin(t),bFct=function(t) cos(t),withBounds=TRUE,Lplus=TRUE)
#'                 c(n=n,Gl=res$Gl[n],Gu=res$Gu[n])
#'               }
#'               )
#'        ),
#'      digits=5)
#'
#'## Left columns pair: b(t) = 0.5
#'round(t(sapply(nLD,
#'               function(n) {
#'                 res <- crossGeneral(tMax=1,h=1/n,cFct=function(t) 1+sin(t),bFct=function(t) rep(0.5,length(t)),withBounds=TRUE,Lplus=TRUE)
#'                 c(n=n,Gl=res$Gl[n],Gu=res$Gu[n])
#'               }
#'               )
#'        ),
#'      digits=5)
#'
#'
#'## Check crossGeneral against an analytical inverse Gaussian
#'## distribution
#'## Define inverse Gaussian parameters
#'mu.true <- 0.075
#'sigma2.true <- 3
#'## Define a function transforming the "drift" (mu.true) and
#'## "noise variance" (sigma2.true) of the default inverse
#'## Gaussian parametrization of STAR into a
#'## linear boundary of an equivalent Wiener process first
#'## passage time problem
#'star2ld <- function(mu,sigma2) c(sqrt(1/sigma2),-sqrt(1/sigma2)/mu)
#'## Get the "equivalent" boundary parameters (y intercept and slope)
#'parB1 <- star2ld(mu.true,sigma2.true)
#'## Plot the "target" inverse Gaussian density
#'xx <- seq(0.001,0.3,0.001)
#'plot(xx,dinvgauss(xx,mu=mu.true,sigma2=sigma2.true),type="l")
#'## Get the numerical estimate of the density using Loader and
#'## Deely Volterra integral equation method
#'igB1 <- crossGeneral(tMax=0.3,h=0.001,cFct=function(t) parB1[1]+parB1[2]*t,withBounds=FALSE)
#'## superpose the numerical estimate to the exact solution
#'## use lines method to do that
#'lines(igB1,"density",col=2)
#'
#'## Use of crossTight and associated function
#'## Get the paramters, a and b, of the "approximate"
#'## tightest boundary: c(t) = a + b*sqrt(t), giving a 
#'## 0.05 probability of exit between 0 and 1
#'## (in fact we are discussing here a pair of symmetrical
#'## bounaries, c(t) and -c(t)). See Kendall et al (2007)
#'## for details
#'## Start by defining the target function
#'target95 <- mkTightBMtargetFct(ci=0.95)
#'## get the optimal log(a) and log(b) using
#'## the values of table 1 of Kendall et al as initial
#'## guesses
#'p95 <- optim(log(c(0.3,2.35)),target95,method="BFGS")
#'## check the convergence of BFGS
#'p95$convergence
#'## check if the parameters changed a lot
#'exp(p95$par)
#'## Get the bounds on G(1) for these optimal parameters
#'d95 <- crossTight(a=exp(p95$par[1]),b=exp(p95$par[2]),withBound=TRUE,logScale=FALSE)
#'## print out the summary
#'summary(d95)
#'## Do the same for the 0.01 probability of first passage
#'target99 <- mkTightBMtargetFct(ci=0.99)
#'p99 <- optim(p95$par,target99,method="BFGS")
#'p99$convergence
#'exp(p99$par)
#'d99 <- crossTight(a=exp(p99$par[1]),b=exp(p99$par[2]),withBound=TRUE,logScale=FALSE)
#'summary(d99) 
#'}
#'
#'
crossGeneral <- function(tMax=1,
                         h=0.001,
                         cFct,
                         cprimeFct,
                         bFct,
                         withBounds=FALSE,
                         Lplus
                         ) {

  ## Check arg. cFct. It should be a function or a numeric
  if (!is.function(cFct)) {
    if (!is.numeric(cFct)) {
      stop("cFct should be a function or a numeric.")
    } else {
      cCst <- cFct[1]
      cFct <- function(t) rep(cCst,length(t))
    } ## End of conditional on !is.numeric(cFct) 
  } ## End of conditional on !is.function(cFct) 

  n <- round(tMax/h)
  tMax <- n*h
  tt <- (1:n)*h
  
  ## Check arg. bFct.
  if (!missing(bFct)) {
    if (!is.function(bFct)) stop("bFct should be a function.")
    if (withBounds) {
      if (any(bFct(tt) < 0)) stop("bFct should be >= 0.")
      if (missing(Lplus) && !missing(cprimeFct)) {
        ## Do checks on bFct
        if (!is.function(cprimeFct)) stop("cprimeFct should be a function.")
        testFct <- function(u,t) 2*cprimeFct(u) - (cFct(t)-cFct(u))/(t-u)
        tests <- sapply(tt[-1],
                        function(t) {
                          b.t <- bFct(t)
                          u <- tt[tt<t]
                          tu <- testFct(u,t)
                          ub <- min(tu)
                          lb <- max(0,max(tu))
                          c(test3p9 = (b.t > ub),
                            test3p12 = (lb > b.t)
                            )
                        }
                        )
        tests <- apply(tests,1,any)
        if (tests[1]) {
          Lplus <- TRUE
        } else {
          if (tests[2]) {
            Lplus <- FALSE
          } else {
            stop("bFct not suitable for bounds calculation.")
          }
        } ## End of conditional of tests[1] 
      } ## End of conditional on !missing(cprimeFct)
    } ## End of conditional on withBounds
  } else {
    if (withBounds && missing(Lplus)) stop("Lplus required")
  } ## End of conditional on !missing(bFct)

  if (missing(bFct)) {
    fFct <- function(t) 2*pnorm(-cFct(t)/sqrt(t))
    kFct <- function(t,u) {
      numerator <- numeric(length(u))
      denominator <- numerator
      ratio <- numerator
      result <- numerator
      largeU <- u >= t
      numerator[!largeU] <- cFct(u[!largeU])-cFct(t)
      denominator[!largeU] <- sqrt(t-u[!largeU])
      ratio[!largeU] <- numerator[!largeU]/denominator[!largeU]
      ratio[!largeU][is.nan(ratio[!largeU])] <- 0
      result[!largeU] <- 2*pnorm(ratio[!largeU])
      result[u==t] <- 1
      result
    }
  } else {
    fFct <- function(t) {
      term1 <- pnorm(-cFct(t)/sqrt(t))
      term2 <- exp(-2*bFct(t)*(cFct(t)-t*bFct(t)))*
        pnorm((-cFct(t)+2*t*bFct(t))/sqrt(t))
      term1+term2
    }
    kFct <- function(t,u) {
      numerator1 <- numeric(length(u))
      numerator2 <- numerator1
      denominator <- numerator1
      ratio1 <- numerator1
      ratio2 <- numerator1
      term1 <- numerator1
      term2 <- numerator1
      largeU <- u >= t
      numerator1[!largeU] <- cFct(u[!largeU])-cFct(t)
      denominator[!largeU] <- sqrt(t-u[!largeU])
      ratio1[!largeU] <- numerator1[!largeU]/denominator[!largeU]
      ## ratio1[!largeU][is.nan(ratio1[!largeU])] <- 0
      term1[!largeU] <- pnorm(ratio1[!largeU])
      term1[u == t] <- 0.5
      numerator2[!largeU] <- cFct(u[!largeU])-cFct(t)+2*(t-u[!largeU])*bFct(t)
      ratio2[!largeU] <- numerator2[!largeU]/denominator[!largeU]
      ## ratio2[!largeU][is.nan(ratio2[!largeU])] <- 0
      term2[!largeU] <- exp(-2*bFct(t)*
                            (cFct(t)-cFct(u[!largeU])-
                             (t-u[!largeU])*bFct(t)))*pnorm(ratio2[!largeU])
      term2[u == t] <- 0.5
      term1+term2
    }
  } ## End of conditional on missing(bFct)

  ## We can now implement the mid-point method of 
  ## Loader & Deely (1987) pp 99-100
  
  mt <- tt-0.5*h
  B <- fFct(tt)
  A <- t(sapply(tt,function(x) kFct(x,mt)))
  x <- forwardsolve(A,B)
  x[x<0] <- 0
  
  result <- list(time=tt,
                 mids=mt,
                 G=cumsum(x),
                 g=x/h,
                 h=h
                 )

  class(result) <- "FirstPassageTime"

  if (withBounds) {
    n <- length(x)
    Gu <- numeric(n)
    Gl <- numeric(n)  
    if (Lplus) {
      ## First case for bounds considered by Loader and Deely
      ## Corresponding to Eq. 3.6 and 3.7
      f <- fFct(h)
      Gu[1] <- f/kFct(h,0)
      Gl[1] <- f

      for (i in 2:n) {
        f <- fFct(i*h)
        sumU <- f
        sumL <- f
        kjm <- kFct(i*h,0)
        kj <- kFct(i*h,h)
        for (j in 1:(i-1)) {
          if (j == i-1) {
            kjp <- 1
          } else {
            kjp <- kFct(i*h,(j+1)*h)
          }
          sumU <- sumU + (kj-kjm)*Gu[j]
          sumL <- sumL + (kjp-kj)*Gl[j]
          kjm <- kj
          kj <- kjp
        } ## End of the loop on j
        Gu[i] <- sumU/kjm
        Gl[i] <- sumL
      } ## End of the loop on i
    } else {
      ## Second case for bounds considered by Loader and Deely
      ## Corresponding to Eq. 3.10 and 3.11
      f <- fFct(h)
      Gu[1] <- f
      kjm <- kFct(h,0)
      kj <- 1
      Gl[1] <- f-Gu[1]*(kjm-kj)
      for (i in 2:n) {
        f <- fFct(i*h)
        sumU <- f
        sumL <- f
        kjm <- kFct(i*h,0)
        kj <- kFct(i*h,h)
        for (j in 1:(i-1)) {
          if (j == i-1) {
            kjp <- 1
          } else {
            kjp <- kFct(i*h,(j+1)*h)
          }
          sumU <- sumU - (kj-kjp)*Gl[j]
          sumL <- sumL - (kjm-kj)*Gu[j]
          kjm <- kj
          kj <- kjp
        } ## End of the loop on j
        Gu[i] <- sumU
        Gl[i] <- sumL
      } ## End of the loop on i
    } ## End of conditional on Lplus
    result$Gu <- Gu
    result$Gl <- Gl
  } ## End of conditional on withBounds
  
  result
  
}

plot.FirstPassageTime <- function(x,y,
                                  which=c("Distribution","density"),
                                  xlab,ylab,...) {
  if (missing(xlab)) xlab <- "Time (s)"
  if (missing(ylab)) {
    ylab <- switch(which[1],
                   Distribution="Probability",
                   density="Density (1/s)"
                   )
  }
  switch(which[1],
         Distribution=plot(x$time,x$G,xlab=xlab,ylab=ylab,type="l",...),
         density=plot(x$mids,x$g,xlab=xlab,ylab=ylab,type="l",...)
         )
  
}

lines.FirstPassageTime <- function(x,
                                   which=c("Distribution","density"),
                                   ...) {
  switch(which[1],
         Distribution=lines(x$time,x$G,...),
         density=lines(x$mids,x$g,...)
         )
  
}

summary.FirstPassageTime <- function(object,digits,...) {
  n <- length(object$time)
  tMax <- object$time[n]
  Gmax <- object$G[n]
  if ("Gu" %in% names(object)) {
    GuMax <- object$Gu[n]
    GlMax <- object$Gl[n]
    range <- GuMax - GlMax
    i <- 1
    rr <- range*10^i
    while (rr < 1) {
      i <- i+1
      rr <- range*10^i
    }
    cat(paste("Prob. of first passage before " ,tMax,": ",round(Gmax,digits=i),sep=""))
    cat(paste(" (bounds: [",round(GlMax,digits=i),",",round(GuMax,digits=i),"])",sep=""))
  } else {
    if (missing(digits)) digits <- 4
    cat(paste("Prob. of first passage before " ,tMax,": ",round(Gmax,digits=digits),sep=""))
  }
  cat("\n")
  cat(paste("Integration time step used: ",object$h,".\n",sep=""))
  
}

print.FirstPassageTime <- function(x,...) {

  n <- length(x$time)
  tMax <- x$time[n]
  Gmax <- x$G[n]
  result <- Gmax
  if ("Gu" %in% names(x)) {
    GuMax <- x$Gu[n]
    GlMax <- x$Gl[n]
    attr(result,"range") <- c(GlMax,GuMax)
  }

  print(result)
  
}

