#'Maximum Likelihood Parameter Estimation of an Inverse Gaussian Model with
#'Possibly Censored Data
#'
#'Estimate inverse Gaussian model parameters by the maximum likelihood method
#'using possibly censored data. Two different parameterizations of the inverse
#'Gaussian distribution can be used.
#'
#'The 2 different parameterizations of the inverse Gaussian distribution are
#'discussed in the manual of \code{\link{dinvgauss}}.
#'
#'In the absence of censored data the ML estimates are available in closed form
#'(Lindsey, 2004, p 212) together with the Hessian matrix at the MLE. In
#'presence of censored data an initial guess for the parameters is obtained
#'using the uncensored data before maximizing the likelihood function to the
#'full data set using \code{\link{optim}} with the \code{BFGS} method. ML
#'estimation is always performed with the \code{"sigma2"} parameterization.
#'Parameters and variance-covariance matrix are transformed at the end if the
#'\code{"boundary"} parameterization is requested.
#'
#'In order to ensure good behavior of the numerical optimization routines,
#'optimization is performed on the log of the parameters (\code{mu} and
#'\code{sigma2}).
#'
#'Standard errors are obtained from the inverse of the observed information
#'matrix at the MLE. They are transformed to go from the log scale used by the
#'optimization routine, when the latter is used (ie, for censored data) to the
#'parameterization requested.
#'
#'@param yi vector of (possibly binned) observations or a \code{spikeTrain}
#'object.
#'@param ni vector of counts for each value of \code{yi}; default:
#'\code{numeric(length(yi))+1}.
#'@param si vector of counts of \emph{uncensored} observations for each value
#'of \code{yi}; default: \code{numeric(length(yi))+1}.
#'@param parameterization parameterization used, \code{"sigma2"} (default) of
#'\code{"boundary"}.
#'@return A list of class \code{durationFit} with the following components:
#'@returnItem estimate the estimated parameters, a named vector.
#'@returnItem se the standard errors, a named vector.
#'@returnItem logLik the log likelihood at maximum.
#'@returnItem r a function returning the log of the relative likelihood
#'function.
#'@returnItem mll a function returning the opposite of the log likelihood
#'function using the log of the parameters.
#'@returnItem call the matched call.
#'@note The returned standard errors (component \code{se}) are valid in the
#'asymptotic limit. You should plot contours using function \code{r} in the
#'returned list and check that the contours are reasonably close to ellipses.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso
#'\code{\link{dinvgauss}},\code{\link{lnormMLE}},\code{\link{gammaMLE}},\code{\link{weibullMLE}},\code{\link{llogisMLE}},\code{\link{rexpMLE}}.
#'@references Lindsey, J.K. (2004) \emph{Introduction to Applied Statistics: A
#'Modelling Approach}. OUP.
#'@keywords distribution ts
#'@examples
#'
#'## Simulate sample of size 100 from an inverse Gaussian
#'## distribution
#'set.seed(1102006,"Mersenne-Twister")
#'sampleSize <- 100
#'mu.true <- 0.075 
#'sigma2.true <- 3
#'sampleSize <- 100
#'sampIG <- rinvgauss(sampleSize,mu=mu.true,sigma2=sigma2.true)
#'## Make a maximum likelihood fit
#'sampIGmleIG <- invgaussMLE(sampIG)
#'## Compare estimates with actual values
#'rbind(est = coef(sampIGmleIG),se = sampIGmleIG$se,true = c(mu.true,sigma2.true))
#'## In the absence of censoring the MLE of the inverse Gaussian is available in a
#'## closed form together with its variance (ie, the observed information matrix)
#'## we can check that we did not screw up at that stage by comparing the observed
#'## information matrix obtained numerically with the analytical one. To do that we
#'## use the MINUS log likelihood function returned by invgaussMLE to get a numerical
#'## estimate
#'detailedFit <- optim(par=as.vector(log(sampIGmleIG$estimate)),
#'                     fn=sampIGmleIG$mll,
#'                     method="BFGS",
#'                     hessian=TRUE)
#'## We should not forget that the "mll" function uses the log of the parameters while
#'## the "se" component of sampIGmleIG list is expressed on the linear scale we must therefore
#'## transform one into the other as follows (Kalbfleisch, 1985, p71):
#'## if x = exp(u) and y = exp(v) and if we have the information matrix in term of
#'## u and v (that's the hessian component of list detailedFit above), we create matrix:
#'##      du/dx du/dy
#'## Q =
#'##      dv/dx dv/dy
#'## and we get I in term of x and y by the following matrix product:
#'## I(x,y) <- t(Q) %*% I(u,v) %*% Q
#'## In the present case:
#'##  du/dx = 1/exp(u), du/dy = 0, dv/dx = 0, dv/dy = 1/exp(v)
#'## Therefore:
#'Q <- diag(1/exp(detailedFit$par))
#'numericalI <- t(Q) %*% detailedFit$hessian %*% Q
#'seComp <- rbind(sampIGmleIG$se, sqrt(diag(solve(numericalI))))
#'colnames(seComp) <- c("mu","sigma2")
#'rownames(seComp) <- c("analytical", "numerical")
#'seComp
#'## We can check the relative differences between the 2
#'apply(seComp,2,function(x) abs(diff(x))/x[1])
#'
#'\dontrun{
#'## Estimate the log relative likelihood on a grid to plot contours
#'Mu <- seq(coef(sampIGmleIG)[1]-4*sampIGmleIG$se[1],
#'          coef(sampIGmleIG)[1]+4*sampIGmleIG$se[1],
#'          sampIGmleIG$se[1]/10)
#'Sigma2 <- seq(coef(sampIGmleIG)[2]-4*sampIGmleIG$se[2],
#'              coef(sampIGmleIG)[2]+4*sampIGmleIG$se[2],
#'              sampIGmleIG$se[2]/10)
#'sampIGmleIGcontour <- sapply(Mu, function(mu) sapply(Sigma2, function(s2) sampIGmleIG$r(mu,s2)))
#'## plot contours using a linear scale for the parameters
#'## draw four contours corresponding to the following likelihood ratios:
#'##  0.5, 0.1, Chi2 with 2 df and p values of 0.95 and 0.99
#'X11(width=12,height=6)
#'layout(matrix(1:2,ncol=2))
#'contour(Mu,Sigma2,t(sampIGmleIGcontour),
#'        levels=c(log(c(0.5,0.1)),-0.5*qchisq(c(0.95,0.99),df=2)),
#'        labels=c("log(0.5)",
#'          "log(0.1)",
#'          "-1/2*P(Chi2=0.95)",
#'          "-1/2*P(Chi2=0.99)"),
#'        xlab=expression(mu),ylab=expression(sigma^2),
#'        main="Log Relative Likelihood Contours"
#'        )
#'points(coef(sampIGmleIG)[1],coef(sampIGmleIG)[2],pch=3)
#'points(mu.true,sigma2.true,pch=16,col=2)
#'## The contours are not really symmetrical about the MLE we can try to
#'## replot them using a log scale for the parameters to see if that improves
#'## the situation
#'contour(log(Mu),log(Sigma2),t(sampIGmleIGcontour),
#'        levels=c(log(c(0.5,0.1)),-0.5*qchisq(c(0.95,0.99),df=2)),
#'        labels="",
#'        xlab=expression(log(mu)),ylab=expression(log(sigma^2)),
#'        main="Log Relative Likelihood Contours",
#'        sub="log scale for the parameters")
#'points(log(coef(sampIGmleIG)[1]),log(coef(sampIGmleIG)[2]),pch=3)
#'points(log(mu.true),log(sigma2.true),pch=16,col=2)
#'
#'## Even with the log scale the contours are not ellipsoidal, so let us compute profiles
#'## For that we are going to use the returned MINUS log likelihood function
#'logMuProfFct <- function(logMu,...) {
#'  myOpt <- optimise(function(x) sampIGmleIG$mll(c(logMu,x))+logLik(sampIGmleIG),...)
#'  as.vector(unlist(myOpt[c("objective","minimum")]))
#'}
#'logMuProfCI <- function(logMu,
#'                        CI,
#'                        a=logS2Seq[1],
#'                        b=logS2Seq[length(logS2Seq)]) logMuProfFct(logMu,c(a,b))[1] - qchisq(CI,1)/2
#'
#'logS2ProfFct <- function(logS2,...) {
#'  myOpt <- optimise(function(x) sampIGmleIG$mll(c(x,logS2))+logLik(sampIGmleIG),...)
#'  as.vector(unlist(myOpt[c("objective","minimum")]))
#'}
#'logS2ProfCI <- function(logS2, CI,
#'                        a=logMuSeq[1],
#'                        b=logMuSeq[length(logMuSeq)]) logS2ProfFct(logS2,c(a,b))[1] - qchisq(CI,1)/2
#'
#'
#'## We compute profiles (on the log scale) eploxing +/- 3 times
#'## the se about the MLE
#'logMuSE <- sqrt(diag(solve(detailedFit$hessian)))[1]
#'logMuSeq <- seq(log(coef(sampIGmleIG)[1])-3*logMuSE,
#'                log(coef(sampIGmleIG)[1])+3*logMuSE,
#'                logMuSE/10)
#'logS2SE <- sqrt(diag(solve(detailedFit$hessian)))[2]
#'logS2Seq <- seq(log(coef(sampIGmleIG)[2])-3*logS2SE,
#'                log(coef(sampIGmleIG)[2])+3*logS2SE,
#'                logS2SE/10)
#'logMuProf <- sapply(logMuSeq,logMuProfFct,
#'                    lower=logS2Seq[1],
#'                    upper=logS2Seq[length(logS2Seq)])
#'## Get 95% and 99% CI
#'logMuCI95 <- c(uniroot(logMuProfCI,
#'                       interval=c(logMuSeq[1],log(coef(sampIGmleIG)[1])),
#'                       CI=0.95)$root,
#'               uniroot(logMuProfCI,
#'                       interval=c(log(coef(sampIGmleIG)[1]),logMuSeq[length(logMuSeq)]),
#'                       CI=0.95)$root
#'               )
#'logMuCI99 <- c(uniroot(logMuProfCI,
#'                       interval=c(logMuSeq[1],log(coef(sampIGmleIG)[1])),
#'                       CI=0.99)$root,
#'               uniroot(logMuProfCI,
#'                       interval=c(log(coef(sampIGmleIG)[1]),logMuSeq[length(logMuSeq)]),
#'                       CI=0.99)$root
#'               )
#'
#'logS2Prof <- sapply(logS2Seq,logS2ProfFct,
#'                    lower=logMuSeq[1],
#'                    upper=logMuSeq[length(logMuSeq)])
#'## Get 95% and 99% CI
#'logS2CI95 <- c(uniroot(logS2ProfCI,
#'                       interval=c(logS2Seq[1],log(coef(sampIGmleIG)[2])),
#'                       CI=0.95)$root,
#'               uniroot(logS2ProfCI,
#'                       interval=c(log(coef(sampIGmleIG)[2]),logS2Seq[length(logS2Seq)]),
#'                       CI=0.95)$root
#'               )
#'logS2CI99 <- c(uniroot(logS2ProfCI,
#'                       interval=c(logS2Seq[1],log(coef(sampIGmleIG)[2])),
#'                       CI=0.99)$root,
#'               uniroot(logS2ProfCI,
#'                       interval=c(log(coef(sampIGmleIG)[2]),logS2Seq[length(logS2Seq)]),
#'                       CI=0.99)$root
#'               )
#'
#'
#'## Add profiles to the previous plot
#'lines(logMuSeq,logMuProf[2,],col=2,lty=2)
#'lines(logS2Prof[2,],logS2Seq,col=2,lty=2)
#'
#'## We can now check the deviations of the (profiled) deviances
#'## from the asymptotic parabolic curves
#'X11()
#'layout(matrix(1:4,nrow=2))
#'oldpar <- par(mar=c(4,4,2,1))
#'logMuSeqOffset <- logMuSeq-log(coef(sampIGmleIG)[1])
#'logMuVar <- diag(solve(detailedFit$hessian))[1]
#'plot(logMuSeq,2*logMuProf[1,],type="l",xlab=expression(log(mu)),ylab="Deviance")
#'lines(logMuSeq,logMuSeqOffset^2/logMuVar,col=2)
#'points(log(coef(sampIGmleIG)[1]),0,pch=3)
#'abline(h=0)
#'abline(h=qchisq(0.95,1),lty=2)
#'abline(h=qchisq(0.99,1),lty=2)
#'lines(rep(logMuCI95[1],2),c(0,qchisq(0.95,1)),lty=2)
#'lines(rep(logMuCI95[2],2),c(0,qchisq(0.95,1)),lty=2)
#'lines(rep(logMuCI99[1],2),c(0,qchisq(0.99,1)),lty=2)
#'lines(rep(logMuCI99[2],2),c(0,qchisq(0.99,1)),lty=2)
#'## We can also "linearize" this last graph
#'plot(logMuSeq,
#'     sqrt(2*logMuProf[1,])*sign(logMuSeqOffset),
#'     type="l",
#'     xlab=expression(log(mu)),
#'     ylab=expression(paste("signed ",sqrt(Deviance)))
#'     )
#'lines(logMuSeq,
#'      sqrt(logMuSeqOffset^2/logMuVar)*sign(logMuSeqOffset),
#'      col=2)
#'points(log(coef(sampIGmleIG)[1]),0,pch=3)
#'
#'logS2SeqOffset <- logS2Seq-log(coef(sampIGmleIG)[2])
#'logS2Var <- diag(solve(detailedFit$hessian))[2]
#'plot(logS2Seq,2*logS2Prof[1,],type="l",xlab=expression(log(sigma^2)),ylab="Deviance")
#'lines(logS2Seq,logS2SeqOffset^2/logS2Var,col=2)
#'points(log(coef(sampIGmleIG)[2]),0,pch=3)
#'abline(h=0)
#'abline(h=qchisq(0.95,1),lty=2)
#'abline(h=qchisq(0.99,1),lty=2)
#'lines(rep(logS2CI95[1],2),c(0,qchisq(0.95,1)),lty=2)
#'lines(rep(logS2CI95[2],2),c(0,qchisq(0.95,1)),lty=2)
#'lines(rep(logS2CI99[1],2),c(0,qchisq(0.99,1)),lty=2)
#'lines(rep(logS2CI99[2],2),c(0,qchisq(0.99,1)),lty=2)
#'## We can also "linearize" this last graph
#'plot(logS2Seq,
#'     sqrt(2*logS2Prof[1,])*sign(logS2SeqOffset),
#'     type="l",
#'     xlab=expression(log(sigma^2)),
#'     ylab=expression(paste("signed ",sqrt(Deviance)))
#'     )
#'lines(logS2Seq,
#'      sqrt(logS2SeqOffset^2/logS2Var)*sign(logS2SeqOffset),
#'      col=2)
#'points(log(coef(sampIGmleIG)[2]),0,pch=3)
#'par(oldpar)
#'
#'## make a parametric boostrap to check the distribution of the deviance
#'nbReplicate <- 1000 #10000
#'sampleSize <- 100
#'system.time(
#'devianceIG100 <- lapply(1:nbReplicate,
#'                        function(idx) {
#'                          if ((idx %% 10 - 1) == 0) cat(paste("Doing now iteration:",idx,"\n"))
#'                          sampIG <- rinvgauss(sampleSize,mu=mu.true,sigma2=sigma2.true)
#'                          sampIGmleIG <- invgaussMLE(sampIG)
#'                          Deviance <- -2*sampIGmleIG$r(mu.true,sigma2.true)
#'                          logPara <- log(coef(sampIGmleIG))
#'                          logParaSE <- sampIGmleIG$se/coef(sampIGmleIG)
#'                          intervalMu <- function(n) c(-n,n)*logParaSE[1]+logPara[1]
#'                          intervalS2 <- function(n) c(-n,n)*logParaSE[2]+logPara[2]
#'                          logMuProfFct <- function(logMu,...) {
#'                            optimise(function(x)
#'                                     sampIGmleIG$mll(c(logMu,x))+logLik(sampIGmleIG),...)$objective
#'                          }
#'                          logMuProfCI <- function(logMu,
#'                                                  CI,
#'                                                  a=intervalS2(4)[1],
#'                                                  b=intervalS2(4)[2])
#'                            logMuProfFct(logMu,c(a,b)) - qchisq(CI,1)/2
#'                          
#'                          logS2ProfFct <- function(logS2,...) {
#'                            optimise(function(x)
#'                                     sampIGmleIG$mll(c(x,logS2))+logLik(sampIGmleIG),...)$objective
#'                          }
#'                          logS2ProfCI <- function(logS2, CI,
#'                                                  a=intervalMu(4)[1],
#'                                                  b=intervalMu(4)[2])
#'                            logS2ProfFct(logS2,c(a,b)) - qchisq(CI,1)/2
#'                          
#'                          factor <- 4
#'                          while((logMuProfCI(intervalMu(factor)[2],0.99) *
#'                                 logMuProfCI(logPara[1],0.99) >= 0) ||
#'                                (logMuProfCI(intervalMu(factor)[1],0.99) *
#'                                 logMuProfCI(logPara[1],0.99) >= 0)
#'                                ) factor <- factor+1
#'                          ##browser()
#'                          logMuCI95 <- c(uniroot(logMuProfCI,
#'                                                 interval=c(intervalMu(factor)[1],logPara[1]),
#'                                                 CI=0.95)$root,
#'                                         uniroot(logMuProfCI,
#'                                                 interval=c(logPara[1],intervalMu(factor)[2]),
#'                                                 CI=0.95)$root
#'                                         )
#'                          logMuCI99 <- c(uniroot(logMuProfCI,
#'                                                 interval=c(intervalMu(factor)[1],logPara[1]),
#'                                                 CI=0.99)$root,
#'                                         uniroot(logMuProfCI,
#'                                                 interval=c(logPara[1],intervalMu(factor)[2]),
#'                                                 CI=0.99)$root
#'                                         )
#'                          factor <- 4
#'                          while((logS2ProfCI(intervalS2(factor)[2],0.99) *
#'                                 logS2ProfCI(logPara[2],0.99) >= 0) ||
#'                                (logS2ProfCI(intervalS2(factor)[1],0.99) *
#'                                 logS2ProfCI(logPara[2],0.99) >= 0)
#'                                ) factor <- factor+1
#'                          logS2CI95 <- c(uniroot(logS2ProfCI,
#'                                                 interval=c(intervalS2(factor)[1],logPara[2]),
#'                                                 CI=0.95)$root,
#'                                         uniroot(logS2ProfCI,
#'                                                    interval=c(logPara[2],intervalS2(factor)[2]),
#'                                                 CI=0.95)$root
#'                                         )
#'                          logS2CI99 <- c(uniroot(logS2ProfCI,
#'                                                 interval=c(intervalS2(factor)[1],logPara[2]),
#'                                                 CI=0.99)$root,
#'                                         uniroot(logS2ProfCI,
#'                                                 interval=c(logPara[2],intervalS2(factor)[2]),
#'                                                 CI=0.99)$root
#'                                         )
#'                          list(deviance=Deviance,
#'                               logMuCI95=logMuCI95,
#'                               logMuNorm95=qnorm(c(0.025,0.975),logPara[1],logParaSE[1]),
#'                               logMuCI99=logMuCI99,
#'                               logMuNorm99=qnorm(c(0.005,0.995),logPara[1],logParaSE[1]),
#'                               logS2CI95=logS2CI95,
#'                               logS2Norm95=qnorm(c(0.025,0.975),logPara[2],logParaSE[2]),
#'                               logS2CI99=logS2CI99,
#'                               logS2Norm99=qnorm(c(0.005,0.995),logPara[2],logParaSE[2]))
#'                        }
#'                        )
#'            )[3]
#'## Find out how many times the true parameters was within the computed CIs
#'nLogMuCI95 <- sum(sapply(devianceIG100,
#'                         function(l) l$logMuCI95[1] <= log(mu.true)  &&
#'                         log(mu.true)<= l$logMuCI95[2]
#'                         )
#'                  )
#'nLogMuNorm95 <- sum(sapply(devianceIG100,
#'                           function(l) l$logMuNorm95[1] <= log(mu.true)  &&
#'                           log(mu.true)<= l$logMuNorm95[2]
#'                           )
#'                    )
#'nLogMuCI99 <- sum(sapply(devianceIG100,
#'                         function(l) l$logMuCI99[1] <= log(mu.true)  &&
#'                         log(mu.true)<= l$logMuCI99[2]
#'                         )
#'                  )
#'nLogMuNorm99 <- sum(sapply(devianceIG100,
#'                           function(l) l$logMuNorm99[1] <= log(mu.true)  &&
#'                           log(mu.true)<= l$logMuNorm99[2]
#'                           )
#'                    )
#'## Check if these counts are compatible with the nominal CIs
#'c(prof95Mu=nLogMuCI95,norm95Mu=nLogMuNorm95)
#'qbinom(c(0.005,0.995),nbReplicate,0.95)
#'c(prof95Mu=nLogMuCI99,norm95Mu=nLogMuNorm99)
#'qbinom(c(0.005,0.995),nbReplicate,0.99)
#'
#'nLogS2CI95 <- sum(sapply(devianceIG100,
#'                         function(l) l$logS2CI95[1] <= log(sigma2.true)  &&
#'                         log(sigma2.true)<= l$logS2CI95[2]
#'                         )
#'                  )
#'nLogS2Norm95 <- sum(sapply(devianceIG100,
#'                           function(l) l$logS2Norm95[1] <= log(sigma2.true)  &&
#'                           log(sigma2.true)<= l$logS2Norm95[2]
#'                           )
#'                    )
#'nLogS2CI99 <- sum(sapply(devianceIG100,
#'                         function(l) l$logS2CI99[1] <= log(sigma2.true)  &&
#'                         log(sigma2.true)<= l$logS2CI99[2]
#'                         )
#'                  )
#'nLogS2Norm99 <- sum(sapply(devianceIG100,
#'                           function(l) l$logS2Norm99[1] <= log(sigma2.true)  &&
#'                           log(sigma2.true)<= l$logS2Norm99[2]
#'                           )
#'                    )
#'## Check if these counts are compatible with the nominal CIs
#'c(prof95S2=nLogS2CI95,norm95S2=nLogS2Norm95)
#'qbinom(c(0.005,0.995),nbReplicate,0.95)
#'c(prof95S2=nLogS2CI99,norm95S2=nLogS2Norm99)
#'qbinom(c(0.005,0.995),nbReplicate,0.99)
#'
#'
#'## Get 95 and 99% confidence intervals for the QQ plot
#'ci <- sapply(1:nbReplicate,
#'                 function(idx) qchisq(qbeta(c(0.005,0.025,0.975,0.995),
#'                                            idx,
#'                                            nbReplicate-idx+1),
#'                                      df=2)
#'             )
#'## make QQ plot
#'X <- qchisq(ppoints(nbReplicate),df=2)
#'Y <- sort(sapply(devianceIG100,function(l) l$deviance))
#'X11()
#'plot(X,Y,type="n",
#'     xlab=expression(paste(chi[2]^2," quantiles")),
#'     ylab="MC quantiles",
#'     main="Deviance with true parameters after ML fit of IG data",
#'     sub=paste("sample size:", sampleSize,"MC replicates:", nbReplicate)
#'     )
#'abline(a=0,b=1)
#'lines(X,ci[1,],lty=2)
#'lines(X,ci[2,],lty=2)
#'lines(X,ci[3,],lty=2)
#'lines(X,ci[4,],lty=2)
#'lines(X,Y,col=2)
#'}
#'
invgaussMLE <- function (yi, ni = numeric(length(yi)) + 1,
                         si = numeric(length(yi)) + 1,
                         parameterization = "sigma2"
                         ) {

  ## check if yi is a spikeTrain object if yes take the "diff"
  if (inherits(yi,"spikeTrain")) yi <- diff(yi)
  if (any(yi < 0)) 
    stop("yi elements must be non-negative")
  yi <- as.numeric(yi)
  if (!identical(length(yi), length(ni))) 
    stop("yi and ni should have the same length")
  if (any(ni < 0)) 
    stop("ni elements must be non-negative")
  if (!identical(class(ni), "integer") && !identical(ni, round(ni))) 
    stop("ni should be a vector of positive integers")
  if (!identical(length(yi), length(si))) 
    stop("yi and si should have the same length")
  if (any(si < 0)) 
    stop("si elements must be non-negative")
  if (!identical(si, round(si))) 
    stop("si should be a vector of positive integers")
  if (any(si > ni)) 
    stop("si elements should not be greater than ni elements")
  minusLogLik <- function(p) {
    if (missing(p)) {
      txt <- paste("This function argument should be a 2 component vector:\n", 
                   "  component 1 is the log of the mu parameter,\n", 
                   "  component 2 is the log of the sigma2 parameter,\n", 
                   "using the 'sigma2' parameterization of the inverse Gaussian.\n")
      cat(txt)
    } else {
      mu <- exp(p[1])
      sigma2 <- exp(p[2])
      -(ifelse(s.dot > 0,
               sum(dinvgauss(yi[si > 0], mu, sigma2, log = TRUE) * si[si > 0], 0)
               ) +
        ifelse(c.dot > 0,
               sum(pinvgauss(yi[ci > 0], mu, sigma2, lower.tail = FALSE,log.p = TRUE) * ci[ci > 0]), 0))
    }
  }
  s.dot <- sum(si)
  n.dot <- sum(ni)
  ci <- ni - si
  c.dot <- sum(ci)
  if (s.dot == n.dot) {
    mu.hat <- weighted.mean(yi, ni)
    inv.mean <- weighted.mean(1/yi, ni)
    sigma2.hat <- inv.mean - 1/mu.hat
    estimate <- c(mu.hat, sigma2.hat)
    observedI <- matrix(c(n.dot/(sigma2.hat * mu.hat^3), 
                          0, 0, n.dot/sigma2.hat^2/2), nrow = 2, byrow = TRUE)
    se <- sqrt(diag(solve(observedI)))
    l <- -minusLogLik(log(estimate))
  } else {
    if (s.dot >= 10) {
      mu.hat <- weighted.mean(yi, si)
      inv.mean <- weighted.mean(1/yi, si)
      sigma2.hat <- inv.mean - 1/mu.hat
    } else {
      mu.hat <- weighted.mean(yi, ni)
      inv.mean <- weighted.mean(1/yi, ni)
      sigma2.hat <- inv.mean - 1/mu.hat
    }
    mleFit <- optim(par=log(c(mu.hat, sigma2.hat)),
                    fn=minusLogLik,
                    method="BFGS",
                    hessian = TRUE)
    estimate <- exp(mleFit$par)
    newVar <- (1/estimate) %o% (1/estimate)
    observedI <- mleFit$hessian * newVar
    se <- sqrt(diag(solve(observedI)))
    l <- -mleFit$value
  }
  if (parameterization == "sigma2") {
    names(estimate) <- c("mu", "sigma2")
    names(se) <- c("mu", "sigma2")
    rFct <- function(mu, sigma2) -minusLogLik(log(c(mu, sigma2))) - l
  } else {
    boundary.hat <- (sigma2.hat)^(-0.5)
    mu.hat <- mu.hat/boundary.hat
    if (s.dot == n.dot) {
      estimate <- c(mu.hat, boundary.hat)
      observedI <- observedI * matrix(c(boundary.hat^2, 
                                        -2/boundary.hat^2, -2/boundary.hat^2, 4/boundary.hat^6), 
                                      nrow = 2, byrow = TRUE)
      se <- se * c(1/boundary.hat, boundary.hat^3/2)
    } else {
      estimate <- c(estimate[1] * sqrt(estimate[2]), 1/sqrt(estimate[2]))
      newVar <- newVar * (c(estimate[2], -2/estimate[2]^3) %o% 
                          c(estimate[2], -2/estimate[2]^3))
      observedI <- mleFit$hessian * newVar
      se <- sqrt(diag(solve(observedI)))
    }
    names(estimate) <- c("mu", "boundary")
    names(se) <- c("mu", "boundary")
    rFct <- function(mu, boundary) -minusLogLik(log(c(mu * boundary, 1/boundary^2))) - l
  }
  result <- list(estimate = estimate,
                 se = se,
                 logLik = l,
                 r = rFct, 
                 mll = minusLogLik,
                 call = match.call()
                 )
  class(result) <- "durationFit"
  return(result)

}



#'Utility Functions for durationFit Objects
#'
#'\code{coef.durationFit} and \code{logLik.durationFit} extract components of a
#'\code{durationFit} object, while \code{is.durationFit} tests if its argument
#'is such an object.
#'
#'Everything is trivial here.
#'
#'@aliases coef.durationFit logLik.durationFit is.durationFit
#'@param object a \code{durationFit} object.
#'@param obj an object to be tested against a \code{durationFit} object.
#'@param \dots see \code{\link{coef}} and \code{\link{logLik}}.
#'@return \code{coef.durationFit} returns the coefficients or the estimates or
#'the fitted parameters of the object: a 2 elements named vector.
#'
#'\code{logLik.durationFit} returns the loglikelihood value.
#'
#'\code{is.durationFit} returns \code{TRUE} if its argument is a
#'\code{durationFit} object and \code{FALSE} otherwise.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{compModels}}, \code{\link{invgaussMLE}},
#'\code{\link{lnormMLE}}, \code{\link{llogisMLE}}, \code{\link{rexpMLE}},
#'\code{\link{gammaMLE}}, \code{\link{weibullMLE}}
#'@keywords distribution ts
#'@examples
#'
#'\dontrun{
#'## load CAL1S data
#'data(CAL1S)
#'## convert the data into spikeTrain objects
#'CAL1S <- lapply(CAL1S,as.spikeTrain)
#'## look at the train of the 1sd neuron
#'CAL1S[["neuron 1"]]
#'## fit a invgauss model to the 1st neuron spike train
#'n1SDFig <- invgaussMLE(CAL1S[["neuron 1"]])
#'is.durationFit(n1SDFig)
#'coef(n1SDFig)
#'logLik(n1SDFig)
#'}
#'
coef.durationFit <- function(object,...) object$estimate
logLik.durationFit <- function(object,...) object$logLik



#'Maximum Likelihood Parameter Estimation of a Gamma Model with Possibly
#'Censored Data
#'
#'Estimate Gamma model parameters by the maximum likelihood method using
#'possibly censored data. Two different parameterizations of the Gamma
#'distribution can be used.
#'
#'There is no closed form expression for the MLE of a Gamma distribution. The
#'numerical method implemented here uses the profile likelihood described by
#'Monahan (2001) pp 210-216.
#'
#'In order to ensure good behavior of the numerical optimization routines,
#'optimization is performed on the log of the parameters (\code{shape} and
#'\code{scale} or \code{rate}).
#'
#'Standard errors are obtained from the inverse of the observed information
#'matrix at the MLE. They are transformed to go from the log scale used by the
#'optimization routine to the parameterization requested.
#'
#'@param yi vector of (possibly binned) observations or a \code{spikeTrain}
#'object.
#'@param ni vector of counts for each value of \code{yi}; default:
#'\code{numeric(length(yi))+1}.
#'@param si vector of counts of \emph{uncensored} observations for each value
#'of \code{yi}; default: \code{numeric(length(yi))+1}.
#'@param scale logical should the scale (\code{TRUE}) or the rate
#'parameterization (\code{FALSE}) be used?
#'@return A list of class \code{durationFit} with the following components:
#'@returnItem estimate the estimated parameters, a named vector.
#'@returnItem se the standard errors, a named vector.
#'@returnItem logLik the log likelihood at maximum.
#'@returnItem r a function returning the log of the relative likelihood
#'function.
#'@returnItem mll a function returning the opposite of the log likelihood
#'function using the log of the parameters.
#'@returnItem call the matched call.
#'@note The returned standard errors (component \code{se}) are valid in the
#'asymptotic limit. You should plot contours using function \code{r} in the
#'returned list and check that the contours are reasonably close to ellipses.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{GammaDist}}, \code{\link{invgaussMLE}},
#'\code{\link{lnormMLE}}
#'@references Monahan, J. F. (2001) \emph{Numerical Methods of Statistics}.
#'CUP.
#'
#'Lindsey, J.K. (2004) \emph{Introduction to Applied Statistics: A Modelling
#'Approach}. OUP.
#'@keywords distribution ts
#'@examples
#'
#'\dontrun{
#'## Simulate sample of size 100 from a gamma distribution
#'set.seed(1102006,"Mersenne-Twister")
#'sampleSize <- 100
#'shape.true <- 6
#'scale.true <- 0.012
#'sampGA <- rgamma(sampleSize,shape=shape.true,scale=scale.true)
#'sampGAmleGA <- gammaMLE(sampGA)
#'rbind(est = sampGAmleGA$estimate,se = sampGAmleGA$se,true = c(shape.true,scale.true))
#'
#'## Estimate the log relative likelihood on a grid to plot contours
#'Shape <- seq(sampGAmleGA$estimate[1]-4*sampGAmleGA$se[1],
#'               sampGAmleGA$estimate[1]+4*sampGAmleGA$se[1],
#'               sampGAmleGA$se[1]/10)
#'Scale <- seq(sampGAmleGA$estimate[2]-4*sampGAmleGA$se[2],
#'             sampGAmleGA$estimate[2]+4*sampGAmleGA$se[2],
#'             sampGAmleGA$se[2]/10)
#'sampGAmleGAcontour <- sapply(Shape, function(sh) sapply(Scale, function(sc) sampGAmleGA$r(sh,sc)))
#'## plot contours using a linear scale for the parameters
#'## draw four contours corresponding to the following likelihood ratios:
#'##  0.5, 0.1, Chi2 with 2 df and p values of 0.95 and 0.99
#'X11(width=12,height=6)
#'layout(matrix(1:2,ncol=2))
#'contour(Shape,Scale,t(sampGAmleGAcontour),
#'        levels=c(log(c(0.5,0.1)),-0.5*qchisq(c(0.95,0.99),df=2)),
#'        labels=c("log(0.5)",
#'          "log(0.1)",
#'          "-1/2*P(Chi2=0.95)",
#'          "-1/2*P(Chi2=0.99)"),
#'        xlab="shape",ylab="scale",
#'        main="Log Relative Likelihood Contours"
#'        )
#'points(sampGAmleGA$estimate[1],sampGAmleGA$estimate[2],pch=3)
#'points(shape.true,scale.true,pch=16,col=2)
#'## The contours are not really symmetrical about the MLE we can try to
#'## replot them using a log scale for the parameters to see if that improves
#'## the situation
#'contour(log(Shape),log(Scale),t(sampGAmleGAcontour),
#'        levels=c(log(c(0.5,0.1)),-0.5*qchisq(c(0.95,0.99),df=2)),
#'        labels="",
#'        xlab="log(shape)",ylab="log(scale)",
#'        main="Log Relative Likelihood Contours",
#'        sub="log scale for the parameters")
#'points(log(sampGAmleGA$estimate[1]),log(sampGAmleGA$estimate[2]),pch=3)
#'points(log(shape.true),log(scale.true),pch=16,col=2)
#'
#'## make a parametric boostrap to check the distribution of the deviance
#'nbReplicate <- 10000
#'sampleSize <- 100
#'system.time(
#'            devianceGA100 <- replicate(nbReplicate,{
#'                             sampGA <- rgamma(sampleSize,shape=shape.true,scale=scale.true)
#'                             sampGAmleGA <- gammaMLE(sampGA)
#'                             -2*sampGAmleGA$r(shape.true,scale.true)
#'                           }
#'                                       )
#'            )[3]
#'
#'## Get 95 and 99% confidence intervals for the QQ plot
#'ci <- sapply(1:nbReplicate,
#'                 function(idx) qchisq(qbeta(c(0.005,0.025,0.975,0.995),
#'                                            idx,
#'                                            nbReplicate-idx+1),
#'                                      df=2)
#'             )
#'## make QQ plot
#'X <- qchisq(ppoints(nbReplicate),df=2)
#'Y <- sort(devianceGA100)
#'X11()
#'plot(X,Y,type="n",
#'     xlab=expression(paste(chi[2]^2," quantiles")),
#'     ylab="MC quantiles",
#'     main="Deviance with true parameters after ML fit of gamma data",
#'     sub=paste("sample size:", sampleSize,"MC replicates:", nbReplicate)
#'     )
#'abline(a=0,b=1)
#'lines(X,ci[1,],lty=2)
#'lines(X,ci[2,],lty=2)
#'lines(X,ci[3,],lty=2)
#'lines(X,ci[4,],lty=2)
#'lines(X,Y,col=2)
#'}
#'
gammaMLE <- function(yi,
                     ni = numeric(length(yi))+1,
                     si = numeric(length(yi))+1,
                     scale = TRUE
                     ) {
  ## check if yi is a spikeTrain object if yes take the "diff"
  if (inherits(yi,"spikeTrain")) yi <- diff(yi)
  ## check that yi elements are positive
  if (any(yi < 0)) stop("yi elements must be non-negative")
  ## coerce yi to vector
  yi <- as.numeric(yi)
  ## check that ni has the same length as yi
  if (!identical(length(yi),length(ni))) stop("yi and ni should have the same length")
  ## check that the elements of ni are non negative and integer
  if (any(ni < 0)) stop("ni elements must be non-negative")
  if (!identical(ni, round(ni))) stop("ni should be a vector of positive integers")
  if (!identical(length(yi),length(si))) stop("yi and si should have the same length")
  ## check that the elements of ni are non negative and integer
  if (any(si < 0)) stop("si elements must be non-negative")
  if (!identical(si, round(si))) stop("si should be a vector of positive integers")
  if (any(si > ni)) stop("si elements should not be greater than ni elements")
  
  ## Create function returning the opposite of the log likelihood
  minusLogLik <- function(p) {
    if (missing(p)) {
      txt <- paste("This function argument should be a 2 component vector:\n",
                   "  component 1 is the log of the shape parameter,\n",
                   "  component 2 is the log of the scale parameter.\n"
                   )
      cat(txt)
    } else {
      shape <- exp(p[1])
      scale <- exp(p[2])
      -(ifelse(s.dot>0,sum(dgamma(yi[si>0],shape=shape,scale=scale,log=TRUE)*si[si>0],0))+
        ifelse(c.dot>0,sum(pgamma(yi[ci>0],shape=shape,
                                  scale=scale,lower.tail=FALSE,log.p=TRUE)*ci[ci>0]),0)
        )
    }
  }

  ## Get the number of uncensored events
  s.dot <- sum(si)
  ## Get the total number of events
  n.dot <- sum(ni)
  ci <- ni-si
  c.dot <- sum(ci)

  ## Define function logProfileLik returning the log of the
  ## profiled likelihood for the shape parameter
  logProfileLik <- function(shape) {
    if (shape <= 0) return(-Inf)
    term1 <- s.dot*shape*(log(shape) - log(scale.hat) -1)
    term2 <- (shape - 1) * sum(ni * log(yi))
    term3 <- - s.dot * lgamma(shape)
    return(term1 + term2 + term3)
  }
  
  if (s.dot == n.dot) {
    ## no censored event
    ## Get the empirical mean which is the MLE of parameter scale
    scale.hat <- weighted.mean(yi, ni)
    s2 <- weighted.mean(yi^2, ni) - scale.hat^2
    shape.hat <- scale.hat^2 / s2
    shape.lower <- floor(shape.hat)
    shape.upper <- ceiling(shape.hat)
    shape.hat <- optimize(logProfileLik,
                          lower = shape.lower,
                          upper = shape.upper,
                          maximum = TRUE)$maximum
    ## mleFit <- nlm(minusLogLik,log(c(shape.hat,scale.hat)),hessian=TRUE)
    mleFit <- optim(par=log(c(shape.hat,scale.hat)),
                    fn=minusLogLik,
                    method="BFGS",
                    hessian=TRUE)
  } else {
    ## some censored events
    ## if more than 10 events are uncesored get inital guess from them
    ## otherwise use all events
    if (s.dot >= 10) {
      scale.hat <- weighted.mean(yi, si)
      s2 <- weighted.mean(yi^2, si) - scale.hat^2
      shape.hat <- scale.hat^2 / s2
      shape.lower <- floor(shape.hat)
      shape.upper <- ceiling(shape.hat)
      shape.hat <- optimize(logProfileLik,
                            lower = shape.lower,
                            upper = shape.upper,
                            maximum = TRUE)$maximum
    } else {
      scale.hat <- weighted.mean(yi, ni)
      s2 <- weighted.mean(yi^2, ni) - scale.hat^2
      shape.hat <- scale.hat^2 / s2
    }
    ## mleFit <- nlm(minusLogLik,log(c(shape.hat,scale.hat)),hessian=TRUE)
    mleFit <- optim(par=log(c(shape.hat,scale.hat)),
                    fn=minusLogLik,
                    method="BFGS",
                    hessian=TRUE)
  } ## End of conditional on s.dot == n.dot
  ## estimate <- exp(mleFit$estimate)
  estimate <- exp(mleFit$par)
  newVar <- (1/estimate) %o% (1/estimate)
  observedI <- mleFit$hessian * newVar
  se <- sqrt(diag(solve(observedI)))
  ## l <- -mleFit$minimum
  l <- -mleFit$value


  if (scale) {
    names(estimate) <- c("shape","scale")
    names(se) <- c("shape","scale")
    rFct <- function(shape,scale) -minusLogLik(log(c(shape,scale))) - l
  } else {
    estimate <- c(estimate[1],1/estimate[2])
    names(estimate) <- c("shape","rate")
    se <- se * c(1,estimate[2]^2)
    names(se) <- c("shape","rate")
    rFct <- function(shape,rate) -minusLogLik(log(c(shape,1/rate))) - l
  }

  result <- list(estimate = estimate,
                 se = se,
                 logLik = l,
                 r = rFct,
                 mll = minusLogLik,
                 call = match.call()
                 )
  class(result) <- "durationFit"
  return(result)
  
}




#'Maximum Likelihood Parameter Estimation of a Log Logistic Model with Possibly
#'Censored Data
#'
#'Estimate log logistic model parameters by the maximum likelihood method using
#'possibly censored data.
#'
#'The MLE for the log logistic is not available in closed formed and is
#'therefore obtained numerically obtained by calling \code{\link{optim}} with
#'the \code{BFGS} method.
#'
#'In order to ensure good behavior of the numerical optimization routines,
#'optimization is performed on the log of parameter \code{scale}.
#'
#'Standard errors are obtained from the inverse of the observed information
#'matrix at the MLE. They are transformed to go from the log scale used by the
#'optimization routine to the requested parameterization.
#'
#'@param yi vector of (possibly binned) observations or a \code{spikeTrain}
#'object.
#'@param ni vector of counts for each value of \code{yi}; default:
#'\code{numeric(length(yi))+1}.
#'@param si vector of counts of \emph{uncensored} observations for each value
#'of \code{yi}; default: \code{numeric(length(yi))+1}.
#'@return A list of class \code{durationFit} with the following components:
#'@returnItem estimate the estimated parameters, a named vector.
#'@returnItem se the standard errors, a named vector.
#'@returnItem logLik the log likelihood at maximum.
#'@returnItem r a function returning the log of the relative likelihood
#'function.
#'@returnItem mll a function returning the opposite of the log likelihood
#'function using the log of parameter \code{sdlog}.
#'@returnItem call the matched call.
#'@note The returned standard errors (component \code{se}) are valid in the
#'asymptotic limit. You should plot contours using function \code{r} in the
#'returned list and check that the contours are reasonably close to ellipses.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{dllogis}}, \code{\link{invgaussMLE}},
#'\code{\link{gammaMLE}}, \code{\link{weibullMLE}}, \code{\link{rexpMLE}},
#'\code{\link{lnormMLE}}
#'@references Lindsey, J.K. (2004) \emph{Introduction to Applied Statistics: A
#'Modelling Approach}. OUP.
#'
#'Lindsey, J.K. (2004) \emph{The Statistical Analysis of Stochastic Processes
#'in Time}. CUP.
#'@keywords distribution ts
#'@examples
#'
#'\dontrun{
#'## Simulate sample of size 100 from a log logisitic
#'## distribution
#'set.seed(1102006,"Mersenne-Twister")
#'sampleSize <- 100
#'location.true <- -2.7
#'scale.true <- 0.025
#'sampLL <- rllogis(sampleSize,location=location.true,scale=scale.true)
#'sampLLmleLL <- llogisMLE(sampLL)
#'rbind(est = sampLLmleLL$estimate,se = sampLLmleLL$se,true = c(location.true,scale.true))
#'
#'## Estimate the log relative likelihood on a grid to plot contours
#'Loc <- seq(sampLLmleLL$estimate[1]-4*sampLLmleLL$se[1],
#'               sampLLmleLL$estimate[1]+4*sampLLmleLL$se[1],
#'               sampLLmleLL$se[1]/10)
#'Scale <- seq(sampLLmleLL$estimate[2]-4*sampLLmleLL$se[2],
#'             sampLLmleLL$estimate[2]+4*sampLLmleLL$se[2],
#'             sampLLmleLL$se[2]/10)
#'sampLLmleLLcontour <- sapply(Loc, function(m) sapply(Scale, function(s) sampLLmleLL$r(m,s)))
#'## plot contours using a linear scale for the parameters
#'## draw four contours corresponding to the following likelihood ratios:
#'##  0.5, 0.1, Chi2 with 2 df and p values of 0.95 and 0.99
#'X11(width=12,height=6)
#'layout(matrix(1:2,ncol=2))
#'contour(Loc,Scale,t(sampLLmleLLcontour),
#'        levels=c(log(c(0.5,0.1)),-0.5*qchisq(c(0.95,0.99),df=2)),
#'        labels=c("log(0.5)",
#'          "log(0.1)",
#'          "-1/2*P(Chi2=0.95)",
#'          "-1/2*P(Chi2=0.99)"),
#'        xlab="Location",ylab="Scale",
#'        main="Log Relative Likelihood Contours"
#'        )
#'points(sampLLmleLL$estimate[1],sampLLmleLL$estimate[2],pch=3)
#'points(location.true,scale.true,pch=16,col=2)
#'## The contours are not really symmetrical about the MLE we can try to
#'## replot them using a log scale for the parameters to see if that improves
#'## the situation
#'contour(Loc,log(Scale),t(sampLLmleLLcontour),
#'        levels=c(log(c(0.5,0.1)),-0.5*qchisq(c(0.95,0.99),df=2)),
#'        labels="",
#'        xlab="log(Location)",ylab="log(Scale)",
#'        main="Log Relative Likelihood Contours",
#'        sub="log scale for parameter: scale")
#'points(sampLLmleLL$estimate[1],log(sampLLmleLL$estimate[2]),pch=3)
#'points(location.true,log(scale.true),pch=16,col=2)
#'
#'## make a parametric boostrap to check the distribution of the deviance
#'nbReplicate <- 10000
#'sampleSize <- 100
#'system.time(
#'            devianceLL100 <- replicate(nbReplicate,{
#'              sampLL <- rllogis(sampleSize,location=location.true,scale=scale.true)
#'              sampLLmleLL <- llogisMLE(sampLL)
#'              -2*sampLLmleLL$r(location.true,scale.true)
#'            }
#'                                       )
#'            )[3]
#'
#'## Get 95 and 99% confidence intervals for the QQ plot
#'ci <- sapply(1:nbReplicate,
#'                 function(idx) qchisq(qbeta(c(0.005,0.025,0.975,0.995),
#'                                            idx,
#'                                            nbReplicate-idx+1),
#'                                      df=2)
#'             )
#'## make QQ plot
#'X <- qchisq(ppoints(nbReplicate),df=2)
#'Y <- sort(devianceLL100)
#'X11()
#'plot(X,Y,type="n",
#'     xlab=expression(paste(chi[2]^2," quantiles")),
#'     ylab="MC quantiles",
#'     main="Deviance with true parameters after ML fit of log logistic data",
#'     sub=paste("sample size:", sampleSize,"MC replicates:", nbReplicate)
#'     )
#'abline(a=0,b=1)
#'lines(X,ci[1,],lty=2)
#'lines(X,ci[2,],lty=2)
#'lines(X,ci[3,],lty=2)
#'lines(X,ci[4,],lty=2)
#'lines(X,Y,col=2)
#'}
#'
llogisMLE <- function(yi,
                      ni = numeric(length(yi))+1,
                      si = numeric(length(yi))+1
                      ) {

  ## check if yi is a spikeTrain object if yes take the "diff"
  if (inherits(yi,"spikeTrain")) yi <- diff(yi)
  ## check that yi elements are positive
  if (any(yi <= 0)) stop("yi elements must be positive or null")
  ## coerce yi to vector
  yi <- as.numeric(yi)
  ## check that ni has the same length as yi
  if (!identical(length(yi),length(ni))) stop("yi and ni should have the same length")
  ## check that the elements of ni are non negative and integer
  if (any(ni < 0)) stop("ni elements must be non-negative")
  if (!identical(ni, round(ni))) stop("ni should be a vector of positive integers")
  if (!identical(length(yi),length(si))) stop("yi and si should have the same length")
  ## check that the elements of ni are non negative and integer
  if (any(si < 0)) stop("si elements must be non-negative")
  if (!identical(si, round(si))) stop("si should be a vector of positive integers")
  if (any(si > ni)) stop("si elements should not be greater than ni elements")

  ## Create function returning the opposite of the log likelihood
  minusLogLik <- function(p) {
    if (missing(p)) {
      txt <- paste("This function argument should be a 2 component vector:\n",
                   "  component 1 is the location parameter,\n",
                   "  component 2 is the log of the scale parameter.\n"
                   )
      cat(txt)
    } else {
      location <- p[1]
      scale <- exp(p[2])
      -(ifelse(s.dot>0,sum(dllogis(yi[si>0],location,scale,log=TRUE)*si[si>0],0))+
        ifelse(c.dot>0,sum(pllogis(yi[ci>0],location,scale,lower.tail=FALSE,log.p=TRUE)*ci[ci>0]),0)
        )
    }
  }

  ## Get the number of uncensored events
  s.dot <- sum(si)
  ## Get the total number of events
  n.dot <- sum(ni)
  ci <- ni-si
  c.dot <- sum(ci)

  lyi <- log(yi)
  if (s.dot == n.dot) {
    ## no censored event
    ## Get the empirical mean which is the MLE of parameter mu
    location.moment <- weighted.mean(lyi, si)
    scale.moment <- sqrt(3*(weighted.mean(lyi^2, si) - location.moment^2))/pi
  } else {
    ## some censored events
    ## if more than 10 events are uncesored get inital guess from them
    ## otherwise use all events
    if (s.dot >= 10) {
      location.moment <- weighted.mean(lyi, si)
      scale.moment <- sqrt(3*(weighted.mean(lyi^2, si) - location.moment^2))/pi
    } else {
      location.moment <- weighted.mean(lyi, ni)
      scale.moment <- sqrt(3*(weighted.mean(lyi^2, ni) - location.moment^2))/pi
    }
  }
  ## mleFit <- nlm(minusLogLik,c(location.moment,log(scale.moment)),hessian=TRUE)
  mleFit <- optim(fn=minusLogLik,
                  par=c(location.moment,log(scale.moment)),
                  method="BFGS",
                  hessian=TRUE)
  ## estimate <- c(mleFit$estimate[1],exp(mleFit$estimate[2]))
  estimate <- c(mleFit$par[1],exp(mleFit$par[2]))
  newVar <- c(1,1/estimate[2]) %o% c(1,1/estimate[2])
  se <- sqrt(diag(solve(mleFit$hessian * newVar)))
  ## l <- -mleFit$minimum
  l <- -mleFit$value
  names(estimate) <- c("location","scale")
  names(se) <- c("location","scale")
  rFct <- function(location,scale) -minusLogLik(c(location,log(scale))) - l
  
  result <- list(estimate = estimate,
                 se = se,
                 logLik = l,
                 r = rFct,
                 mll = minusLogLik,
                 call = match.call()
                 )
  class(result) <- "durationFit"
  return(result)

}




#'Maximum Likelihood Parameter Estimation of a Log Normal Model with Possibly
#'Censored Data
#'
#'Estimate log normal model parameters by the maximum likelihood method using
#'possibly censored data.
#'
#'In the absence of censored data the ML estimates are available in closed form
#'together with the Hessian matrix at the MLE. In presence of censored data an
#'initial guess for the parameters is obtained using the uncensored data before
#'maximizing the likelihood function to the full data set using
#'\code{\link{optim}} with the \code{BFGS} method.
#'
#'In order to ensure good behavior of the numerical optimization routines,
#'optimization is performed on the log of parameter \code{sdlog}.
#'
#'Standard errors are obtained from the inverse of the observed information
#'matrix at the MLE. They are transformed to go from the log scale used by the
#'optimization routine, when the latter is used (ie, for censored data) to the
#'parameterization requested.
#'
#'@param yi vector of (possibly binned) observations or a \code{spikeTrain}
#'object.
#'@param ni vector of counts for each value of \code{yi}; default:
#'\code{numeric(length(yi))+1}.
#'@param si vector of counts of \emph{uncensored} observations for each value
#'of \code{yi}; default: \code{numeric(length(yi))+1}.
#'@return A list of class \code{durationFit} with the following components:
#'@returnItem estimate the estimated parameters, a named vector.
#'@returnItem se the standard errors, a named vector.
#'@returnItem logLik the log likelihood at maximum.
#'@returnItem r a function returning the log of the relative likelihood
#'function.
#'@returnItem mll a function returning the opposite of the log likelihood
#'function using the log of parameter \code{sdlog}.
#'@returnItem call the matched call.
#'@note The returned standard errors (component \code{se}) are valid in the
#'asymptotic limit. You should plot contours using function \code{r} in the
#'returned list and check that the contours are reasonably close to ellipses.
#'@author Christophe Pouzat \email{christophe.pouzat@@univ-paris5.fr}
#'@seealso \code{\link{Lognormal}},\code{\link{invgaussMLE}}
#'@references Lindsey, J.K. (2004) \emph{Introduction to Applied Statistics: A
#'Modelling Approach}. OUP.
#'@keywords distribution ts
#'@examples
#'
#'## Simulate sample of size 100 from a log normal
#'## distribution
#'set.seed(1102006,"Mersenne-Twister")
#'sampleSize <- 100
#'meanlog.true <- -2.4
#'sdlog.true <- 0.4
#'sampLN <- rlnorm(sampleSize,meanlog.true,sdlog.true)
#'sampLNmleLN <- lnormMLE(sampLN)
#'rbind(est = sampLNmleLN$estimate,se = sampLNmleLN$se,true = c(meanlog.true,sdlog.true))
#'## In the absence of censoring the MLE of the log normal is available in a
#'## closed form together with its variance (ie, the observed information matrix)
#'## we can check that we did not screw up at that stage by comparing the observed
#'## information matrix obtained numerically with the analytical one. To do that we
#'## use the MINUS log likelihood function returned by lnormMLE to get a numerical
#'## estimate
#'detailedFit <- optim(fn=sampLNmleLN$mll,
#'                     par=as.vector(c(sampLNmleLN$estimate[1],log(sampLNmleLN$estimate[2]))),
#'                     method="BFGS",
#'                     hessian=TRUE)
#'## We should not forget that the "mll" function uses the log of the sdlog parameter while
#'## the "se" component of sampLNmleLN list is expressed on the linear scale we must therefore
#'## transform one into the other as follows (Kalbfleisch, 1985, p71):
#'## if x = u and y = exp(v) and if we have the information matrix in term of
#'## u and v (that's the hessian component of list detailedFit above), we create matrix:
#'##      du/dx du/dy
#'## Q =
#'##      dv/dx dv/dy
#'## and we get I in term of x and y by the following matrix product:
#'## I(x,y) <- t(Q) %*% I(u,v) %*% Q
#'## In the present case:
#'##  du/dx = 1, du/dy = 0, dv/dx = 0, dv/dy = 1/exp(v)
#'## Therefore:
#'Q <- diag(c(1,1/exp(detailedFit$par[2])))
#'numericalI <- t(Q) %*% detailedFit$hessian %*% Q
#'seComp <- rbind(sampLNmleLN$se, sqrt(diag(solve(numericalI))))
#'colnames(seComp) <- c("meanlog","sdlog")
#'rownames(seComp) <- c("analytical", "numerical")
#'seComp
#'## We can check the relative differences between the 2
#'apply(seComp,2,function(x) abs(diff(x))/x[1])
#'
#'\dontrun{
#'## Estimate the log relative likelihood on a grid to plot contours
#'MeanLog <- seq(sampLNmleLN$estimate[1]-4*sampLNmleLN$se[1],
#'               sampLNmleLN$estimate[1]+4*sampLNmleLN$se[1],
#'               sampLNmleLN$se[1]/10)
#'SdLog <- seq(sampLNmleLN$estimate[2]-4*sampLNmleLN$se[2],
#'             sampLNmleLN$estimate[2]+4*sampLNmleLN$se[2],
#'             sampLNmleLN$se[2]/10)
#'sampLNmleLNcontour <- sapply(MeanLog, function(mu) sapply(SdLog, function(s) sampLNmleLN$r(mu,s)))
#'## plot contours using a linear scale for the parameters
#'## draw four contours corresponding to the following likelihood ratios:
#'##  0.5, 0.1, Chi2 with 2 df and p values of 0.95 and 0.99
#'X11(width=12,height=6)
#'layout(matrix(1:2,ncol=2))
#'contour(MeanLog,SdLog,t(sampLNmleLNcontour),
#'        levels=c(log(c(0.5,0.1)),-0.5*qchisq(c(0.95,0.99),df=2)),
#'        labels=c("log(0.5)",
#'          "log(0.1)",
#'          "-1/2*P(Chi2=0.95)",
#'          "-1/2*P(Chi2=0.99)"),
#'        xlab=expression(mu),ylab=expression(sigma),
#'        main="Log Relative Likelihood Contours"
#'        )
#'points(sampLNmleLN$estimate[1],sampLNmleLN$estimate[2],pch=3)
#'points(meanlog.true,sdlog.true,pch=16,col=2)
#'## The contours are not really symmetrical about the MLE we can try to
#'## replot them using a log scale for the parameters to see if that improves
#'## the situation
#'contour(MeanLog,log(SdLog),t(sampLNmleLNcontour),
#'        levels=c(log(c(0.5,0.1)),-0.5*qchisq(c(0.95,0.99),df=2)),
#'        labels="",
#'        xlab=expression(mu),ylab=expression(log(sigma)),
#'        main="Log Relative Likelihood Contours",
#'        sub=expression(paste("log scale for parameter: ",sigma)))
#'points(sampLNmleLN$estimate[1],log(sampLNmleLN$estimate[2]),pch=3)
#'points(meanlog.true,log(sdlog.true),pch=16,col=2)
#'
#'## make a parametric boostrap to check the distribution of the deviance
#'nbReplicate <- 10000
#'sampleSize <- 100
#'system.time(
#'            devianceLN100 <- replicate(nbReplicate,{
#'                             sampLN <- rlnorm(sampleSize,meanlog=meanlog.true,sdlog=sdlog.true)
#'                             sampLNmleLN <- lnormMLE(sampLN)
#'                             -2*sampLNmleLN$r(meanlog.true,sdlog.true)
#'                           }
#'                                       )
#'            )[3]
#'
#'## Get 95 and 99% confidence intervals for the QQ plot
#'ci <- sapply(1:nbReplicate,
#'                 function(idx) qchisq(qbeta(c(0.005,0.025,0.975,0.995),
#'                                            idx,
#'                                            nbReplicate-idx+1),
#'                                      df=2)
#'             )
#'## make QQ plot
#'X <- qchisq(ppoints(nbReplicate),df=2)
#'Y <- sort(devianceLN100)
#'X11()
#'plot(X,Y,type="n",
#'     xlab=expression(paste(chi[2]^2," quantiles")),
#'     ylab="MC quantiles",
#'     main="Deviance with true parameters after ML fit of logNorm data",
#'     sub=paste("sample size:", sampleSize,"MC replicates:", nbReplicate)
#'     )
#'abline(a=0,b=1)
#'lines(X,ci[1,],lty=2)
#'lines(X,ci[2,],lty=2)
#'lines(X,ci[3,],lty=2)
#'lines(X,ci[4,],lty=2)
#'lines(X,Y,col=2)
#'}
#'
lnormMLE <- function(yi,
                     ni = numeric(length(yi))+1,
                     si = numeric(length(yi))+1
                     ) {

  ## check if yi is a spikeTrain object if yes take the "diff"
  if (inherits(yi,"spikeTrain")) yi <- diff(yi)
  ## check that yi elements are positive
  if (any(yi <= 0)) stop("yi elements must be strictly positive")
  ## coerce yi to vector
  yi <- as.numeric(yi)
  ## check that ni has the same length as yi
  if (!identical(length(yi),length(ni))) stop("yi and ni should have the same length")
  ## check that the elements of ni are non negative and integer
  if (any(ni < 0)) stop("ni elements must be non-negative")
  if (!identical(ni, round(ni))) stop("ni should be a vector of positive integers")
  if (!identical(length(yi),length(si))) stop("yi and si should have the same length")
  ## check that the elements of ni are non negative and integer
  if (any(si < 0)) stop("si elements must be non-negative")
  if (!identical(si, round(si))) stop("si should be a vector of positive integers")
  if (any(si > ni)) stop("si elements should not be greater than ni elements")

  ## Get the number of uncensored events
  s.dot <- sum(si)
  ## Get the total number of events
  n.dot <- sum(ni)
  ci <- ni-si
  c.dot <- sum(ci)
  ## Create function returning the opposite of the log likelihood
  minusLogLik <- function(p) {
    if (missing(p)) {
      txt <- paste("This function argument should be a 2 component vector:\n",
                   "  component 1 is the meanlog parameter,\n",
                   "  component 2 is the log of the sdlog parameter.\n"
                   )
      cat(txt)
    } else {
      mean <- p[1]
      sd <- exp(p[2])
      -(ifelse(s.dot>0,sum(dlnorm(yi[si>0],mean,sd,log=TRUE)*si[si>0],0))+
        ifelse(c.dot>0,sum(plnorm(yi[ci>0],mean,sd,lower.tail=FALSE,log.p=TRUE)*ci[ci>0]),0)
        )
    }
  }
  
  ## Get the empirical mean which is the MLE of parameter mu
  if (s.dot >= 2) {
    mu.hat <- weighted.mean(log(yi), si)
    s2 <- weighted.mean(log(yi)^2, si) - mu.hat^2
  } else {
    mu.hat <- weighted.mean(log(yi), ni)
    s2 <- weighted.mean(log(yi)^2, ni) - mu.hat^2
  }
  if (s.dot == n.dot) {
    ## all events are uncensored
    estimate <- c(mu.hat,sqrt(s2))
    se <- c(sqrt(s2/n.dot),sqrt(s2/(2*n.dot)))
    l <- -minusLogLik(c(estimate[1],log(estimate[2])))
  } else {
    ## numerical fit needed
    ## mleFit <- nlm(minusLogLik,c(mu.hat,0.5*log(s2)),hessian=TRUE)
    mleFit <- optim(fn=minusLogLik,
                    par=c(mu.hat,0.5*log(s2)),
                    method="BFGS",
                    hessian=TRUE)
    ## estimate <- c(mleFit$estimate[1],exp(mleFit$estimate[2]))
    estimate <- c(mleFit$par[1],exp(mleFit$par[2]))
    newVar <- c(1,1/estimate[2]) %o% c(1,1/estimate[2])
    se <- sqrt(diag(solve(mleFit$hessian * newVar)))
    ## l <- -mleFit$minimum
    l <- -mleFit$value
  }
  names(estimate) <- names(formals(dlnorm))[2:3]
  names(se) <- names(formals(dlnorm))[2:3]
  result <- list(estimate = estimate,
                 se = se,
                 logLik = l,
                 r = function(meanlog,sdlog) -minusLogLik(c(meanlog,log(sdlog))) - l,
                 mll = minusLogLik,
                 call = match.call()
                 )
  class(result) <- "durationFit"
  return(result)

}



#'Maximum Likelihood Parameter Estimation of a Refractory Exponential Model
#'with Possibly Censored Data
#'
#'Estimate refractory exponential model parameters by the maximum likelihood
#'method using possibly censored data.
#'
#'The MLE are available in closed form even in the censored case for this
#'model. The likelihood function cannot be differentiated with respect to the
#'\code{rp} (refractory period) parameter at the maximum. COnfidence intervals
#'for this parameter are therefore not available.
#'
#'@param yi vector of (possibly binned) observations or a \code{spikeTrain}
#'object.
#'@param ni vector of counts for each value of \code{yi}; default:
#'\code{numeric(length(yi))+1}.
#'@param si vector of counts of \emph{uncensored} observations for each value
#'of \code{yi}; default: \code{numeric(length(yi))+1}.
#'@return A list of class \code{durationFit} with the following components:
#'@returnItem estimate the estimated parameters, a named vector.
#'@returnItem se the standard errors, a named vector.
#'@returnItem logLik the log likelihood at maximum.
#'@returnItem r a function returning the log of the relative likelihood
#'function.
#'@returnItem mll a function returning the opposite of the log likelihood
#'function using the log of the parameters.
#'@returnItem call the matched call.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{drexp}}, \code{\link{invgaussMLE}},
#'\code{\link{lnormMLE}}, \code{\link{gammaMLE}}, \code{\link{weibullMLE}}
#'@keywords distribution ts
#'@examples
#'
#'\dontrun{
#'## Simulate sample of size 100 from a refractory exponential distribution
#'set.seed(1102006,"Mersenne-Twister")
#'sampleSize <- 100
#'rate.true <- 20
#'rp.true <- 0.01
#'sampRE <- rrexp(sampleSize,rate=rate.true,rp=rp.true)
#'sampREmleRE <- rexpMLE(sampRE)
#'rbind(est = sampREmleRE$estimate,se = sampREmleRE$se,true = c(rate.true,rp.true))
#'
#'## make a parametric boostrap to check the distribution of the deviance
#'nbReplicate <- 10000
#'system.time(
#'            devianceRE100 <- replicate(nbReplicate,{
#'              sampRE <- rrexp(sampleSize,rate=rate.true,rp=rp.true)
#'              sampREmleRE <- rexpMLE(sampRE)
#'              -2*sampREmleRE$r(rate.true,rp.true)
#'            }
#'                                       )
#'            )[3]
#'
#'## Get 95 and 99% confidence intervals for the QQ plot
#'ci <- sapply(1:nbReplicate,
#'                 function(idx) qchisq(qbeta(c(0.005,0.025,0.975,0.995),
#'                                            idx,
#'                                            nbReplicate-idx+1),
#'                                      df=2)
#'             )
#'## make QQ plot
#'X <- qchisq(ppoints(nbReplicate),df=2)
#'Y <- sort(devianceRE100)
#'X11()
#'plot(X,Y,type="n",
#'     xlab=expression(paste(chi[2]^2," quantiles")),
#'     ylab="MC quantiles",
#'     main="Deviance with true parameters after ML fit of refractory Poisson data",
#'     sub=paste("sample size:", sampleSize,"MC replicates:", nbReplicate)
#'     )
#'abline(a=0,b=1)
#'lines(X,ci[1,],lty=2)
#'lines(X,ci[2,],lty=2)
#'lines(X,ci[3,],lty=2)
#'lines(X,ci[4,],lty=2)
#'lines(X,Y,col=2)
#'}
#'
rexpMLE <- function(yi,
                    ni = numeric(length(yi))+1,
                    si = numeric(length(yi))+1
                    ) {

  ## check if yi is a spikeTrain object if yes take the "diff"
  if (inherits(yi,"spikeTrain")) yi <- diff(yi)
  ## check that yi elements are positive
  if (any(yi <= 0)) stop("yi elements must be positive or null")
  ## coerce yi to vector
  yi <- as.numeric(yi)
  ## check that ni has the same length as yi
  if (!identical(length(yi),length(ni))) stop("yi and ni should have the same length")
  ## check that the elements of ni are non negative and integer
  if (any(ni < 0)) stop("ni elements must be non-negative")
  if (!identical(ni, round(ni))) stop("ni should be a vector of positive integers")
  if (!identical(length(yi),length(si))) stop("yi and si should have the same length")
  ## check that the elements of ni are non negative and integer
  if (any(si < 0)) stop("si elements must be non-negative")
  if (!identical(si, round(si))) stop("si should be a vector of positive integers")
  if (any(si > ni)) stop("si elements should not be greater than ni elements")

  ## Create function returning the opposite of the log likelihood
  minusLogLik <- function(p) {
    if (missing(p)) {
      txt <- paste("This function argument should be a 2 component vector:\n",
                   "  component 1 is the log of the rate parameter,\n",
                   "  component 2 is the log of the rp (refractory period) parameter.\n"
                   )
      cat(txt)
    } else {
      rate <- exp(p[1])
      rp <- exp(p[2])
      -(ifelse(s.dot>0,sum(drexp(yi[si>0],rate=rate,rp=rp,log=TRUE)*si[si>0],0))+
        ifelse(c.dot>0,sum(prexp(yi[ci>0],rate=rate,
                                    rp=rp,lower.tail=FALSE,log.p=TRUE)*ci[ci>0]),0)
        )
    }
  }
  
  ## Get the number of uncensored events
  s.dot <- sum(si)
  ## Get the total number of events
  n.dot <- sum(ni)
  ci <- ni-si
  c.dot <- sum(ci)
  if (s.dot == 0) stop("No uncensored events")
  rp.hat <- min(yi[si==1]) - .Machine$double.eps
  mu.hat <- weighted.mean(yi-rp.hat, si)
  rate.hat <- 1/mu.hat
  estimate <- c(rate.hat,rp.hat)
  se <- c(rate.hat/sqrt(s.dot),NA)
  l <- -minusLogLik(log(c(rate.hat,rp.hat)))

  names(estimate) <- c("rate","rp")
  names(se) <- c("rate","rp")
  rFct <- function(rate,rp) -minusLogLik(log(c(rate,rp))) - l
  
  result <- list(estimate = estimate,
                 se = se,
                 logLik = l,
                 r = rFct,
                 mll = minusLogLik,
                 call = match.call()
                 )
  class(result) <- "durationFit"
  return(result)

}




#'Maximum Likelihood Parameter Estimation of a Weibull Model with Possibly
#'Censored Data
#'
#'Estimate Weibull model parameters by the maximum likelihood method using
#'possibly censored data.
#'
#'There is no closed form expression for the MLE of a Weibull distribution. The
#'numerical method implemented here uses the profile likelihood described by
#'Kalbfleisch (1985) pp 56-58.
#'
#'In order to ensure good behavior of the numerical optimization routines,
#'optimization is performed on the log of the parameters (\code{shape} and
#'\code{scale}).
#'
#'Standard errors are obtained from the inverse of the observed information
#'matrix at the MLE. They are transformed to go from the log scale used by the
#'optimization routine to the parameterization requested.
#'
#'@param yi vector of (possibly binned) observations or a \code{spikeTrain}
#'object.
#'@param ni vector of counts for each value of \code{yi}; default:
#'\code{numeric(length(yi))+1}.
#'@param si vector of counts of \emph{uncensored} observations for each value
#'of \code{yi}; default: \code{numeric(length(yi))+1}.
#'@param shape.min numeric, the inital guess of the minimal possible value of
#'the \code{shape} parameter, used by \code{optimise}.
#'@param shape.max numeric, the inital guess of the maximal possible value of
#'the \code{shape} parameter, used by \code{optimise}.
#'@return A list of class \code{durationFit} with the following components:
#'@returnItem estimate the estimated parameters, a named vector.
#'@returnItem se the standard errors, a named vector.
#'@returnItem logLik the log likelihood at maximum.
#'@returnItem r a function returning the log of the relative likelihood
#'function.
#'@returnItem mll a function returning the opposite of the log likelihood
#'function using the log of the parameters.
#'@returnItem call the matched call.
#'@note The returned standard errors (component \code{se}) are valid in the
#'asymptotic limit. You should plot contours using function \code{r} in the
#'returned list and check that the contours are reasonably close to ellipses.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{Weibull}}, \code{\link{invgaussMLE}},
#'\code{\link{lnormMLE}}, \code{\link{gammaMLE}}
#'@references Kalbfleisch, J. G. (1985) \emph{Probability and Statistical
#'Inference. Volume 2: Statistical Inference}. Springer-Verlag.
#'
#'Lindsey, J.K. (2004) \emph{Introduction to Applied Statistics: A Modelling
#'Approach}. OUP.
#'@keywords distribution ts
#'@examples
#'
#'\dontrun{
#'## Simulate sample of size 100 from a weibull distribution
#'set.seed(1102006,"Mersenne-Twister")
#'sampleSize <- 100
#'shape.true <- 2.5
#'scale.true <- 0.085
#'sampWB <- rweibull(sampleSize,shape=shape.true,scale=scale.true)
#'sampWBmleWB <- weibullMLE(sampWB)
#'rbind(est = sampWBmleWB$estimate,se = sampWBmleWB$se,true = c(shape.true,scale.true))
#'
#'## Estimate the log relative likelihood on a grid to plot contours
#'Shape <- seq(sampWBmleWB$estimate[1]-4*sampWBmleWB$se[1],
#'               sampWBmleWB$estimate[1]+4*sampWBmleWB$se[1],
#'               sampWBmleWB$se[1]/10)
#'Scale <- seq(sampWBmleWB$estimate[2]-4*sampWBmleWB$se[2],
#'             sampWBmleWB$estimate[2]+4*sampWBmleWB$se[2],
#'             sampWBmleWB$se[2]/10)
#'sampWBmleWBcontour <- sapply(Shape, function(sh) sapply(Scale, function(sc) sampWBmleWB$r(sh,sc)))
#'## plot contours using a linear scale for the parameters
#'## draw four contours corresponding to the following likelihood ratios:
#'##  0.5, 0.1, Chi2 with 2 df and p values of 0.95 and 0.99
#'X11(width=12,height=6)
#'layout(matrix(1:2,ncol=2))
#'contour(Shape,Scale,t(sampWBmleWBcontour),
#'        levels=c(log(c(0.5,0.1)),-0.5*qchisq(c(0.95,0.99),df=2)),
#'        labels=c("log(0.5)",
#'          "log(0.1)",
#'          "-1/2*P(Chi2=0.95)",
#'          "-1/2*P(Chi2=0.99)"),
#'        xlab="shape",ylab="scale",
#'        main="Log Relative Likelihood Contours"
#'        )
#'points(sampWBmleWB$estimate[1],sampWBmleWB$estimate[2],pch=3)
#'points(shape.true,scale.true,pch=16,col=2)
#'## The contours are not really symmetrical about the MLE we can try to
#'## replot them using a log scale for the parameters to see if that improves
#'## the situation
#'contour(log(Shape),log(Scale),t(sampWBmleWBcontour),
#'        levels=c(log(c(0.5,0.1)),-0.5*qchisq(c(0.95,0.99),df=2)),
#'        labels="",
#'        xlab="log(shape)",ylab="log(scale)",
#'        main="Log Relative Likelihood Contours",
#'        sub="log scale for the parameters")
#'points(log(sampWBmleWB$estimate[1]),log(sampWBmleWB$estimate[2]),pch=3)
#'points(log(shape.true),log(scale.true),pch=16,col=2)
#'
#'## make a parametric boostrap to check the distribution of the deviance
#'nbReplicate <- 10000
#'sampleSize <- 100
#'system.time(
#'            devianceWB100 <- replicate(nbReplicate,{
#'              sampWB <- rweibull(sampleSize,shape=shape.true,scale=scale.true)
#'              sampWBmleWB <- weibullMLE(sampWB)
#'              -2*sampWBmleWB$r(shape.true,scale.true)
#'            }
#'                                       )
#'            )[3]
#'
#'## Get 95 and 99% confidence intervals for the QQ plot
#'ci <- sapply(1:nbReplicate,
#'                 function(idx) qchisq(qbeta(c(0.005,0.025,0.975,0.995),
#'                                            idx,
#'                                            nbReplicate-idx+1),
#'                                      df=2)
#'             )
#'## make QQ plot
#'X <- qchisq(ppoints(nbReplicate),df=2)
#'Y <- sort(devianceWB100)
#'X11()
#'plot(X,Y,type="n",
#'     xlab=expression(paste(chi[2]^2," quantiles")),
#'     ylab="MC quantiles",
#'     main="Deviance with true parameters after ML fit of gamma data",
#'     sub=paste("sample size:", sampleSize,"MC replicates:", nbReplicate)
#'     )
#'abline(a=0,b=1)
#'lines(X,ci[1,],lty=2)
#'lines(X,ci[2,],lty=2)
#'lines(X,ci[3,],lty=2)
#'lines(X,ci[4,],lty=2)
#'lines(X,Y,col=2)
#'}
#'
weibullMLE <- function(yi,
                       ni = numeric(length(yi))+1,
                       si = numeric(length(yi))+1,
                       shape.min = 0.05,
                       shape.max = 5) {
                       
  #####################################################
  ## define internal function makeProfileWeibullLogLik
  #####################################################
  makeProfileWeibullLogLik <- function(data) {

    ## check that data is strictly positive
    if ( any(data <= 0) )
      stop("data should be a vector of positive numbers")
    ## if data is somthing else than a vector (eg a matrix)
    ## coerce it to a vector
    dim(data) <- NULL
    ## get the sample size
    nDot <- length(data)
    ## get the sum of the logs of the variate
    sumLog <- sum(log(data))
    
    logLik <- function(shape) {
      ## check that shape is strictly positive
      if (shape <= 0) {
        return(-Inf)
      } else {
        term1 <- nDot * log(shape)
        term2 <- (shape - 1) * sumLog
        term3 <- - nDot * (1 + log(sum(data^shape)/nDot))
        return( term1 + term2 + term3)
      }
    }

    return(logLik)
  
  }
  ###############################################################
  ## End of internal function makeProfileWeibullLogLik definition
  ###############################################################

  ## check if yi is a spikeTrain object if yes take the "diff"
  if (inherits(yi,"spikeTrain")) yi <- diff(yi)
  ## check that yi elements are positive
  if (any(yi <= 0)) stop("yi elements must be strictly positive")
  ## coerce yi to vector
  yi <- as.numeric(yi)
  ## check that ni has the same length as yi
  if (!identical(length(yi),length(ni))) stop("yi and ni should have the same length")
  ## check that the elements of ni are non negative and integer
  if (any(ni < 0)) stop("ni elements must be non-negative")
  if (!identical(ni, round(ni))) stop("ni should be a vector of positive integers")
  if (!identical(length(yi),length(si))) stop("yi and si should have the same length")
  ## check that the elements of ni are non negative and integer
  if (any(si < 0)) stop("si elements must be non-negative")
  if (!identical(si, round(si))) stop("si should be a vector of positive integers")
  if (any(si > ni)) stop("si elements should not be greater than ni elements")

  ## Create function returning the opposite of the log likelihood
  minusLogLik <- function(p) {
    if (missing(p)) {
      txt <- paste("This function argument should be a 2 component vector:\n",
                   "  component 1 is the log of the shape parameter,\n",
                   "  component 2 is the log of the scale parameter.\n"
                   )
      cat(txt)
    } else {
      shape <- exp(p[1])
      scale <- exp(p[2])
      -(ifelse(s.dot>0,sum(dweibull(yi[si>0],shape=shape,scale=scale,log=TRUE)*si[si>0],0))+
        ifelse(c.dot>0,sum(pweibull(yi[ci>0],shape=shape,
                                    scale=scale,lower.tail=FALSE,log.p=TRUE)*ci[ci>0]),0)
        )
    }
  }
  
  ## Get the number of uncensored events
  s.dot <- sum(si)
  ## Get the total number of events
  n.dot <- sum(ni)
  ci <- ni-si
  c.dot <- sum(ci)

  if (s.dot == n.dot) {
    ## no censored event
    mu.hat <- weighted.mean(yi, si)
    s2 <- weighted.mean(yi^2, si) - mu.hat^2
    stat1 <- mu.hat^2/(s2 + mu.hat^2)
    shapeFct <- function(shape) (gamma(1 + 1/shape))^2 / gamma(1 + 2/shape) - stat1
    myRoot <- uniroot(shapeFct, lower = shape.min, upper = shape.max)$root
    shapeInterval <- c(floor(myRoot) - .Machine$double.eps,ceiling(myRoot) + .Machine$double.eps)
    data <- rep(yi[si>0],each=ni[si>0])
    myLogLik <- makeProfileWeibullLogLik(data)
    shape.hat <- optimize(myLogLik, interval = shapeInterval, maximum = TRUE)$maximum
    scale.hat <- (sum(data^shape.hat)/s.dot)^(1/shape.hat)
  } else {
    ## some censored events
    ## if more than 10 events are uncesored get inital guess from them
    ## otherwise use all events
    if (s.dot >= 10) {
      mu.hat <- weighted.mean(yi, si)
      s2 <- weighted.mean(yi^2, si) - mu.hat^2
      stat1 <- mu.hat^2/(s2 + mu.hat^2)
      shapeFct <- function(shape) (gamma(1 + 1/shape))^2 / gamma(1 + 2/shape) - stat1
      myRoot <- uniroot(shapeFct, lower = shape.min, upper = shape.max)$root
      shapeInterval <- c(floor(myRoot) - .Machine$double.eps,ceiling(myRoot) + .Machine$double.eps)
      data <- rep(yi[si>0],each=ni[si>0])
      myLogLik <- makeProfileWeibullLogLik(data)
      shape.hat <- optimize(myLogLik, interval = shapeInterval, maximum = TRUE)$maximum
      scale.hat <- (sum(data^shape.hat)/s.dot)^(1/shape.hat)
    } else {
      mu.hat <- weighted.mean(yi, ni)
      s2 <- weighted.mean(yi^2, ni) - mu.hat^2
      stat1 <- mu.hat^2/(s2 + mu.hat^2)
      shapeFct <- function(shape) (gamma(1 + 1/shape))^2 / gamma(1 + 2/shape) - stat1
      myRoot <- uniroot(shapeFct, lower = shape.min, upper = shape.max)$root
      shapeInterval <- c(floor(myRoot) - .Machine$double.eps,ceiling(myRoot) + .Machine$double.eps)
      data <- rep(yi,each=ni)
      myLogLik <- makeProfileWeibullLogLik(data)
      shape.hat <- optimize(myLogLik, interval = shapeInterval, maximum = TRUE)$maximum
      scale.hat <- (sum(data^shape.hat)/n.dot)^(1/shape.hat)
    }
  } ## End of conditional on s.dot == n.dot
  ## mleFit <- nlm(minusLogLik,log(c(shape.hat,scale.hat)),hessian=TRUE)
  mleFit <- optim(fn=minusLogLik,
                  par=log(c(shape.hat,scale.hat)),
                  method="BFGS",
                  hessian=TRUE)
  ## estimate <- exp(mleFit$estimate)
  estimate <- exp(mleFit$par)
  newVar <- (1/estimate) %o% (1/estimate)
  observedI <- mleFit$hessian * newVar
  se <- sqrt(diag(solve(observedI)))
  ## l <- -mleFit$minimum
  l <- -mleFit$value

  names(estimate) <- c("shape","scale")
  names(se) <- c("shape","scale")
  rFct <- function(shape,scale) -minusLogLik(log(c(shape,scale))) - l
  
  result <- list(estimate = estimate,
                 se = se,
                 logLik = l,
                 r = rFct,
                 mll = minusLogLik,
                 call = match.call()
                 )
  class(result) <- "durationFit"
  return(result)

}

is.durationFit <- function(obj) {

  if (!("durationFit" %in% class(obj))) return(FALSE)
  
  expectedNames <- c("estimate",
                     "se",
                     "logLik",
                     "r",
                     "mll",
                     "call")
  if (!all(expectedNames %in% names(obj))) return(FALSE)
  TRUE
    
}



#'Quantile-Quantile Plot For Fitted Duration Distributions
#'
#'Produces a QQ plot of empirical against theoretical quantiles of one of the
#'following duration distributions: inverse Gaussian, log normal, log logistic,
#'refractory exponential, gamma, weibull.
#'
#'If the data to which the model was fitted have censored events, the latter
#'are not used to build the empirical quantiles.
#'
#'@param durationFit a \code{durationFit} object, that is, a list returned by
#'one of these functions: \code{\link{invgaussMLE}}, \code{\link{lnormMLE}},
#'\code{\link{llogisMLE}}, \code{\link{rexpMLE}}, \code{\link{gammaMLE}},
#'\code{\link{weibullMLE}}.
#'@param CI a numeric vector with at most tow components, the confidence
#'intervals to be drawn. If \code{NULL}, intervals are not drawn.
#'@param type,xlab,ylab,main,sub,ylim see \code{\link{plot}}, default values
#'are provided if arguments are missing.
#'@param dataLwd non negative integer, the width of the line used to draw the
#'data.
#'@param ablineCol color of the diagonal.
#'@param \dots additional arguments passed to \code{\link{plot}}.
#'@return Nothing is returned, the function is used for its side effect, a plot
#'is generated.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{compModels}}, \code{\link{invgaussMLE}},
#'\code{\link{lnormMLE}}, \code{\link{llogisMLE}}, \code{\link{rexpMLE}},
#'\code{\link{gammaMLE}}, \code{\link{weibullMLE}}
#'@keywords distribution ts
#'@examples
#'
#'\dontrun{
#'## Simulate a sample with 100 events from an inverse Gaussian
#'set.seed(1102006,"Mersenne-Twister")
#'mu.true <- 0.075
#'sigma2.true <- 3
#'sampleSize <- 100
#'sampIG <- rinvgauss(sampleSize,mu=mu.true,sigma2=sigma2.true)
#'## Fit it with an inverse Gaussian Model
#'sampIGmleIG <- invgaussMLE(sampIG)
#'## draw the QQ plot on a log scale
#'qqDuration(sampIGmleIG,log="xy")
#'## Fit it with a log normal Model
#'sampIGmleLN <- lnormMLE(sampIG)
#'## draw the QQ plot on a log scale
#'qqDuration(sampIGmleLN,log="xy")
#'## Fit it with a gamma Model
#'sampIGmleGA <- gammaMLE(sampIG)
#'## draw the QQ plot on a log scale
#'qqDuration(sampIGmleGA,log="xy")
#'## Fit it with a Weibull Model
#'sampIGmleWB <- weibullMLE(sampIG)
#'## draw the QQ plot on a log scale
#'qqDuration(sampIGmleWB,log="xy")
#'## Fit it with a refractory exponential Model
#'sampIGmleRE <- rexpMLE(sampIG)
#'## draw the QQ plot on a log scale
#'qqDuration(sampIGmleRE,log="xy")
#'## Fit it with a log logisitc Model
#'sampIGmleLL <- llogisMLE(sampIG)
#'## draw the QQ plot on a log scale
#'qqDuration(sampIGmleLL,log="xy")
#'}
#'
qqDuration <- function(durationFit,
                       CI=c(0.95,0.99),
                       type="l",
                       xlab,
                       ylab,
                       main,
                       sub,
                       ylim,
                       dataLwd=2,
                       ablineCol=2,
                       ...
                       ) {
  ## Check that durationFit is a durationFit object
  if (!is.durationFit(durationFit))
    stop(paste(deparse(substitute(durationFit)),
               "should be a durationFit object."
               )
         )

  ## check CI
  if (!is.null(CI)) {
    ## make sure that CI is at most of length 2 otherwise
    ## keep the first 2 components
    if (length(CI) > 2) CI <- CI[1:2]
    ## Check that each component of CI is between 0 and 1
    if (any(CI>=1 | CI<=0))
      stop(paste(deparse(substitute(CI)),
                 "components should be in (0,1)")
           )
  } ## End of conditional on !is.null(CI)
  
  ## get the uncensored events
  ## The next code line causes an unimportant warning during
  ## automatic code check by "R CMD check"
  x <- evalq(rep(yi[si>0],each=si[si>0]),envir=environment(durationFit$r))
  ## get the y values of the QQ plot
  y <- sort(x)
  ## get the theoretical x values
  distribution <- switch(deparse(durationFit$call[[1]]),
                         invgaussMLE = qinvgauss,
                         gammaMLE = qgamma,
                         weibullMLE = qweibull,
                         lnormMLE = qlnorm,
                         rexpMLE = qrexp,
                         llogisMLE = qllogis)

  if (is.null(distribution))
    stop(paste("Unknown distribution in",
               deparse(substitute(durationFit))
               )
         )

  sampleSize <- length(y)
  pp <- ppoints(sampleSize)
  x.theo <- do.call(distribution,
                    c(list(p=pp),
                      as.list(durationFit$estimate)
                      )
                    )


  if (!is.null(CI)) {
    ## Get the confidence intervals
    ci <- sapply(1:sampleSize,
                 function(idx)
                 as.vector(sapply(CI,
                                  function(l) {
                                    thePs <- c((1-l)/2,1-(1-l)/2)
                                    thePs <- qbeta(thePs,
                                                   idx,
                                                   sampleSize-idx+1)
                                    do.call(distribution,c(list(p=thePs),
                                                           as.list(durationFit$estimate))
                                            )
                                    
                                  }
                                  )
                           )
                 )
  } ## End of conditional on !is.null(CI)
  
  modelName <- switch(deparse(durationFit$call[[1]]),
                      invgaussMLE = "invgauss",
                      gammaMLE = "gamma",
                      weibullMLE = "weibull",
                      lnormMLE = "lnorm",
                      rexpMLE = "rexp",
                      llogisMLE = "llogis")
  sampleName <- deparse(durationFit$call[["yi"]])

  if (missing(xlab)) xlab <- paste("Quantiles of",
                                   modelName)
  if (missing(ylab)) ylab <- "Sample quantiles"
  if (missing(main)) main <- paste("QQ plot of",
                                   sampleName,
                                   "vs fitted",
                                   modelName)
  if (missing(sub)) {
    sub <- paste(sampleSize,
                 "uncensored intervals used.")
    if (!is.null(CI))
      sub <- paste(sub,
                   "CI at:",
                   paste(CI,collapse=",")
                   )
  }

  if (missing(ylim)) ylim <- c(min(min(x.theo),min(y)),
                               max(max(x.theo),max(y))
                               )
  
  plot(x.theo,y,
       xlab=xlab,
       ylab=ylab,
       main=main,
       sub=sub,
       ylim=ylim,
       type="n",
       ...
       )
  abline(a=0,b=1,col=ablineCol)
  if (!is.null(CI)) apply(ci,1,function(x) lines(x.theo,x,lty=2))
  lines(x.theo,y,type=type,lwd=dataLwd)
  
}



#'Compare Duration Models on a Specific Data Set
#'
#'Fit duration models with the maximum likelihood method to a given duration
#'data set. The data can be censored. The models should be among the following
#'list: inverse Gaussian, log normal, log logistic, gamma, Weibull, refractory
#'exponential. The Akaike information criterion (AIC) is used to produce a
#'numerical output. Diagnostic QQ or survival plots can also be generated.
#'
#'Fits are performed by maximizing the likelihood.
#'
#'@param yi vector of (possibly binned) observations or a \code{spikeTrain}
#'object.
#'@param ni vector of counts for each value of \code{yi}; default:
#'\code{numeric(length(yi))+1}.
#'@param si vector of counts of \emph{uncensored} observations for each value
#'of \code{yi}; default: \code{numeric(length(yi))+1}.
#'@param models a character vector whose elements are selected among:
#'\code{"invgauss"}, \code{"lnorm"}, \code{"gamma"}, \code{"weibull"},
#'\code{"llogis"}, \code{"rexp"}.
#'@param type should a QQ plot (\code{"qq"}) or a survival plot (\code{"s"}) be
#'generated?
#'@param log should a log scale be used?
#'@param plot should a plot be generated?
#'@return A vector whose component are nammed according to the model used and
#'ordered along increasing AIC values.
#'
#'if argument \code{plot} is set to \code{TRUE} (the default), a plot is
#'generated as a side effect.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{qqDuration}}, \code{\link{invgaussMLE}},
#'\code{\link{lnormMLE}}, \code{\link{llogisMLE}}, \code{\link{rexpMLE}},
#'\code{\link{gammaMLE}}, \code{\link{weibullMLE}}
#'@references Lindsey, J.K. (2004) \emph{The Statistical Analysis of Stochastic
#'Processes in Time}. CUP.
#'@keywords distribution ts
#'@examples
#'
#'\dontrun{
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
#'## It does not look too bad so let fit simple models
#'compModels(CAL1S[["neuron 1"]])
#'
#'## Simulate a sample with 100 events from an inverse Gaussian
#'set.seed(1102006,"Mersenne-Twister")
#'mu.true <- 0.075
#'sigma2.true <- 3
#'sampleSize <- 100
#'sampIG <- rinvgauss(sampleSize,mu=mu.true,sigma2=sigma2.true)
#'
#'## Compare models and display QQ plot
#'compModels(sampIG,type="qq")
#'
#'## Compare models and display survival plot
#'compModels(sampIG,type="s")
#'
#'
#'## Generate a censored sample using an exponential distribution
#'sampEXP <- rexp(sampleSize,1/(2*mu.true))
#'sampIGtime <- pmin(sampIG,sampEXP)
#'sampIGstatus <- as.numeric(sampIG <= sampEXP)
#'## Compare models and display QQ plot
#'## WARNING with censored data like here the QQ plot is misleading
#'compModels(yi=sampIGtime,si=sampIGstatus,type="qq")
#'## Compare models and display survival plot
#'compModels(yi=sampIGtime,si=sampIGstatus,type="s")
#'}
#'
compModels <- function(yi,
                       ni = numeric(length(yi))+1,
                       si = numeric(length(yi))+1,
                       models = c("invgauss","lnorm","gamma","weibull","llogis","rexp"),
                       type = c("qq","s"),
                       log = TRUE,
                       plot=TRUE
                       ) {

  if (is.spikeTrain(yi)) yi <- diff(yi)
  ## Check that the models are known
  knownModels <- c("invgauss",
                   "lnorm",
                   "gamma",
                   "weibull",
                   "llogis",
                   "rexp")

  if (!all(models %in% knownModels))
    stop("Some requested models are not implemented.")

  theFits <- lapply(models,
                    function(m) switch(m,
                                       invgauss = invgaussMLE(yi,ni,si),
                                       lnorm = lnormMLE(yi,ni,si),
                                       gamma = gammaMLE(yi,ni,si),
                                       weibull = weibullMLE(yi,ni,si),
                                       llogis = llogisMLE(yi,ni,si),
                                       rexp = rexpMLE(yi,ni,si)
                                       )
                    )

  if (plot) {
    myNrow <- (length(models)-1) %/% 2 + 1
    layout(matrix(1:(myNrow*2),nrow=myNrow))
    oldpar <- par(mar=c(4,3,2,1))
    on.exit(par(oldpar))
    if (type[1] == "qq") {
      sapply(1:length(models),
             function(idx) {
               if (log) {
                 qqDuration(theFits[[idx]],
                          main="",
                            ylab="",
                            sub="",
                            log="xy")
               } else {
                 qqDuration(theFits[[idx]],
                            main="",
                            ylab="",
                            sub="")
               }
             }
             )
      
    } else {
      
      isi <- rep(yi[ni>0],each=ni[ni>0])
      status <- unlist(lapply(1:length(ni),
                              function(idx) {
                                if (ni[idx]==0) return(numeric())
                                else {
                                  result <- numeric(ni[idx])
                                  if (si[idx]>0) result[1:si[idx]] <- 1
                                  return(result)
                                }
                              }
                              )
                       )
      kmFit <- survfit(Surv(isi,status) ~1)
      sapply(1:length(models),
             function(idx) {
               plot(kmFit,log=log,mark.time=FALSE,main=models[idx])
               survFct <- switch(models[idx],
                                 invgauss = pinvgauss,
                                 lnorm = plnorm,
                                 gamma = pgamma,
                                 weibull = pweibull,
                                 llogis = pllogis,
                                 rexp = prexp)
               isi <- sort(isi)
               y <- do.call(survFct,c(list(q=isi,lower.tail=FALSE),as.list(theFits[[idx]]$estimate)))
               lines(isi,y,col=2)
             }
             )
    }
  } ## End of conditional on plot
  aic <- sapply(theFits, function(f) -2*f$logLik+4)
  ordered <- sort.int(aic,index.return=TRUE)
  result <- ordered$x
  names(result) <- models[ordered$ix]
  result

}



#'The Inverse Gaussian Distribution
#'
#'Density, distribution function, quantile function, and random generation for
#'the inverse Gaussian.
#'
#'With the default, \code{"sigma2"}, parameterization (\code{mu = m, sigma2 =
#'s^2}) the inverse Gaussian distribution has density: \deqn{% }{% f(x) =
#'1/sqrt(2*pi*x^3*s^2) * exp(-0.5*(x-mu)^2/(x*s^2*m^2)) }\deqn{
#'f(x)=\frac{1}{\sqrt{2 \, \pi \, \sigma^2 \, x^3}} \, \exp% }{% f(x) =
#'1/sqrt(2*pi*x^3*s^2) * exp(-0.5*(x-mu)^2/(x*s^2*m^2)) }\deqn{
#'(-\frac{1}{2}\frac{(x-\mu)^2}{x \, \sigma^2 \, \mu^2}) }{% f(x) =
#'1/sqrt(2*pi*x^3*s^2) * exp(-0.5*(x-mu)^2/(x*s^2*m^2)) }\deqn{ }{% f(x) =
#'1/sqrt(2*pi*x^3*s^2) * exp(-0.5*(x-mu)^2/(x*s^2*m^2)) } with \eqn{\sigma^2 >
#'0}{s^2 > 0}.  The theoretical mean is: \eqn{\mu}{m} and the theoretical
#'variance is: \eqn{\mu^3 \sigma^2}{m^3*s^2}.  With the default,
#'\code{"boundary"}, parameterization (\code{mu = m, boundary = b})the inverse
#'Gaussian distribution has density: \deqn{% }{% f(x) = (b/sqrt(2*pi*x^3)) *
#'exp(-0.5*(x-b*mu)^2/(x*m^2)) }\deqn{ f(x)=\frac{b}{\sqrt{2 \, \pi \, x^3}} \,
#'\exp% }{% f(x) = (b/sqrt(2*pi*x^3)) * exp(-0.5*(x-b*mu)^2/(x*m^2)) }\deqn{
#'(-\frac{1}{2}\frac{(x-b \, \mu)^2}{x \, \mu^2}) }{% f(x) = (b/sqrt(2*pi*x^3))
#'* exp(-0.5*(x-b*mu)^2/(x*m^2)) }\deqn{ }{% f(x) = (b/sqrt(2*pi*x^3)) *
#'exp(-0.5*(x-b*mu)^2/(x*m^2)) } with \eqn{\sigma^2 > 0}{s^2 > 0}.  The
#'theoretical mean is: \eqn{\mu \, b}{m * b} and the theoretical variance is:
#'\eqn{\mu^3 \sigma^2}{m^3*b}.  The latent Brownian motion is described in
#'Lindsey (2004) pp 209-213, Whitemore and Seshadri (1987), Aalen and Gjessing
#'(2001) and Gerstein and Mandelbrot (1964).
#'
#'The expression for the distribution function is given in Eq. 4 of Whitemore
#'and Seshadri (1987).
#'
#'Initial guesses for the inversion of the distribution function used in
#'\code{qinvgauss} are obtained with the transformation of Whitemore and
#'Yalovsky (1978).
#'
#'Random variates are obtained with the method of Michael et al (1976) which is
#'also described by Devroye (1986, p 148) and Gentle (2003, p 193).
#'
#'@aliases dinvgauss qinvgauss pinvgauss rinvgauss
#'@param x,q vector of quantiles.
#'@param p vector of probabilities.
#'@param n number of observations. If \code{length(n) > 1}, the length is taken
#'to be the number required.
#'@param mu mean value of the distribution in the default parameterization,
#'\code{mean value / boundary} otherwise. Can also be viewed as the inverse of
#'the drift of the latent Brownian motion.
#'@param sigma2 variance of the latent Brownian motion. When this
#'parameterization is used (the default) the distance between the "starting"
#'point and the boundary ("absorbing barrier") is set to 1.
#'@param boundary distance between the starting point and the "absorbing
#'barrier" of the latent Brownian motion. When this parameterization is used
#'the Brownian motion variance is set to 1.
#'@param lower.tail logical; if \code{TRUE} (default), probabilities are
#'\code{P[X <= x]}, otherwise, \code{P[X > x]}.
#'@param log,log.p logical; if \code{TRUE}, probabilities p are given as
#'log(p).
#'@return \code{dinvgauss} gives the density, \code{pinvgauss} gives the
#'distribution function, \code{qinvgauss} gives the quantile function and
#'\code{rinvgauss} generates random deviates.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{invgaussMLE}}, \code{\link{Lognormal}},
#'\code{\link{hinvgauss}}
#'@references Gerstein, George L. and Mandelbrot, Benoit (1964) Random Walk
#'Models for the Spike Activity of a Single Neuron. \emph{Biophys J.} \bold{4}:
#'41--68.
#'\url{http://www.pubmedcentral.nih.gov/articlerender.fcgi?tool=pubmed&pubmedid=14104072}.
#'
#'Whitemore, G. A. and Yalovsky, M. (1978) A normalizing logarithmic
#'transformation for inverse Gaussian random variables. \emph{Technometrics}
#'\bold{20}: 207--208.
#'
#'Whitmore, G. A. and Seshadri, V. (1987) A Heuristic Derivation of the Inverse
#'Gaussian Distribution. \emph{The American Statistician} \bold{41}: 280--281.
#'
#'Aalen, Odd O. and Gjessing, Hakon K. (2001) Understanding the Shape of the
#'Hazard Rate: A Process Point of View. \emph{Statistical Science} \bold{16}:
#'1--14.
#'
#'Lindsey, J.K. (2004) \emph{Introduction to Applied Statistics: A Modelling
#'Approach}. OUP.
#'
#'Michael, J. R., Schucany, W. R. and Haas, R. W. (1976) Generating random
#'variates using transformations with multiple roots. \emph{The American
#'Statistician} \bold{30}: 88--90.
#'
#'Devroye, L. (1986) \emph{Non-Uniform Random Variate Generation}.
#'Springer-Verlag. \url{http://cg.scs.carleton.ca/~luc/rnbookindex.html}.
#'
#'Gentle, J. E. (2003) \emph{Random Number Generation and Monte Carlo Methods}.
#'Springer.
#'@keywords distribution ts
#'@examples
#'
#'\dontrun{
#'## Start with the inverse Gauss
#'## Define standard mu and sigma
#'mu.true <- 0.075 ## a mean ISI of 75 ms
#'sigma2.true <- 3
#'## Define a sequence of points on the time axis
#'X <- seq(0.001,0.3,0.001)
#'## look at the density
#'plot(X,dinvgauss(X,mu.true,sigma2.true),type="l",xlab="ISI (s)",ylab="Density")
#'
#'## Generate a sample of 100 ISI from this distribution
#'sampleSize <- 100
#'sampIG <- rinvgauss(sampleSize,mu=mu.true,sigma2=sigma2.true)
#'## check out the empirical survival function (obtained with the Kaplan-Meyer
#'## estimator) against the true one
#'library(survival)
#'sampIG.KMfit <- survfit(Surv(sampIG,1+numeric(length(sampIG))) ~1)
#'plot(sampIG.KMfit,log=TRUE)
#'lines(X,pinvgauss(X,mu.true,sigma2.true,lower.tail=FALSE),col=2)
#'
#'## Get a ML fit
#'sampIGmleIG <- invgaussMLE(sampIG)
#'## compare true and estimated parameters
#'rbind(est = sampIGmleIG$estimate,se = sampIGmleIG$se,true = c(mu.true,sigma2.true))
#'## plot contours of the log relative likelihood function
#'Mu <- seq(sampIGmleIG$estimate[1]-3*sampIGmleIG$se[1],
#'          sampIGmleIG$estimate[1]+3*sampIGmleIG$se[1],
#'          sampIGmleIG$se[1]/10)
#'Sigma2 <- seq(sampIGmleIG$estimate[2]-7*sampIGmleIG$se[2],
#'              sampIGmleIG$estimate[2]+7*sampIGmleIG$se[2],
#'              sampIGmleIG$se[2]/10)
#'sampIGmleIGcontour <- sapply(Mu, function(mu) sapply(Sigma2, function(s2) sampIGmleIG$r(mu,s2))) 
#'contour(Mu,Sigma2,t(sampIGmleIGcontour),
#'        levels=c(log(c(0.5,0.1)),-0.5*qchisq(c(0.95,0.99),df=2)),
#'        labels=c("log(0.5)",
#'          "log(0.1)",
#'          "-1/2*P(Chi2=0.95)",
#'          "-1/2*P(Chi2=0.99)"),
#'        xlab=expression(mu),ylab=expression(sigma^2))
#'points(mu.true,sigma2.true,pch=16,col=2)
#'## We can see that the contours are more parabola like on a log scale
#'contour(log(Mu),log(Sigma2),t(sampIGmleIGcontour),
#'        levels=c(log(c(0.5,0.1)),-0.5*qchisq(c(0.95,0.99),df=2)),
#'        labels=c("log(0.5)",
#'          "log(0.1)",
#'          "-1/2*P(Chi2=0.95)",
#'          "-1/2*P(Chi2=0.99)"),
#'        xlab=expression(log(mu)),ylab=expression(log(sigma^2)))
#'points(log(mu.true),log(sigma2.true),pch=16,col=2)
#'## make a deviance test for the true parameters
#'pchisq(-2*sampIGmleIG$r(mu.true,sigma2.true),df=2)
#'## check fit with a QQ plot
#'qqDuration(sampIGmleIG,log="xy")
#'
#'## Generate a censored sample using an exponential distribution
#'sampEXP <- rexp(sampleSize,1/(2*mu.true))
#'sampIGtime <- pmin(sampIG,sampEXP)
#'sampIGstatus <- as.numeric(sampIG <= sampEXP)
#'## fit the censored sample
#'sampIG2mleIG <- invgaussMLE(sampIGtime,,sampIGstatus)
#'## look at the results
#'rbind(est = sampIG2mleIG$estimate,
#'      se = sampIG2mleIG$se,
#'      true = c(mu.true,sigma2.true))
#'pchisq(-2*sampIG2mleIG$r(mu.true,sigma2.true),df=2)
#'## repeat the survival function estimation
#'sampIG2.KMfit <- survfit(Surv(sampIGtime,sampIGstatus) ~1)
#'plot(sampIG2.KMfit,log=TRUE)
#'lines(X,pinvgauss(X,sampIG2mleIG$estimate[1],sampIG2mleIG$estimate[2],lower.tail=FALSE),col=2)
#'}
#'
dinvgauss <- function(x,
                      mu = 1,
                      sigma2 = 1,
                      boundary = NULL,
                      log = FALSE
                      ){
  
  if( any(x <= 0) ) stop("y must contain positive values")
  if( any(mu <= 0) ) stop("mu must be positive")
  if (all(!is.null(sigma2)))
    if( any(sigma2 <= 0) ) stop("sigma2 must be positive")
  if (all(!is.null(boundary)))
    if( any(boundary <= 0) ) stop("boundary must be positive")
  if ( all(!is.null(sigma2)) && all(!is.null(boundary)) )
    stop("One of sigma2 or boundary must be specified, not both")

  if (all(!is.null(boundary))) {
    ## We work with the parameterization in term of boundary (ie, sigma2 = 1)
    ## We convert it in term of mu and sigma2
    sigma2 <- (1/boundary)^2
    mu <- boundary*mu
  }

  tmp <- -(x-mu)^2/(2*x*sigma2*mu^2)-(log(2*pi*sigma2)+3*log(x))/2
  if(!log) tmp <- exp(tmp)
  tmp

}



#'The Log Logistic Distribution
#'
#'Density, distribution function, quantile function, and random generation for
#'the log logistic.
#'
#'If \code{location} or \code{scale} are omitted, they assume the default
#'values of 0 and 1 respectively.
#'
#'The log-Logistic distribution with \code{location = m} and \code{scale = s}
#'has distribution function
#'
#'\deqn{\mathrm{F}(x) = \frac{1}{1+ \exp(-\frac{\log (x) - m}{s})}}{F(x) = 1 /
#'(1 + exp(-(log(x)-m)/s))}
#'
#'and density
#'
#'\deqn{\mathrm{f}(x)=\frac{1}{s \, x} \frac{\exp (-\frac{\log (x) -
#'m}{s})}{(1+ \exp(-\frac{\log (x) - m}{s}))^2}}{f(x) = 1/(s*x)
#'exp(-(log(x)-m)/s) (1 + exp(-(log(x)-m)/s))^-2.}
#'
#'@aliases dllogis pllogis qllogis rllogis
#'@param x,q vector of quantiles.
#'@param p vector of probabilities.
#'@param n number of observations. If \code{length(n) > 1}, the length is taken
#'to be the number required.
#'@param location,scale location and scale parameters (non-negative numeric).
#'@param lower.tail logical; if \code{TRUE} (default), probabilities are
#'\code{P[X <= x]}, otherwise, \code{P[X > x]}.
#'@param log,log.p logical; if \code{TRUE}, probabilities p are given as
#'log(p).
#'@return \code{dllogis} gives the density, \code{pllogis} gives the
#'distribution function, \code{qllogis} gives the quantile function and
#'\code{rllogis} generates random deviates.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{llogisMLE}}, \code{\link{Lognormal}},
#'\code{\link{hllogis}}
#'@references Lindsey, J.K. (2004) \emph{Introduction to Applied Statistics: A
#'Modelling Approach}. OUP.
#'
#'Lindsey, J.K. (2004) \emph{The Statistical Analysis of Stochastic Processes
#'in Time}. CUP.
#'@keywords distribution ts
#'@examples
#'
#'\dontrun{
#'tSeq <- seq(0.001,0.6,0.001)
#'location.true <- -2.7
#'scale.true <- 0.025
#'Yd <- dllogis(tSeq, location.true, scale.true)
#'Yh <- hllogis(tSeq, location.true, scale.true)
#'max.Yd <- max(Yd)
#'max.Yh <- max(Yh)
#'Yd <- Yd / max.Yd
#'Yh <- Yh / max.Yh
#'oldpar <- par(mar=c(5,4,4,4))
#'plot(tSeq, Yd, type="n", axes=FALSE, ann=FALSE,
#'     xlim=c(0,0.6), ylim=c(0,1))
#'axis(2,at=seq(0,1,0.2),labels=round(seq(0,1,0.2)*max.Yd,digits=2))
#'mtext("Density (1/s)", side=2, line=3)
#'axis(1,at=pretty(c(0,0.6)))
#'mtext("Time (s)", side=1, line=3)
#'axis(4, at=seq(0,1,0.2), labels=round(seq(0,1,0.2)*max.Yh,digits=2))
#'mtext("Hazard (1/s)", side=4, line=3, col=2)
#'mtext("Log Logistic Density and Hazard Functions", side=3, line=2,cex=1.5)
#'lines(tSeq,Yd)
#'lines(tSeq,Yh,col=2)
#'par(oldpar)
#'}
#'
dllogis <- function(x,
                    location = 0,
                    scale = 1,
                    log = FALSE) {
  
  ## Check that x elements are strictly positive
  if (any(x <= 0)) stop("x elements must be strictly positive")

  if (!log) dlogis(log(x),location,scale)/x
  else dlogis(log(x),location,scale,log=TRUE) - log(x)
  
}



#'The Refractory Exponential Distribution
#'
#'Density, distribution function, quantile function, and random generation for
#'the refractory exponential.
#'
#'The refractory exponential distribution with \code{rate}, r, and
#'\code{refractory period}, rp, has density:
#'
#'f(x) = r exp(- r (x-rp))
#'
#'for \code{x >= rp}.
#'
#'@aliases drexp prexp qrexp rrexp
#'@param x,q vector of quantiles.
#'@param p vector of probabilities.
#'@param n number of observations. If \code{length(n) > 1}, the length is taken
#'to be the number required.
#'@param lower.tail logical; if \code{TRUE} (default), probabilities are
#'\code{P[X <= x]}, otherwise, \code{P[X > x]}.
#'@param log,log.p logical; if \code{TRUE}, probabilities p are given as
#'log(p).
#'@param rate rate parameter (non-negative numeric).
#'@param rp refractory period parameter (non-negative numeric).
#'@return \code{drexp} gives the density, \code{prexp} gives the distribution
#'function, \code{qrexp} gives the quantile function and \code{rrexp} generates
#'random deviates.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{rexpMLE}}
#'@references Johnson, D. H. and Swami, A. (1983) The transmission of signals
#'by auditory-nerve fiber discharge patterns. \emph{J. Acoust. Soc. Am.}
#'\bold{74}: 493--501.
#'@keywords distribution ts
#'@examples
#'
#'\dontrun{
#'tSeq <- seq(0.001,0.6,0.001)
#'rate.true <- 20
#'rp.true <- 0.01
#'Yd <- drexp(tSeq, rate.true, rp.true)
#'Yh <- hrexp(tSeq, rate.true, rp.true)
#'max.Yd <- max(Yd)
#'max.Yh <- max(Yh)
#'Yd <- Yd / max.Yd
#'Yh <- Yh / max.Yh
#'oldpar <- par(mar=c(5,4,4,4))
#'plot(tSeq, Yd, type="n", axes=FALSE, ann=FALSE,
#'     xlim=c(0,0.6), ylim=c(0,1))
#'axis(2,at=seq(0,1,0.2),labels=round(seq(0,1,0.2)*max.Yd,digits=2))
#'mtext("Density (1/s)", side=2, line=3)
#'axis(1,at=pretty(c(0,0.6)))
#'mtext("Time (s)", side=1, line=3)
#'axis(4, at=seq(0,1,0.2), labels=round(seq(0,1,0.2)*max.Yh,digits=2))
#'mtext("Hazard (1/s)", side=4, line=3, col=2)
#'mtext("Refractory Exponential Density and Hazard Functions", side=3, line=2,cex=1.5)
#'lines(tSeq,Yd)
#'lines(tSeq,Yh,col=2)
#'par(oldpar)
#'}
#'
drexp <- function(x,
                  rate = 10,
                  rp = 0.005,
                  log = FALSE) {

  ## Check that rate is strictly positive
  if (rate <= 0) stop("rate should be strictly positive")
  ## Check that rp is positive or null
  if (rp < 0) stop("rp should be positive or null")
  x <- x-rp
  dexp(x,rate=rate,log=log)
}




#'Hazard Functions for Some Common Duration Distributions
#'
#'Hazard functions for the gamma, weibull, lognormal, inverse Gaussian, log
#'logistic and refractory exponential distributions
#'
#'These functions are simply obtained by deviding the density by the survival
#'fucntion.
#'
#'@aliases hgamma hweibull hlnorm hinvgauss hllogis hrexp
#'@param x vector of quantiles.
#'@param shape,scale,rate,sdlog strictly positive parameters. See corresponding
#'distributions for detail.
#'@param mu,sigma2,boundary parameters associated with the inverse Gaussian
#'distribution.
#'@param meanlog parameter associated with the log normal distribution.
#'@param location,rp parameters of the log logistic and refratory exponential.
#'@param log should the log hazard be returned? \code{FALSE} by default.
#'@return A vector of hazard rates.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{dinvgauss}}, \code{\link{dllogis}}, \code{\link{drexp}}
#'@references Lindsey, J.K. (2004) \emph{Introduction to Applied Statistics: A
#'Modelling Approach}. OUP.
#'
#'Lindsey, J.K. (2004) \emph{The Statistical Analysis of Stochastic Processes
#'in Time}. CUP.
#'@keywords distribution ts
#'@examples
#'
#'\dontrun{
#'## use a few plots to compare densities and hazard functions
#'
#'## lognormal
#'tSeq <- seq(0.001,0.6,0.001)
#'meanlog.true <- -2.4
#'sdlog.true <- 0.4
#'Yd <- dlnorm(tSeq,meanlog.true,sdlog.true)
#'Yh <- hlnorm(tSeq,meanlog.true,sdlog.true)
#'max.Yd <- max(Yd)
#'max.Yh <- max(Yh)
#'Yd <- Yd / max.Yd
#'Yh <- Yh / max.Yh
#'oldpar <- par(mar=c(5,4,4,4))
#'plot(tSeq, Yd, type="n", axes=FALSE, ann=FALSE,
#'     xlim=c(0,0.6), ylim=c(0,1))
#'axis(2,at=seq(0,1,0.2),labels=round(seq(0,1,0.2)*max.Yd,digits=2))
#'mtext("Density (1/s)", side=2, line=3)
#'axis(1,at=pretty(c(0,0.6)))
#'mtext("Time (s)", side=1, line=3)
#'axis(4, at=seq(0,1,0.2), labels=round(seq(0,1,0.2)*max.Yh,digits=2))
#'mtext("Hazard (1/s)", side=4, line=3, col=2)
#'mtext("Lognormal Density and Hazard Functions", side=3, line=2,cex=1.5)
#'lines(tSeq,Yd)
#'lines(tSeq,Yh,col=2)
#'par(oldpar)
#'
#'## inverse Gaussian
#'tSeq <- seq(0.001,0.6,0.001)
#'mu.true <- 0.075
#'sigma2.true <- 3
#'Yd <- dinvgauss(tSeq,mu.true,sigma2.true)
#'Yh <- hinvgauss(tSeq,mu.true,sigma2.true)
#'max.Yd <- max(Yd)
#'max.Yh <- max(Yh)
#'Yd <- Yd / max.Yd
#'Yh <- Yh / max.Yh
#'oldpar <- par(mar=c(5,4,4,4))
#'plot(tSeq, Yd, type="n", axes=FALSE, ann=FALSE,
#'     xlim=c(0,0.6), ylim=c(0,1))
#'axis(2,at=seq(0,1,0.2),labels=round(seq(0,1,0.2)*max.Yd,digits=2))
#'mtext("Density (1/s)", side=2, line=3)
#'axis(1,at=pretty(c(0,0.6)))
#'mtext("Time (s)", side=1, line=3)
#'axis(4, at=seq(0,1,0.2), labels=round(seq(0,1,0.2)*max.Yh,digits=2))
#'mtext("Hazard (1/s)", side=4, line=3, col=2)
#'mtext("Inverse Gaussian Density and Hazard Functions", side=3, line=2,cex=1.5)
#'lines(tSeq,Yd)
#'lines(tSeq,Yh,col=2)
#'par(oldpar)
#'
#'## gamma
#'tSeq <- seq(0.001,0.6,0.001)
#'shape.true <- 6
#'scale.true <- 0.012
#'Yd <- dgamma(tSeq, shape=shape.true, scale=scale.true)
#'Yh <- hgamma(tSeq, shape=shape.true, scale=scale.true)
#'max.Yd <- max(Yd)
#'max.Yh <- max(Yh)
#'Yd <- Yd / max.Yd
#'Yh <- Yh / max.Yh
#'oldpar <- par(mar=c(5,4,4,4))
#'plot(tSeq, Yd, type="n", axes=FALSE, ann=FALSE,
#'     xlim=c(0,0.6), ylim=c(0,1))
#'axis(2,at=seq(0,1,0.2),labels=round(seq(0,1,0.2)*max.Yd,digits=2))
#'mtext("Density (1/s)", side=2, line=3)
#'axis(1,at=pretty(c(0,0.6)))
#'mtext("Time (s)", side=1, line=3)
#'axis(4, at=seq(0,1,0.2), labels=round(seq(0,1,0.2)*max.Yh,digits=2))
#'mtext("Hazard (1/s)", side=4, line=3, col=2)
#'mtext("Gamma Density and Hazard Functions", side=3, line=2,cex=1.5)
#'lines(tSeq,Yd)
#'lines(tSeq,Yh,col=2)
#'par(oldpar)
#'
#'## Weibull
#'tSeq <- seq(0.001,0.6,0.001)
#'shape.true <- 2.5
#'scale.true <- 0.085
#'Yd <- dweibull(tSeq, shape=shape.true, scale=scale.true)
#'Yh <- hweibull(tSeq, shape=shape.true, scale=scale.true)
#'max.Yd <- max(Yd)
#'max.Yh <- max(Yh)
#'Yd <- Yd / max.Yd
#'Yh <- Yh / max.Yh
#'oldpar <- par(mar=c(5,4,4,4))
#'plot(tSeq, Yd, type="n", axes=FALSE, ann=FALSE,
#'     xlim=c(0,0.6), ylim=c(0,1))
#'axis(2,at=seq(0,1,0.2),labels=round(seq(0,1,0.2)*max.Yd,digits=2))
#'mtext("Density (1/s)", side=2, line=3)
#'axis(1,at=pretty(c(0,0.6)))
#'mtext("Time (s)", side=1, line=3)
#'axis(4, at=seq(0,1,0.2), labels=round(seq(0,1,0.2)*max.Yh,digits=2))
#'mtext("Hazard (1/s)", side=4, line=3, col=2)
#'mtext("Weibull Density and Hazard Functions", side=3, line=2,cex=1.5)
#'lines(tSeq,Yd)
#'lines(tSeq,Yh,col=2)
#'par(oldpar)
#'
#'## refractory exponential
#'tSeq <- seq(0.001,0.6,0.001)
#'rate.true <- 20
#'rp.true <- 0.01
#'Yd <- drexp(tSeq, rate.true, rp.true)
#'Yh <- hrexp(tSeq, rate.true, rp.true)
#'max.Yd <- max(Yd)
#'max.Yh <- max(Yh)
#'Yd <- Yd / max.Yd
#'Yh <- Yh / max.Yh
#'oldpar <- par(mar=c(5,4,4,4))
#'plot(tSeq, Yd, type="n", axes=FALSE, ann=FALSE,
#'     xlim=c(0,0.6), ylim=c(0,1))
#'axis(2,at=seq(0,1,0.2),labels=round(seq(0,1,0.2)*max.Yd,digits=2))
#'mtext("Density (1/s)", side=2, line=3)
#'axis(1,at=pretty(c(0,0.6)))
#'mtext("Time (s)", side=1, line=3)
#'axis(4, at=seq(0,1,0.2), labels=round(seq(0,1,0.2)*max.Yh,digits=2))
#'mtext("Hazard (1/s)", side=4, line=3, col=2)
#'mtext("Refractory Exponential Density and Hazard Functions", side=3, line=2,cex=1.5)
#'lines(tSeq,Yd)
#'lines(tSeq,Yh,col=2)
#'par(oldpar)
#'
#'## log logistic
#'tSeq <- seq(0.001,0.6,0.001)
#'location.true <- -2.7
#'scale.true <- 0.025
#'Yd <- dllogis(tSeq, location.true, scale.true)
#'Yh <- hllogis(tSeq, location.true, scale.true)
#'max.Yd <- max(Yd)
#'max.Yh <- max(Yh)
#'Yd <- Yd / max.Yd
#'Yh <- Yh / max.Yh
#'oldpar <- par(mar=c(5,4,4,4))
#'plot(tSeq, Yd, type="n", axes=FALSE, ann=FALSE,
#'     xlim=c(0,0.6), ylim=c(0,1))
#'axis(2,at=seq(0,1,0.2),labels=round(seq(0,1,0.2)*max.Yd,digits=2))
#'mtext("Density (1/s)", side=2, line=3)
#'axis(1,at=pretty(c(0,0.6)))
#'mtext("Time (s)", side=1, line=3)
#'axis(4, at=seq(0,1,0.2), labels=round(seq(0,1,0.2)*max.Yh,digits=2))
#'mtext("Hazard (1/s)", side=4, line=3, col=2)
#'mtext("Log Logistic Density and Hazard Functions", side=3, line=2,cex=1.5)
#'lines(tSeq,Yd)
#'lines(tSeq,Yh,col=2)
#'par(oldpar)
#'}
#'
hgamma <- function(x,
                   shape,
                   rate = 1,
                   scale = 1/rate,
                   log = FALSE) {

  if (!log) dgamma(x,shape,,scale)/pgamma(x,shape,,scale,lower.tail=FALSE)
  else dgamma(x,shape,,scale,log=TRUE) -
    pgamma(x,shape,,scale,lower.tail=FALSE,log.p=TRUE)

}

hinvgauss <- function(x,
                      mu = 1,
                      sigma2 = 1,
                      boundary = NULL,
                      log = FALSE) {

  if( any(mu <= 0) ) stop("mu must be positive")
  if (all(!is.null(sigma2)))
    if( any(sigma2 <= 0) ) stop("sigma2 must be positive")
  if (all(!is.null(boundary)))
    if( any(boundary <= 0) ) stop("boundary must be positive")
  if ( all(!is.null(sigma2)) && all(!is.null(boundary)) )
    stop("One of sigma2 or boundary must be specified, not both")

  if (all(!is.null(boundary))) {
    ## We work with the parameterization in term of boundary (ie, sigma2 = 1)
    ## We convert it in term of mu and sigma2
    sigma2 <- (1/boundary)^2
    mu <- boundary*mu
  }
  
  t <- x/mu
  v <- sqrt(x*sigma2)
  cutOff <- 1e-12
  bigEnough <- pinvgauss(x, mu, sigma2, lower.tail = FALSE) > cutOff
  logIntensity <- -( (t[bigEnough]-1)^2/(x[bigEnough]*sigma2) + log(2*sigma2*pi*x[bigEnough]^3) )/2 -
    log(1-pnorm((t[bigEnough]-1)/v[bigEnough])-exp(2/(mu*sigma2))*pnorm(-(t[bigEnough]+1)/v[bigEnough])
        )
  if (sum(bigEnough) < length(x)) {
    cutOffQ <- qinvgauss(1 - cutOff, mu, sigma2)
    cutOffValue <- dinvgauss(cutOffQ, mu, sigma2, log = TRUE) - log(cutOff)
    logIntensity <- c(logIntensity,rep(cutOffValue,length(x) - sum(bigEnough))) 
  }
  if (!log) return(exp(logIntensity))
  else return(logIntensity)
  
}

hllogis <- function(x,
                    location = 0,
                    scale = 1,
                    log = FALSE) {
  
  ## Check that x elements are strictly positive
  if (any(x <= 0)) stop("x elements must be strictly positive")

  if (!log) pllogis(x,location,scale)/(scale*x)
  else pllogis(x,location,scale,log.p=TRUE) - log(scale) - log(x)

}

hlnorm <- function(x,
                   meanlog = 0,
                   sdlog = 1,
                   log = FALSE) {
  ## Check that y elements are strictly positive
  if (any(x <= 0))
    stop(paste("The elements of",
               deparse(substitute(x)),
               "should be strictly positive.")
         )

  if (!log) dlnorm(x,meanlog,sdlog)/plnorm(x,meanlog,sdlog,lower.tail=FALSE)
  else dlnorm(x,meanlog,sdlog,log=TRUE)-
    plnorm(x,meanlog,sdlog,lower.tail=FALSE,log.p=TRUE)
  
}

hrexp <- function(x,
                  rate = 10,
                  rp = 0.005,
                  log = FALSE) {

  ## Check that rate is strictly positive
  if (rate <= 0) stop("rate should be strictly positive")
  ## Check that rp is positive or null
  if (rp < 0) stop("rp should be positive or null")
  
  if (!log) ifelse(x >= rp, rate, 0)
  else ifelse(x >= rp, log(rate), -Inf)
}


hweibull <- function(x,
                     shape,
                     scale = 1,
                     log = FALSE) {

  if (!log) dweibull(x, shape, scale) / pweibull(x,shape,scale,lower.tail=FALSE)
  else dweibull(x,shape,scale,log=TRUE) - pweibull(x,shape,scale,lower.tail=FALSE,log.p=TRUE)
  
}



#'ISI Histogram With Fitted Model and CI
#'
#'Fits a duration model to isis from a spike train. Confidence intervals are
#'also drawn.
#'
#'Assuming that the train is reasonably well described by a renewal process, a
#'\code{model} distribution is fitted to the inter-spike intervals (isis)
#'obtained from \code{spikeTrain}. The fitted distribution is then used to set
#'the histogram breaks such that a uniform bin count would be expected if the
#'fitted distribution was the true one. Confidence segments are also obtained
#'from the binomial distribution. The histogram is build and the fitted density
#'together with confidence intervals are drawn.
#'
#'@param spikeTrain a \code{spikeTrain} object or a numeric vector that can be
#'coerced to such an object.
#'@param model a character vector whose elements are selected among:
#'\code{"invgauss"}, \code{"lnorm"}, \code{"gamma"}, \code{"weibull"},
#'\code{"llogis"}, \code{"rexp"}.
#'@param nbins the number of bins to use.
#'@param CI the confidence coefficient.
#'@param \dots additional arguments passed to \code{hist}, see
#'\code{\link{hist}}.
#'@return Nothing returned, \code{isiHistFit} is used for its side effect, a
#'plot is generated on the current graphic device.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{compModels}}, \code{\link{hist}}
#'@keywords distribution ts
#'@examples
#'
#'\dontrun{
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
#'## It does not look too bad so let fit simple models
#'compModels(CAL1S[["neuron 1"]])
#'## the best one is the invgauss. Let's look at
#'## it in detail
#'isiHistFit(CAL1S[["neuron 1"]],"invgauss",xlim=c(0,0.5))
#'}
#'
isiHistFit <- function(spikeTrain,
                       model,
                       nbins = 10,
                       CI=0.95,
                       ...) {

  spikeTrainName <- deparse(substitute(spikeTrain))
  if (!is.spikeTrain(spikeTrain)) spikeTrain <- as.spikeTrain(spikeTrain)

  isi <- diff(spikeTrain)
  knownModels <- c("invgauss", "lnorm", "gamma", "weibull", "llogis", "rexp")
  model <- model[1]
  if (!model %in% knownModels) 
    stop(paste(deparse(substitute(model)),
               "is not implemented.")
         )

  CI <- CI[1]
  if (CI <= 0 | CI >= 1)
    stop(paste("A CI of", CI, "does not make sense."))

  theFit <- switch(model,
                   invgauss = invgaussMLE(isi),
                   lnorm = lnormMLE(isi),
                   gamma = gammaMLE(isi),
                   weibull = weibullMLE(isi),
                   llogis = llogisMLE(isi),
                   rexp = rexpMLE(isi)
                   )

  pSeq <- (1:(nbins-1))/nbins

  theEstimates <- theFit$estimate
  qFct <- switch(model,
                 invgauss = function(p) qinvgauss(p,theEstimates[1],theEstimates[2]),
                 lnorm = function(p) qlnorm(p,theEstimates[1],theEstimates[2]),
                 gamma = function(p) qgamma(p,theEstimates[1],theEstimates[2]),
                 weibull = function(p) qweibull(p,theEstimates[1],theEstimates[2]),
                 llogis = function(p) qllogis(p,theEstimates[1],theEstimates[2]),
                 rexp = function(p) qrexp(p,theEstimates[1],theEstimates[2])
                 )

  breaks <- c(0,qFct(pSeq),max(isi)+.Machine$double.eps)
  mids <- breaks[-(nbins+1)] + diff(breaks)/2
  main <- paste("Isi histogram and fitted",
                model,
                "distribution for",
                spikeTrainName
                )
  sub <- paste("CI at ",
               CI*100,
               "%. Sample size: ",
               length(isi),".",
               sep=""
               )
  hist(isi,
       breaks=breaks,
       main=main,
       sub=sub,
       xlab="isi (s)",
       ...)

  X <- seq(.Machine$double.eps,breaks[nbins+1],length.out=501)

  dFct <- switch(model,
                 invgauss = function(x) dinvgauss(x,theEstimates[1],theEstimates[2]),
                 lnorm = function(x) dlnorm(x,theEstimates[1],theEstimates[2]),
                 gamma = function(x) dgamma(x,theEstimates[1],theEstimates[2]),
                 weibull = function(x) dweibull(x,theEstimates[1],theEstimates[2]),
                 llogis = function(x) dllogis(x,theEstimates[1],theEstimates[2]),
                 rexp = function(x) drexp(x,theEstimates[1],theEstimates[2])
                 )
  Y <- dFct(X)
  lines(X,Y,col=2)
  ci <- qbinom(c((1-CI)/2,1-(1-CI)/2),size=length(isi),prob=1/nbins)
  l <- ci[1]/(length(isi)*diff(breaks))
  u <- ci[2]/(length(isi)*diff(breaks))
  invisible(sapply(1:nbins,function(idx) segments(mids[idx],l[idx],mids[idx],u[idx],col=2)))
  
}

pinvgauss <- function(q,
                      mu = 1,
                      sigma2 = 1,
                      boundary = NULL,
                      lower.tail = TRUE,
                      log.p = FALSE
                      ){
  
  if( any(q <= 0) ) stop("q must contain positive values")
  if( any(mu <= 0) ) stop("mu must be positive")
  if (all(!is.null(sigma2)))
    if( any(sigma2 <= 0) ) stop("sigma2 must be positive")
  if (all(!is.null(boundary)))
    if( any(boundary <= 0) ) stop("boundary must be positive")
  if ( all(!is.null(sigma2)) && all(!is.null(boundary)) )
    stop("One of sigma2 or boundary must be specified, not both")

  if (all(!is.null(boundary))) {
    ## We work with the parameterization in term of boundary (ie, sigma2 = 1)
    ## We convert it in term of mu and sigma2
    sigma2 <- (1/boundary)^2
    mu <- boundary*mu
  }
  t <- q/mu
  v <- sqrt(q*sigma2)

  ## Use Eq. 4 of Whitemore GA and Seshadri V (1987)
  ## The American Statistician 41:280-281
  if (lower.tail & !log.p)
    return(pnorm((t-1)/v)+exp(2/(mu*sigma2))*pnorm(-(t+1)/v))
  if (!lower.tail & !log.p)
    return(1 - (pnorm((t-1)/v)+exp(2/(mu*sigma2))*pnorm(-(t+1)/v)))
  if (lower.tail & log.p)
    return(log(pnorm((t-1)/v)+exp(2/(mu*sigma2))*pnorm(-(t+1)/v)))
  if (!lower.tail & log.p)
    return(log(1 - (pnorm((t-1)/v)+exp(2/(mu*sigma2))*pnorm(-(t+1)/v))))
  
}

pllogis <- function(q,
                    location = 0,
                    scale = 1,
                    lower.tail = TRUE,
                    log.p = FALSE) {
  
  plogis(log(q),location,scale,lower.tail,log.p)
  
}

prexp <- function(q,
                  rate = 10,
                  rp = 0.005,
                  lower.tail = TRUE,
                  log.p = FALSE) {

  ## Check that rate is strictly positive
  if (rate <= 0) stop("rate should be strictly positive")
  ## Check that rp is positive or null
  if (rp < 0) stop("rp should be positive or null")
  q <- q-rp
  pexp(q,rate,lower.tail,log.p)
}

qinvgauss <- function(p,
                      mu = 1,
                      sigma2 = 1,
                      boundary = NULL
                      ){

  if( any(p < 0 | p > 1) ) stop("p must lie between 0 and 1")
  if( any(mu <= 0) ) stop("mu must be positive")
  if (all(!is.null(sigma2)))
    if( any(sigma2 <= 0) ) stop("sigma2 must be positive")
  if (all(!is.null(boundary)))
    if( any(boundary <= 0) ) stop("boundary must be positive")
  if ( all(!is.null(sigma2)) && all(!is.null(boundary)) )
    stop("One of sigma2 or boundary must be specified, not both")

  if (all(!is.null(boundary))) {
    ## We work with the parameterization in term of boundary (ie, sigma2 = 1)
    ## We convert it in term of mu and sigma2
    sigma2 <- (1/boundary)^2
    mu <- boundary*mu
  }
  
  len <- max(length(p),length(mu),length(sigma2))

  if(length(p) != len) {
    if(length(p) == 1) p <- rep(p,len)
    else stop("length of p incorrect")
  }
  if(length(mu) != len) {
    if(length(mu) == 1) mu <- rep(mu,len)
    else stop("length of m incorrect")
  }
  if(length(sigma2) != len) {
    if(length(sigma2) == 1) sigma2 <- rep(sigma2,len)
    else stop("length of sigma2 incorrect")
  }

  ## Use Whitemore and Yalovky (1978, Technometrics, 20:207-208)
  ## approximation to get starting value for the numerical
  ## inversion of the cumulative distribution function.
  theta <- 1/mu/sigma2
  approx <- mu * exp(qnorm(p)*sqrt(1/theta)-0.5/theta)
  sapply(1:len, function(idx) {
    if (identical(p[idx],0)) return(0)
    if (identical(p[idx],1)) return(Inf)
    interval <- approx[idx]*c(0.95,1.05)
    h <- function(q) pinvgauss(q, mu[idx], sigma2[idx]) - p[idx]
    while (h(interval[1])*h(interval[2]) > 0)
      interval <- interval*c(0.9,1.1)
    uniroot(h,interval)$root
  }
         )
  
}

qllogis <- function(p,
                    location = 0,
                    scale = 1,
                    lower.tail = TRUE,
                    log.p = FALSE) {
  
  exp(qlogis(p,location,scale,lower.tail,log.p))
  
}

qrexp <- function(p,
                  rate = 10,
                  rp = 0.005,
                  lower.tail = TRUE,
                  log.p = FALSE) {

  ## Check that rate is strictly positive
  if (rate <= 0) stop("rate should be strictly positive")
  ## Check that rp is positive or null
  if (rp < 0) stop("rp should be positive or null")
  rp + qexp(p,rate,lower.tail,log.p)

}


rinvgauss <- function(n = 1,
                      mu = 1,
                      sigma2 = 1,
                      boundary = NULL
                      ){

  if( any(mu <= 0) ) stop("mu must be positive")
  if (all(!is.null(sigma2)))
    if( any(sigma2 <= 0) ) stop("sigma2 must be positive")
  if (all(!is.null(boundary)))
    if( any(boundary <= 0) ) stop("boundary must be positive")
  if ( all(!is.null(sigma2)) && all(!is.null(boundary)) )
    stop("One of sigma2 or boundary must be specified, not both")

  if (all(!is.null(boundary))) {
    ## We work with the parameterization in term of boundary (ie, sigma2 = 1)
    ## We convert it in term of mu and sigma2
    sigma2 <- (1/boundary)^2
    mu <- boundary*mu
  }

  ## Use method of Michael JR, Schucany WR and Haas RW (1976)
  ## The American Statistician, 30:88-90
  v0 <- rchisq(n,1)
  x1 <- mu + 0.5*mu^2*sigma2*v0 - 0.5*mu*sigma2*sqrt(4*mu*v0/sigma2+mu^2*v0^2)
  ifelse(rbinom(length(x1),1,mu/(mu+x1)) == 1,x1,mu^2/x1)
}


rllogis <- function(n,
                    location = 0,
                    scale = 1) {
  
  exp(rlogis(n,location,scale))

}

rrexp <- function(n,
                  rate = 10,
                  rp = 0.005) {
  ## Check that rate is strictly positive
  if (rate <= 0) stop("rate should be strictly positive")
  ## Check that rp is positive or null
  if (rp < 0) stop("rp should be positive or null")
  rp+rexp(n,rate)
}
