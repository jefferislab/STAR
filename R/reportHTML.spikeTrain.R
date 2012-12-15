#'Performs Basic Spike Train Analysis and Generates a Report in HTML Format
#'from a spikeTrain Object
#'
#'Performs a "standard" analysis on a \code{spikeTrain} object, computing some
#'cross-correlation statistics if additional \code{spikeTrain} objects are
#'provided, writes results to disk and generates a report in html format.
#'
#'A spike train plot (\code{\link{plot.spikeTrain}}) is performed first. The
#'summary (\code{\link{summary.spikeTrain}}) is computed next and part of its
#'output is written to the html file. The renewal tests are then carried out
#'and their results added (\code{\link{renewalTestPlot}}). The six duration
#'distributions are fitted (\code{\link{compModels}} with argument \code{plot}
#'set to \code{FALSE}) and the best one is used to apply a time transformation
#'to \code{spikeTrain}. The Ogata's tests are applied
#'(\code{\link{summary.transformedTrain}}) and if they are all within the 99\%
#'confidence interval, the result of the transformation is plotted
#'(\code{\link{plot.transformedTrain}}) as well as all the Q-Q plots of
#'\code{\link{compModels}}. If \code{forceTT} is set to \code{TRUE} (default),
#'then these last two plots are added even if the best model does not pass the
#'tests.
#'
#'If other \code{spikeTrain} objects are provided as a named list via argument
#'\code{otherST}, then cross-correlation/cross-intensity functions are
#'estimated; Two estimations methods are available, the classical histogram and
#'a smooth version of it. Argument \code{cch} controls if a single estimation
#'is performed or if both are performed. If the smooth version is requested a
#'summary of the \code{\link[gss]{gssanova}}, \code{\link[gss]{gssanova0}} or
#'\code{\link[mgcv]{gam}} fit is printed (depending on the chosen value for
#'argument \code{method}). Moreover if argument \code{doGamCheck} is set to
#'\code{TRUE} (and if \code{method} is set to \code{\link{gamlockedTrain}})
#'then check plots (\code{\link[mgcv]{gam.check}}) are added to the report.
#'
#'A \code{R} data file (\code{filename.rda}) is also generated with the
#'following objects: \itemize{ \item \code{cm}: the result of
#'\code{\link{compModels}}.  \item \code{bestFit}: the
#'\code{durationDistribution} object returned obtained by fitting the best
#'model among the 6.  \item \code{Lambda}: the integrated intensity of
#'\code{spikeTrain} with the best model.  \item \code{fct}: the matched call.
#'\item \code{cchL}: if other trains were provided and if argument \code{cch}
#'was set to \code{"both"} or to \code{"cch"}. A list with as many components
#'as the \code{otherST} argument. Each component is the a
#'\code{hist.lockedTrain} object.  \item \code{scchL}: if other trains were
#'provided and if argument \code{cch} was set to \code{"both"} or to
#'\code{"scch"}. A list with as many components as the \code{otherST} argument.
#'Each component is the a \code{gsslockedTrain}, \code{gsslockedTrain0} or
#'\code{gamlockedTrain} object.  }
#'
#'@param object a \code{spikeTrain} object.
#'@param filename a character string. The generic name of all the files (html,
#'png as well as \code{R} data files which will be generated. See also
#'\code{\link[R2HTML]{HTMLInitFile}}.
#'@param extension see \code{\link[R2HTML]{HTMLInitFile}}.
#'@param directory the full or relative path to the directory where the results
#'are going to be stored. See also \code{\link[R2HTML]{HTMLInitFile}}.
#'@param Title See \code{\link[R2HTML]{HTMLInitFile}}. If missing a default
#'value baed on \code{filename} is provided.
#'@param forceTT Should a time transformation be performed and the
#'\code{\link{compModels}} plots be generated even if none of the six renewal
#'models fits the data?
#'@param timeUnit,digits see \code{\link{summary.spikeTrain}}.
#'@param otherST a named list of \code{spikeTrain} objects from simultaneously
#'recorded neurons or nothing.
#'@param laglim see \code{\link{lockedTrain}}.
#'@param cch if \code{otherST} is given (ie, not missing) cross-intensity plots
#'will be made using the neuron of \code{spikeTrain} as a reference. Should
#'smooth version of the cross-intensity be computed (\code{"scch"}), a
#'"classical" one (\code{"cch"}) or both (\code{"both"}). Only the first
#'element of \code{cch} is used.
#'@param method A character string, the name of the function used to generate
#'the smooth cross-correlation histograms, one of:
#'\code{\link{gsslockedTrain}}, \code{\link{gsslockedTrain0}},
#'\code{\link{gamlockedTrain}}.
#'@param doGamCheck if smooth estimates are requested and \code{method} is set
#'to \code{\link{gamlockedTrain}}, should function
#'\code{\link[mgcv]{gam.check}} be used on them?
#'@param k see \code{\link{gamlockedTrain}}.
#'@param bs see \code{\link{gamlockedTrain}}.
#'@param nbEvtPerBin a number of event per bin used in a way similar to the
#'argument with the same name in \code{\link{jpsth}} when a bining is used for
#'pre-processing.
#'@param \dots Passed to \code{\link{gsslockedTrain}},
#'\code{\link{gsslockedTrain0}}, \code{\link{gamlockedTrain}}.
#'@return Nothing is returned, an html file and figures in png format are
#'written to disk together with the \code{R} variables generated during the
#'analysis.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{as.spikeTrain}}, \code{\link{plot.spikeTrain}},
#'\code{\link{summary.spikeTrain}}, \code{\link{renewalTestPlot}},
#'\code{\link{plot.spikeTrain}}, \code{\link{compModels}},
#'\code{\link{transformedTrain}}, \code{\link{plot.transformedTrain}},
#'\code{\link{summary.transformedTrain}}, \code{\link[gss]{gssanova}},
#'\code{\link[gss]{gssanova0}}, \code{\link[mgcv]{gam}},
#'\code{\link[mgcv]{gam.check}}, \code{\link{lockedTrain}},
#'\code{\link{gsslockedTrain}}, \code{\link{gsslockedTrain0}},
#'\code{\link{gamlockedTrain}}
#'@keywords models smooth regression
#'@examples
#'
#'\dontrun{
#'## load e070528spont data set
#'data(e070528spont)
#'## perform a standard analysis on neuron 1, looking for cross-correlations
#'## with the 3 other neurons up to lag +/- 250 ms.
#'## Store the results under the generic name: e070528spontN1
#'reportHTML(e070528spont[["neuron 1"]],"e070528spontN1",otherST=e070528spont[-1],laglim=c(-1,1)*0.25,forceTT=FALSE)
#'## Neuron 1 of e070528spont is exceptional in that it can be well
#'## described by a renewal process...}
#'
reportHTML.spikeTrain <- function(object,
                                  filename,
                                  extension="html",
                                  directory=getwd(),
                                  Title,
                                  forceTT=TRUE,
                                  digits=3,
                                  timeUnit="s",
                                  otherST,
                                  laglim=c(-0.1,0.1),
                                  cch=c("both","scch","cch"),
                                  method=c("gsslockedTrain0","gsslockedTrain","gamlockedTrain"),
                                  doGamCheck=FALSE,
                                  k=100,
                                  bs="tp",
                                  nbEvtPerBin=10,
                                  ...
                                  )
### reportHTML method for spikeTrain objects.
### A spike train plot of object is generated first
### by a call to plot.spikeTrain. A short numerical
### summary is generated next (partial output of a call
### to summary.spikeTrain). A renewal test is performed
### by a call to renewalTestPlot using default arguments
### of the latter. The six duration distributions with
### two paramters are fitted with a call to compModels with
### the plot argument set to FALSE and a call to xMLE where
### "x" stands for the best model (in the AIC sense) among the
### the 6. A summary (best value and se of the two model 
### parameters) of the best fit is printed. A time transformation using:
### cumsum(px(diff(y))
### where "x" as the same meaning as above and y stands for the
### time transformed (TT) version of the original train using the
### best among the 6 models. If this TT version passes the "tests
### of Ogata" at the 99% level (these tests are returned by a call
### to summary.transformedTrain), then a plot is generated by a
### new call to compModels displaying the QQ plot on a log scale
### for each of the 6 fitted models as well as the Ogata's test
### plots. If arg forceTT is set to TRUE, these plots are generated
### even if the tests are not passed.
### Args "digits" and "timeUnit" are the corresponding args of
### summary.spikeTrain
{

  objectN <- deparse(substitute(object))
  
  ## check is object is a spikeTrain object
  if (!is.spikeTrain(object)) object <- as.spikeTrain(object)

  if (missing(filename))
    filename <- paste(objectN,"analysis")
  
  if (missing(Title))
    Title <- filename

  
  HTMLInitFile(outdir=directory,
               filename=filename,
               extension=extension,
               Title=Title)

  fullName <- paste(directory,"/",
                    filename,".",
                    extension,sep="")

  saveName <- paste(directory,"/",
                    filename,".rda",
                    sep="")
  
  HTML.title(filename,
             file=fullName,
             HR=2)
  
  ## add a spike train plot 
  HTML.title(paste("Spike train plot of ",
                   filename,sep=""),
             file=fullName,HR=3)

  
  stFigName <- paste(filename,"_st.png",sep="")
  figFname <- paste(directory,"/",stFigName,sep="")
  png(figFname,width=500,height=500)
  plot(object)
  dev.off()
  HTMLInsertGraph(stFigName,
                  file=fullName,
                  WidthHTML=500,
                  HeightHTML=500)

  HTMLbr(2,file=fullName)
  HTML.title(paste("Short summary of ",
                   filename,sep=""),
             file=fullName,HR=3)
  stRange <- range(object)
  stNb <- length(object)
  isi <- diff(object)
  stStat1 <- c(mean(isi), sd(isi))
  stStat2 <- c(mean(log(isi)), sd(log(isi)))
  cat(paste("<p>A spike train with ", stNb, " events, starting at: ", 
            round(stRange[1], digits = digits), " and ending at: ", 
            round(stRange[2], digits = digits), " (", timeUnit, ").</p>", 
            "<p>The mean ISI is: ", round(stStat1[1], digits = digits), 
            " and its SD is: ", round(stStat1[2], digits = digits), 
            " (", timeUnit, ").</p>", "<p>The shortest interval is: ", 
            round(min(isi), digits = digits), " and the longest is: ", 
            round(max(isi), digits = digits), " (", timeUnit, ").</p>", 
            sep = ""),
      file=fullName,
      append=TRUE)
  
  ## add a renewal test plot if there are more than 50 events
  if (length(object) > 50) {
    HTMLhr(file=fullName)
    HTML.title(paste("Renewal test of ",
                     filename,sep=""),
               file=fullName,HR=3)
    
    rtFigName <- paste(filename,"_rt.png",sep="")
    figFname <- paste(directory,"/",rtFigName,sep="")
    png(figFname,width=800,height=800)
    renewalTestPlot(object)
    dev.off()
    HTMLInsertGraph(rtFigName,
                    file=fullName,
                    WidthHTML=800,
                    HeightHTML=800)
  } ## End of conditional on length(object) > 50

  ## compare models
  cm <- compModels(object,plot=FALSE)
  ## get best model
  bestM <- names(cm)[1]
  HTMLhr(file=fullName)
  cat(paste("<p>The best model with 2 parameters is: ",
            bestM,
            ".</p>",sep=""),
      file=fullName,
      append=TRUE
      )
  cat("<p>Its estimated parameters values and associated se are:</p>",
      file=fullName,
      append=TRUE
      )
  
  bestFit <- switch(bestM,
                    invgauss=invgaussMLE(object),
                    lnorm=lnormMLE(object),
                    gamma=gammaMLE(object),
                    weibull=weibullMLE(object),
                    llogis=llogisMLE(object),
                    rexp=rexpMLE(object)
                    )
  
  toPrint <- rbind(bestFit$estimate,bestFit$se)
  rownames(toPrint) <- c("mean","se")
  HTML(toPrint,file=fullName)
  
  Lambda <- switch(bestM,
                   invgauss=pinvgauss(isi,bestFit$estimate[1],
                     bestFit$estimate[2],log.p=TRUE,lower.tail=FALSE),
                   lnorm=plnorm(isi,bestFit$estimate[1],
                     bestFit$estimate[2],log.p=TRUE,lower.tail=FALSE),
                   gamma=pgamma(isi,bestFit$estimate[1],
                     bestFit$estimate[2],log.p=TRUE,lower.tail=FALSE),
                   weibull=pweibull(isi,bestFit$estimate[1],
                     bestFit$estimate[2],log.p=TRUE,lower.tail=FALSE),
                   llogis=pllogis(isi,bestFit$estimate[1],
                     bestFit$estimate[2],log.p=TRUE,lower.tail=FALSE),
                   rexp=prexp(isi,bestFit$estimate[1],
                     bestFit$estimate[2],log.p=TRUE,lower.tail=FALSE)
                   )
  Lambda <- -cumsum(Lambda)
  class(Lambda) <- c("transformedTrain","spikeTrain")
  
  if (max(Lambda)/50 > 2) {
    LambdaS <- summary(Lambda)

    otTRUE <- LambdaS[[1]][2] &&
    LambdaS[[2]][2] &&
    (0.005 <= LambdaS[[3]][2]) &&
    (0.995 >= LambdaS[[3]][2])
  } else {
    otTRUE <- FALSE
  } ## End of conditional on max(Lambda)/50 > 2

  HTMLbr(1,file=fullName)
  HTML.title(paste("Model comparison for ",
                   filename,sep=""),
             file=fullName,HR=3)
  cmFigName <- paste(filename,"_cm.png",sep="")
  figFname <- paste(directory,"/",cmFigName,sep="")
  png(figFname,width=800,height=800)
  compModels(object)
  dev.off()
  HTMLInsertGraph(cmFigName,
                  file=fullName,
                  WidthHTML=800,
                  HeightHTML=800)
  
  if (forceTT) {

    HTMLbr(3,file=fullName)
    HTML.title(paste("Ogata's tests after time transformation of ",
                     filename,
                     " using a ",bestM," model",
                     sep=""),
             file=fullName,HR=3)
  
    otFigName <- paste(filename,"_TTot.png",sep="")
    figFname <- paste(directory,"/",otFigName,sep="")
    png(figFname,width=800,height=800)
    if (max(Lambda)/50 > 2) {
      plot(Lambda,
           which=c(1,2,4,5),
           ask=FALSE)
    } else {
      plot(Lambda,
           which=c(1,2,4),
           ask=FALSE)
    }
    dev.off()
    HTMLInsertGraph(otFigName,
                    file=fullName,
                    WidthHTML=800,
                    HeightHTML=800)
  } ## End of conditional on otTRUE 

  ## check if otherST is given
  if (!missing(otherST)) {

    keepGoing <- TRUE
    ## check that otherST is a named list of spikeTrain objects
    if (!inherits(otherST,"list") || is.null(names(otherST))) {
      warning("otherST should be a named list of spikeTrain objects.")
      keepGoing <- FALSE
    }
    
    if (keepGoing) {

      if (cch[1] %in% c("both","scch")) {
        scchL <- lapply(otherST,
                        function(st) {
                          lt <- lockedTrain(object,st,laglim=laglim)
                          testF <- lt$nbTestSpikes/lt$obsTime
                          nRef <- lt$nbRefSpikes
                          theBW <- max(0.001,round(nbEvtPerBin/testF/nRef,digits=3))
                          k <- min(k,diff(lt$laglim)%/%theBW)
                          switch(method[1],
                                 gsslockedTrain=gsslockedTrain(lt,bw=theBW,...),
                                 gsslockedTrain0=gsslockedTrain0(lt,bw=theBW,...),
                                 gamlockedTrain=gamlockedTrain(lt,bw=theBW,bs=bs,k=k,...)
                                 )
                        }
                        )
      } ## End of conditional on cch[1] %in% c("both","scch")

      if (cch[1] %in% c("both","cch")) {
        cchL <- lapply(otherST,
                       function(st) {
                         lt <- lockedTrain(object,st,laglim=laglim)
                         testF <- lt$nbTestSpikes/lt$obsTime
                         nRef <- lt$nbRefSpikes
                         theBW <- round(nbEvtPerBin*1.5/testF/nRef,digits=4)
                         hist(lt,bw=theBW,plot=FALSE)
                        }
                        )
      } ## End of conditional on cch[1] %in% c("both","cch")

      sapply(seq(length(otherST)),
             function(trainIdx) {

               ## Write general header
               HTMLbr(1,file=fullName)
               HTMLhr(file=fullName)
               HTML.title(paste("Cross-Intensity with ",
                                names(otherST)[trainIdx],
                                " as a test train",sep=""),
                          file=fullName,HR=3)
               
               if (cch[1] %in% c("both","scch")) {
                 ## A smoothed cross-intensity was evaluated
                 ## write down its summary
                 scch.txt <- switch(method[1],
                                    gsslockedTrain="gssanova fit summary",
                                    gsslockedTrain0="gssanova0 fit summary",
                                    gamlockedTrain="gam fit summary"
                                    )
                 scch.txt <- paste(scch.txt,
                                   " (pre-binning bin width: ",
                                   scchL[[trainIdx]][["bwV"]][1],
                                   "):",sep="")
      
                 HTML.title(scch.txt,
                            file=fullName,HR=4)
                 fitS <- summary(scchL[[trainIdx]])
                 HTML(fitS,file=fullName)

                 if (method[1] == "gamlockedTrain" && doGamCheck) {  
                   HTMLbr(1,file=fullName)
                   HTML.title("GAM goodness of fit diagnostics:",,
                              file=fullName,HR=4)
                   gcFigName <- paste(filename,
                                      "_",
                                      paste(strsplit(names(otherST)[trainIdx]," ")[[1]],collapse="_"),
                                      "_gc.png",sep="")
                   figFname <- paste(directory,"/",gcFigName,sep="")
                   png(figFname,width=800,height=800)
                   gam.check(scchL[[trainIdx]][["gamFit"]])
                   dev.off()
                   HTMLInsertGraph(gcFigName,
                                   file=fullName,
                                   WidthHTML=800,
                                   HeightHTML=800)
                 } ## End of conditional on doGamCheck
                 
                 if (cch[1] %in% c("both")) {
                   ## plot both cross-intensity estiamtes
                   HTML.title("Smoothed and \"classical\" cross-intensity plots:",
                              file=fullName,HR=4)
                   ciFigName <- paste(filename,
                                      "_",
                                      paste(strsplit(names(otherST)[trainIdx]," ")[[1]],collapse="_"),
                                      "_ci.png",sep="")
                   figFname <- paste(directory,"/",ciFigName,sep="")
                   png(figFname,width=1000,height=700)
                   layout(matrix(1:2,1,2))
                   plot(scchL[[trainIdx]],main="")
                   plot(cchL[[trainIdx]],main="")
                   dev.off()
                   HTMLInsertGraph(ciFigName,
                                   file=fullName,
                                   WidthHTML=1000,
                                   HeightHTML=700)
                 } else {
                   ## plot only the smoothed intensity estimate
                   HTML.title("Smoothed cross-intensity plot:",
                              file=fullName,HR=4)
                   ciFigName <- paste(filename,
                                      "_",
                                      paste(strsplit(names(otherST)[trainIdx]," ")[[1]],collapse="_"),
                                      "_ci.png",sep="")
                   figFname <- paste(directory,"/",ciFigName,sep="")
                   png(figFname,width=500,height=500)
                   plot(scchL[[trainIdx]],main="")
                   dev.off()
                   HTMLInsertGraph(ciFigName,
                                   file=fullName,
                                   WidthHTML=500,
                                   HeightHTML=500)
                 } ## End of conditional on cch[1] %in% c("both")
               } else {
                 ## plot the classical cross-intensity estimate
                 HTML.title("\"Classical\" cross-intensity plots:",
                            file=fullName,HR=4)
                 ciFigName <- paste(filename,
                                    "_",names(otherST)[trainIdx],
                                    "_ci.png",sep="")
                 figFname <- paste(directory,"/",ciFigName,sep="")
                 png(figFname,width=500,height=500)
                 plot(cchL[[trainIdx]],main="")
                 dev.off()
                 HTMLInsertGraph(ciFigName,
                                 file=fullName,
                                 WidthHTML=500,
                                 HeightHTML=500)
               } ## End of conditional on cch[1] %in% c("both","scch")
             }
             )
      
    } ## End of conditional on keepGoing
    
  } ## End of conditional on !missing(otherST)
  
  HTMLEndFile()

  fctCall <- match.call()

  if (!exists("scchL") && !exists("cchL")) { 
    save(cm,bestFit,Lambda,fctCall,file=saveName)
  } else {
    if (exists("scchL") && exists("cchL")) {
      save(cm,bestFit,Lambda,fctCall,scchL,cchL,file=saveName)
    } else {
      if (exists("scchL")) {
        save(cm,bestFit,Lambda,fctCall,scchL,file=saveName)
      } else {
        save(cm,bestFit,Lambda,fctCall,cchL,file=saveName)
      } ## End of conditional on exists("scchL")
    } ## End of conditional on exists("scchL") && exists("cchL")
  } ## End of conditional on !exists("scchL") && !exists("cchL")
  
}
