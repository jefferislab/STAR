#'Generates a Report in HTML Format from a STAR gam Object
#'
#'Writes the result of a \code{gam} fit in an html file.
#'
#'A summary (\code{\link[mgcv]{summary.gam}}) of \code{object} is added to the
#'report. A plot of the spike train after time transformation
#'\code{\link{transformedTrain}} comes next followed by a renewal test plot
#'(\code{\link{renewalTestPlot}}) of the spike train on the time transformed
#'scale. The "usual" Ogata's tests plots (\code{\link{plot.transformedTrain}})
#'are added. Then if other trains are provided as a named list via argument
#'\code{neuronEvts}, interactions plots (\code{\link{plot.frt}}) are built
#'showing both the survivor function and the Berman's test. The report ends
#'with the \code{call} which generated \code{object}.
#'
#'@param object an object returned by \code{\link[mgcv]{gam}}.
#'@param filename a character string. The generic name of all the files (html,
#'png as well as \code{R} data files which will be generated. See also
#'\code{\link[R2HTML]{HTMLInitFile}}.
#'@param extension see \code{\link[R2HTML]{HTMLInitFile}}.
#'@param directory the full or relative path to the directory where the results
#'are going to be stored. See also \code{\link[R2HTML]{HTMLInitFile}}.
#'@param Title See \code{\link[R2HTML]{HTMLInitFile}}. If missing a default
#'value baed on \code{filename} is provided.
#'@param neuron a character string describing to which the analysis refers and
#'used for the titles of the interaction plots (see \code{\link{plot.frt}}).
#'@param neuronEvts a named list with the \code{event} variable from the data
#'frame returned by \code{\link{mkGLMdf}} and corresponding to the other
#'neurons recorded simultaneously. One list element per neuron.
#'@param \dots Not used, only there for compatibilty with the generic method
#'definition.
#'@return Nothing is returned, an html file and figures in png format are
#'written to disk.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{mkGLMdf}}, \code{\link[mgcv]{gam}},
#'\code{\link[mgcv]{gam.check}}, \code{\link{frt}},
#'\code{\link{transformedTrain}}, \code{\link{plot.transformedTrain}},
#'\code{\link{summary.transformedTrain}}
#'@keywords models smooth regression
#'@examples
#'
#'\dontrun{
#'## load e070528spont data set
#'data(e070528spont)
#'## make a data frame for gam using a 2 ms bin width 
#'spontDF <- mkGLMdf(e070528spont,0.002,0,60)
#'## make data frames specific of each neuron
#'n1.spontDF <- spontDF[spontDF$neuron=="1",]
#'n2.spontDF <- spontDF[spontDF$neuron=="2",]
#'n3.spontDF <- spontDF[spontDF$neuron=="3",]
#'n4.spontDF <- spontDF[spontDF$neuron=="4",]
#'## save space by removing the now redundant spontDF
#'rm(spontDF)
#'## fit neuron 1 using the gam representation of a
#'## renewal process and a binomial model
#'n1.spontFit1 <- gam(event ~ s(lN.1,k=25,bs="cr"),data=n1.spontDF,family=binomial())
#'## create a list with the discretized spike times of the 3 other neurons
#'preN1 <- list(n2=with(n2.spontDF,event),n3=with(n3.spontDF,event),n4=with(n4.spontDF,event))
#'## generate the report
#'reportHTML(n1.spontFit1,"e070528spontN1gFit",neuron="1",neuronEvts=preN1)
#'}
#'
reportHTML.gam <- function(object,
                           filename,
                           extension="html",
                           directory=getwd(),
                           Title,
                           neuron,
                           neuronEvts,
                           ...) {
  
  objectN <- deparse(substitute(object))
  
  ## if neuronEvts given check that it is a named list
  if (!missing(neuronEvts)) {
    if (!inherits(neuronEvts,"list")) neuronEvts <- list(neuronEvts)
    if (is.null(names(neuronEvts)))
      stop("neuronEvts should be a named list.")
  } ## End of conditional on !missing(neuronEvts)

  ## if neuron is missing try to figure out its value from the
  ## neuron variable of the "data" component of object
  if (missing(neuron)) {
    if (!("neuron" %in% names(object$data)))
      stop("neuron is not specified and cannot be obtained from object.")
    neuron <- with(object$data,unique(neuron)[1])
  } ## End of conditional on missing(neuron)
  
  if (missing(filename))
    filename <- paste(objectN,"GAM analysis")
  
  if (missing(Title))
    Title <- filename

  
  HTMLInitFile(outdir=directory,
               filename=filename,
               extension=extension,
               Title=Title)

  fullName <- paste(directory,"/",
                    filename,".",
                    extension,sep="")

  ## Write object's summary
  HTML.title("GAM fit summary:",file=fullName,HR=3)
  objectS <- summary(object) 
  HTML(objectS,file=fullName)

  
  HTMLhr(file=fullName)

  ## add a spike train plot after time transformation
  HTML.title("Spike train after time transformation:",
             file=fullName,HR=3)

  Lambda <- transformedTrain(object)
  st <- Lambda
  class(st) <- "spikeTrain"

  stFigName <- paste(filename,"_TTst.png",sep="")
  figFname <- paste(directory,"/",stFigName,sep="")
  png(figFname,width=500,height=500)
  plot(st,
       xlab="Transformed time",
       main=""
       )
  dev.off()
  HTMLInsertGraph(stFigName,
                  file=fullName,
                  WidthHTML=500,
                  HeightHTML=500)

  ## add a renewal test plot after time transformation
  HTMLhr(file=fullName)
  HTML.title("Renewal test after time transformation:",
             file=fullName,HR=3)
  
  rtFigName <- paste(filename,"_TTrt.png",sep="")
  figFname <- paste(directory,"/",rtFigName,sep="")
  png(figFname,width=800,height=800)
  renewalTestPlot(Lambda)
  dev.off()
  HTMLInsertGraph(rtFigName,
                  file=fullName,
                  WidthHTML=800,
                  HeightHTML=800)

  ## add a Ogata tests after time transformation
  HTMLhr(file=fullName)
  HTML.title("Ogata's tests after time transformation:",
             file=fullName,HR=3)
  
  otFigName <- paste(filename,"_TTot.png",sep="")
  figFname <- paste(directory,"/",otFigName,sep="")
  png(figFname,width=800,height=800)
  plot(Lambda,
       which=c(1,2,4,5),
       ask=FALSE)
  dev.off()
  HTMLInsertGraph(otFigName,
                  file=fullName,
                  WidthHTML=800,
                  HeightHTML=800)

  ## Check if the Ogata tests are passed and if yes and if
  ## neurons is given add interaction plots
  LambdaS <- summary(Lambda)
  otTRUE <- LambdaS[[1]][2] &&
  LambdaS[[2]][2] &&
  (0.005 <= LambdaS[[3]][2]) &&
  (0.995 >= LambdaS[[3]][2])
  
  if (otTRUE && !missing(neuronEvts)) {
    HTMLhr(file=fullName)
    HTML.title("Interaction tests after time transformation:",
               file=fullName,HR=3)

    ## load other neurons data
    pre.Lambda <- lapply(neuronEvts,
                         function(l) transformedTrain(object,l)
                         )
    names(pre.Lambda) <- names(neuronEvts)
    
    itFigName <- paste(filename,"_TTit.png",sep="")
    figFname <- paste(directory,"/",itFigName,sep="")
    png(figFname,width=300*length(neuronEvts),height=1000)
    layout(matrix(1:(2*length(neuronEvts)),
                  nrow=2)
           )

    sapply(seq(neuronEvts),
           function(nIdx) {
             theFRT <- pre.Lambda[[nIdx]] %frt% Lambda
             plot(theFRT,
                  main=paste(names(neuronEvts)[nIdx],"/ N",neuron),
                  ask=FALSE,which=1)
             plot(theFRT,
                  main=paste(names(neuronEvts)[nIdx],"/ N",neuron),
                  ask=FALSE,which=2)
           }
           )
    dev.off()
    HTMLInsertGraph(itFigName,
                    file=fullName,
                    WidthHTML=300*length(neuronEvts),
                    HeightHTML=1000)
  } ## end of conditional on otTRUE && !missing(neuronEvts)

  HTMLbr(1,file=fullName)
  HTMLhr(file=fullName)
  
  ## End with the call used to generate the object
  HTML.title(paste("Call used to generate object ",
                   objectN, ":",sep=""),
             file=fullName,HR=3)
  thisC <- object$call
  HTML(thisC,file=fullName)
  HTMLbr(1,file=fullName)
  
  HTMLEndFile()

}
