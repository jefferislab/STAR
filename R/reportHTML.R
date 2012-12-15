#'Generic Function for Automatic HTML Report Generation
#'
#'When a standard analysis is applied to some object it is useful to keep all
#'the plots and summaries related to that analysis in a single place where they
#'can be easily accessed and visualized. An html file containing the report of
#'this analysis is ideally suited for that. The methods \code{reportHTML}
#'generate such reports.
#'
#'
#'@param object an object from which the report is going to be generated,
#'perhaps following some standard analysis procedure.
#'@param filename a character string. The generic name of all the files (html,
#'png as well as \code{R} data files which will be generated. See also
#'\code{\link[R2HTML]{HTMLInitFile}}.
#'@param extension see \code{\link[R2HTML]{HTMLInitFile}}.
#'@param directory the full or relative path to the directory where the results
#'are going to be stored. See also \code{\link[R2HTML]{HTMLInitFile}}.
#'@param Title See \code{\link[R2HTML]{HTMLInitFile}}. If missing a default
#'value baed on \code{filename} is provided.
#'@param \dots additional parameters passed to the functions internally called
#'by the actual methods.
#'@return Nothing is returned, an html file and figures in png format are
#'written to disk together with the \code{R} variables generated during the
#'analysis , if an analysis was performed.
#'@author Christophe Pouzat \email{christophe.pouzat@@gmail.com}
#'@seealso \code{\link{reportHTML.spikeTrain}},
#'\code{\link{reportHTML.repeatedTrain}}, \code{\link{reportHTML.gam}}
#'@keywords print file IO
#'@examples
#'
#'##
#'
reportHTML <- function(object,
                       filename,
                       extension,
                       directory,
                       Title,
                       ...)
### Generic reportHTML method definition
### object should contain data like a spikeTrain object
### or a repeatedTrain object on whih some standard analysis
### should be carried out first before exporting the results,
### that is, numerical and graphical summaries as well as,
### if necessary actual R objects onto the disk. The summaries
### are formated in html for easy inspection.
### Arguments filename, extension, directory and Title are passed
### to HTMLInitFile function (from R2HTML), directory corresponds
### to argument outdir of the latter.
{

  UseMethod("reportHTML")

}
