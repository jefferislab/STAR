

#'Spike Train Analysis with R
#'
#'Functions to analyze neuronal spike trains
#'
#'\tabular{ll}{ Package: \tab STAR\cr Type: \tab Package\cr Version: \tab
#'0.2-1\cr Date: \tab 2008-11-02\cr Depends: \tab gss, survival, mgcv,
#'R2HTML\cr Suggests: \tab lattice, HiddenMarkov, snow, rstream\cr License:
#'\tab GPL version 2 or newer\cr URL: \tab
#'http://sites.google.com/site/spiketrainanalysiswithr\cr }
#'
#'@name STAR-package
#'@aliases STAR-package STAR
#'@docType package
#'@author Christophe Pouzat
#'
#'Maintainer: Christophe Pouzat <christophe.pouzat@@gmail.com>
#'@keywords package
NULL





#'Shallow Shocks (M >= 6.0) in OFF Tohoku Area for 1885-1980
#'
#'Earthquakes data used by Yosihiko Ogata in his 1988 JASA paper.
#'
#'Quakes 213 and 214 were given exactly the same dates in Ogata (1988). Quake
#'214 has here been delayed by 1 minute.
#'
#'@name ShallowShocks
#'@docType data
#'@format A \code{data.frame} with the following variables:
#'
#'\tabular{ll}{ year: \tab year of occurrence.\cr month: \tab month of
#'occurrence.\cr day: \tab day of occurrence.\cr hour: \tab hour of
#'occurrence.\cr minute: \tab minute of occurrence.\cr magnitude: \tab
#'magnitude on Richter's scale.\cr type: \tab type of earthquake: \code{main}
#'(shock), \code{foreshock}, \code{aftershock}; according to Utsu.\cr Date:
#'\tab date in days starting from January 1st 1885.\cr energy.sqrt: \tab square
#'root of the energy expressed in erg.\cr }
#'@references Ogata, Yosihiko (1988) Statistical Models for Earthquake
#'Occurrences and Residual Analysis for Point Processes. \emph{Journal of the
#'American Statistical Association} \bold{83}: 9-27.
#'@source Ogata (1988) Table 1, pp 14-15.
#'@keywords datasets
#'@examples
#'
#'data(ShallowShocks)
#'## Reproduce Fig. 2 of Ogata 1988
#'layout(matrix(1:3, nrow = 3))
#'plot(ShallowShocks$Date,
#'     cumsum(ShallowShocks$energy.sqrt) / 10^13,
#'     type ="l",
#'     xlab = "",
#'     ylab = "",
#'     main = "Cumulative square root of energy")
#'plot(ShallowShocks$Date,
#'     cumsum(1+numeric(dim(ShallowShocks)[1])),
#'     type ="l",
#'     xlab = "",
#'     ylab = "",
#'     main = "Cumulative number of shocks")
#'plot(ShallowShocks$Date,
#'     ShallowShocks$magnitude,
#'     type = "h",
#'     ylim = c(5,9),
#'     xlab = "Time (days)",
#'     ylab = "",
#'     main = "Magnitude vs Occurrence time")
#'
NULL





#'Spike Trains of several Cockroach Antennal Lobe Neurons Recorded from Six
#'Animals
#'
#'Four (\code{CAL1S} and \code{CAL1V}), three (\code{CAL2S} and \code{CAL2C}),
#'three (\code{e060517spont} and \code{e060517ionon}), three
#'(\code{e060817spont}, \code{e060817terpi}, \code{e060817citron} and
#'\code{e060817mix}), two (\code{e060824spont} and \code{e060824citral}) and
#'four (\code{e070528spont} and \code{e070528citronellal}) Cockroach
#'(\emph{Periplaneta americana}) antennal lobe neurons (putative projection
#'neurons) were recorded simultaneously and extracellularly during spontaneous
#'activity and odors (vanilin, citral, citronellal, terpineol, beta-ionon)
#'responses from six different animals. The data sets contain the sorted spike
#'trains of the neurons.
#'
#'Every \code{repeatedTrain} object of these data sets has an \code{attribute}
#'named \code{stimTimeCourse} containing the openng and closing times of the
#'odor delivery valve.
#'
#'The data were recorded from neighboring sites on a \emph{NeuroNexus}
#'(\url{http://neuronexustech.com/}) silicon probe. Sorting was done with
#'\code{SpikeOMatic} with superposition resolution which can AND DOES lead to
#'artifcats on cross-correlograms.
#'
#'@name cockroachAlData
#'@aliases CAL1S CAL1V CAL2S CAL2C e060517spont e060517ionon e060817spont
#'e060817terpi e060817citron e060817mix e060824spont e060824citral e070528spont
#'e070528citronellal
#'@docType data
#'@format \code{CAL1S} is a named list with 4 components (\code{"neuron 1"},
#'\code{"neuron 2"}, \code{"neuron 3"}, \code{"neuron 4"}). Each component
#'contains the spike train (ie, action potentials occurrence times) of one
#'neuron recorded during 30 s of spontaneous activity. \emph{Times are
#'expressed in seconds}.
#'
#'\code{CAL1V} is a named list with 4 components (\code{"neuron 1"},
#'\code{"neuron 2"}, \code{"neuron 3"}, \code{"neuron 4"}).  Each component is
#'a named list with 20 components: \code{"stim. 1"}, ..., \code{"stim. 20"}.
#'Each sub-list contains the spike train of one neuron during 1 stimulation
#'(odor puff) with \emph{vanillin}
#'(\url{http://en.wikipedia.org/wiki/Vanillin}). Each acquisition was 10 s
#'long. The command to the odor delivery valve was on between sec 4.49 and sec
#'4.99.
#'
#'\code{CAL2S} is a named list with 3 components (\code{"neuron 1"},
#'\code{"neuron 2"}, \code{"neuron 3"}). Each component contains the spike
#'train (ie, action potentials occurrence times) of one neuron recorded during
#'1 mn of spontaneous activity. \emph{Times are expressed in seconds}.
#'
#'\code{CAL2C} is a named list with 3 components (\code{"neuron 1"},
#'\code{"neuron 2"}, \code{"neuron 3"}).  Each component is a named list with
#'20 components: \code{"stim. 1"}, ..., \code{"stim. 20"}. Each sub-list
#'contains the spike train of one neuron during 1 stimulation (odor puff) with
#'\emph{citral} (\url{http://en.wikipedia.org/wiki/Citral}). Each acquisition
#'was 14 s long. The command to the odor delivery valve was on between sec 5.87
#'and sec 6.37.
#'
#'\code{e060517spont} is a named list of with 3 components (\code{"neuron 1"},
#'\code{"neuron 2"}, \code{"neuron 3"}). Each component is a \code{spikeTrain}
#'object (ie, action potentials occurrence times) of one neuron recorded during
#'61 s of spontaneous activity. \emph{Times are expressed in seconds}.
#'
#'\code{e060517ionon} is a named list with 3 components (\code{"neuron 1"},
#'\code{"neuron 2"}, \code{"neuron 3"}).  Each component is a
#'\code{repeatedTrain} object with 19 \code{spikeTrain} objects: \code{"stim.
#'1"}, ..., \code{"stim. 19"}. Each \code{spikeTrain} contains the spike train
#'of one neuron during 1 stimulation (odor puff) with \emph{beta-ionon}
#'(\url{http://commons.wikimedia.org/wiki/Image:Beta-Ionon.svg}). Each
#'acquisition was 15 s long. The command to the odor delivery valve was on
#'between sec 6.07 and sec 6.57.
#'
#'\code{e060817spont} is a named list of with 3 components (\code{"neuron 1"},
#'\code{"neuron 2"}, \code{"neuron 3"}). Each component is a \code{spikeTrain}
#'object (ie, action potentials occurrence times) of one neuron recorded during
#'60 s of spontaneous activity. \emph{Times are expressed in seconds}.
#'
#'\code{e060817terpi} is a named list with 3 components (\code{"neuron 1"},
#'\code{"neuron 2"}, \code{"neuron 3"}).  Each component is a
#'\code{repeatedTrain} object with 20 \code{spikeTrain} objects: \code{"stim.
#'1"}, ..., \code{"stim. 20"}. Each \code{spikeTrain} contains the spike train
#'of one neuron during 1 stimulation (odor puff) with \emph{terpineol}
#'(\url{http://en.wikipedia.org/wiki/Terpineol}). Each acquisition was 15 s
#'long. The command to the odor delivery valve was on between sec 6.03 and sec
#'6.53.
#'
#'\code{e060817citron} is a named list with 3 components (\code{"neuron 1"},
#'\code{"neuron 2"}, \code{"neuron 3"}).  Each component is a
#'\code{repeatedTrain} object with 20 \code{spikeTrain} objects: \code{"stim.
#'1"}, ..., \code{"stim. 20"}. Each \code{spikeTrain} contains the spike train
#'of one neuron during 1 stimulation (odor puff) with \emph{citronellal}
#'(\url{http://en.wikipedia.org/wiki/Citronellal}). Each acquisition was 15 s
#'long. The command to the odor delivery valve was on between sec 5.99 and sec
#'6.49.
#'
#'\code{e060817mix} is a named list with 3 components (\code{"neuron 1"},
#'\code{"neuron 2"}, \code{"neuron 3"}).  Each component is a
#'\code{repeatedTrain} object with 20 \code{spikeTrain} objects: \code{"stim.
#'1"}, ..., \code{"stim. 20"}. Each \code{spikeTrain} contains the spike train
#'of one neuron during 1 stimulation (odor puff) with a mixture of
#'\emph{terpinaol} and \emph{citronellal} (the sum of the two previous stim.).
#'Each acquisition was 15 s long. The command to the odor delivery valve was on
#'between sec 6.01 and sec 6.51.
#'
#'\code{e060824spont} is a named list of with 2 components (\code{"neuron 1"},
#'\code{"neuron 2"}). Each component is a \code{spikeTrain} object (ie, action
#'potentials occurrence times) of one neuron recorded during 59 s of
#'spontaneous activity. \emph{Times are expressed in seconds}.
#'
#'\code{e060824citral} is a named list with 2 components (\code{"neuron 1"},
#'\code{"neuron 2"}).  Each component is a named list with 20 components:
#'\code{"stim. 1"}, ..., \code{"stim. 20"}. Each sub-list contains the spike
#'train of one neuron during 1 stimulation (odor puff) with \emph{citral}
#'(\url{http://en.wikipedia.org/wiki/Citral}). Each acquisition was 15 s long.
#'The command to the odor delivery valve was on between sec 6.01 and sec 6.51.
#'
#'\code{e070528spont} is a named list of with 4 components (\code{"neuron 1"},
#'\code{"neuron 2"}, \code{"neuron 3"}, \code{"neuron 4"}). Each component is a
#'\code{spikeTrain} object (ie, action potentials occurrence times) of one
#'neuron recorded during 60 s of spontaneous activity. \emph{Times are
#'expressed in seconds}.
#'
#'\code{e070528citronellal} is a named list with 4 components (\code{"neuron
#'1"}, \code{"neuron 2"}, \code{"neuron 3"}, \code{"neuron 4"}).  Each
#'component is a \code{repeatedTrain} object with 15 \code{spikeTrain} objects:
#'\code{"stim. 1"}, ..., \code{"stim. 15"}. Each \code{spikeTrain} contains the
#'spike train of one neuron during 1 stimulation (odor puff) with
#'\emph{citronellal} (\url{http://en.wikipedia.org/wiki/Citronellal}). Each
#'acquisition was 13 s long. The command to the odor delivery valve was on
#'between sec 6.14 and sec 6.64.
#'@references
#'\url{http://www.biomedicale.univ-paris5.fr/physcerv/C_Pouzat/Doc/ChaffiolEtAl_FENS2006.pdf}
#'@source Recording and spike sorting performed by Antoine Chaffiol
#'\email{antoine.chaffiol@@univ-paris5.fr} at the Cerebral Physiology Lab, CNRS
#'UMR 8118:
#'\url{http://www.biomedicale.univ-paris5.fr/physcerv/physiologie_cerebrale.htm}.
#'@keywords datasets
#'@examples
#'
#'## load CAL1S data
#'data(CAL1S)
#'## convert the data into spikeTrain objects
#'CAL1S <- lapply(CAL1S,as.spikeTrain)
#'## look at the train of the 1sd neuron
#'CAL1S[["neuron 1"]]
#'## fit the 6 different renewal models to the 1st neuron spike train
#'compModels(CAL1S[["neuron 1"]])
#'## look at the ISI distribution with the fitted invgauss dist for
#'## this 1st neuron
#'isiHistFit(CAL1S[["neuron 1"]],model="invgauss")
#'
#'## load CAL1V data
#'data(CAL1V)
#'## convert them to repeatedTrain objects
#'CAL1V <- lapply(CAL1V, as.repeatedTrain)
#'## look at the raster of the 1st neuron
#'CAL1V[["neuron 1"]]
#'
#'## load e070528spont data
#'data(e070528spont)
#'## look at the spike train of the 1st neuron
#'e070528spont[["neuron 1"]]
#'
#'## load e070528citronellal data
#'data(e070528citronellal)
#'## Get the stimulus time course
#'attr(e070528citronellal[["neuron 1"]],"stimTimeCourse")
#'## look at the raster of the 1st neuron
#'plot(e070528citronellal[["neuron 1"]],stim=c(6.14,6.64))
#'
#'\dontrun{
#'## A "detailed" analysis of e060817 were 2 odors as well as there mixtures
#'## were used.
#'## Load the terpineol, citronellal and mixture response data
#'data(e060817terpi)
#'data(e060817citron)
#'data(e060817mix)
#'## get smooth psths with gsspsth0
#'e060817terpiN1PSTH <- gsspsth0(e060817terpi[["neuron 1"]])
#'e060817terpiN2PSTH <- gsspsth0(e060817terpi[["neuron 2"]])
#'e060817terpiN3PSTH <- gsspsth0(e060817terpi[["neuron 3"]])
#'e060817citronN1PSTH <- gsspsth0(e060817citron[["neuron 1"]])
#'e060817citronN2PSTH <- gsspsth0(e060817citron[["neuron 2"]])
#'e060817citronN3PSTH <- gsspsth0(e060817citron[["neuron 3"]])
#'e060817mixN1PSTH <- gsspsth0(e060817mix[["neuron 1"]])
#'e060817mixN2PSTH <- gsspsth0(e060817mix[["neuron 2"]])
#'e060817mixN3PSTH <- gsspsth0(e060817mix[["neuron 3"]])
#'## look at them
#'## Neuron 1
#'plot(e060817terpiN1PSTH,stimTimeCourse=attr(e060817terpi[["neuron 1"]],"stimTimeCourse"),colCI=2)
#'plot(e060817citronN1PSTH,stimTimeCourse=attr(e060817citron[["neuron 1"]],"stimTimeCourse"),colCI=2)
#'plot(e060817mixN1PSTH,stimTimeCourse=attr(e060817mix[["neuron 1"]],"stimTimeCourse"),colCI=2)
#'## Neuron 2
#'plot(e060817terpiN2PSTH,stimTimeCourse=attr(e060817terpi[["neuron 2"]],"stimTimeCourse"),colCI=2)
#'plot(e060817citronN2PSTH,stimTimeCourse=attr(e060817citron[["neuron 2"]],"stimTimeCourse"),colCI=2)
#'plot(e060817mixN2PSTH,stimTimeCourse=attr(e060817mix[["neuron 2"]],"stimTimeCourse"),colCI=2)
#'## Neuron 3
#'plot(e060817terpiN3PSTH,stimTimeCourse=attr(e060817terpi[["neuron 3"]],"stimTimeCourse"),colCI=2)
#'plot(e060817citronN3PSTH,stimTimeCourse=attr(e060817citron[["neuron 3"]],"stimTimeCourse"),colCI=2)
#'plot(e060817mixN3PSTH,stimTimeCourse=attr(e060817mix[["neuron 3"]],"stimTimeCourse"),colCI=2)
#'
#'## Make now fancier plots with superposed psths ####
#'## Take into account the fact that the stimuli onsets are not identical
#'
#'## Neuron 1
#'plot(e060817mixN1PSTH$mids-0.02,e060817mixN1PSTH$ciUp,type="n",ylim=c(0,max(e060817mixN1PSTH$ciUp)),xlim=c(5,14),xlab="Time (s)",ylab="Firing rate (Hz)",main="Neuron 1 e060817")
#'rect(5.99,0,6.49,max(e060817mixN1PSTH$ciUp),col="grey80",border=NA)
#'abline(h=0)
#'polygon(c(e060817mixN1PSTH$mids-0.02,rev(e060817mixN1PSTH$mids-0.02)),c(e060817mixN1PSTH$ciLow,rev(e060817mixN1PSTH$ciUp)),col=rgb(1,0,1,0.5),border=NA)
#'polygon(c(e060817citronN1PSTH$mids,rev(e060817citronN1PSTH$mids)),c(e060817citronN1PSTH$ciLow,rev(e060817citronN1PSTH$ciUp)),col=rgb(1,0,0,0.5),border=NA)
#'polygon(c(e060817terpiN1PSTH$mids-0.04,rev(e060817terpiN1PSTH$mids-0.04)),c(e060817terpiN1PSTH$ciLow,rev(e060817terpiN1PSTH$ciUp)),col=rgb(0,0,1,0.5),border=NA)
#'lines(e060817terpiN1PSTH$mids-0.04,e060817terpiN1PSTH$freq,col=rgb(0,0,1),lwd=2)
#'lines(e060817citronN1PSTH$mids,e060817citronN1PSTH$freq,col=rgb(1,0,0),lwd=2)
#'lines(e060817mixN1PSTH$mids-0.02,e060817mixN1PSTH$freq,col=rgb(0,0,0),lwd=2)
#'legend(8,0.9*max(e060817mixN1PSTH$ciUp),c("Terpineol","Citronellal","Mixture"),col=c(4,2,1),lwd=2)
#'
#'## Neuron 2
#'plot(e060817mixN2PSTH$mids-0.02,e060817mixN2PSTH$ciUp,type="n",ylim=c(0,max(e060817mixN2PSTH$ciUp)),xlim=c(5,14),xlab="Time (s)",ylab="Firing rate (Hz)",main="Neuron 2 e060817")
#'rect(5.99,0,6.49,max(e060817mixN2PSTH$ciUp),col="grey80",border=NA)
#'abline(h=0)
#'polygon(c(e060817mixN2PSTH$mids-0.02,rev(e060817mixN2PSTH$mids-0.02)),c(e060817mixN2PSTH$ciLow,rev(e060817mixN2PSTH$ciUp)),col=rgb(1,0,1,0.5),border=NA)
#'polygon(c(e060817citronN2PSTH$mids,rev(e060817citronN2PSTH$mids)),c(e060817citronN2PSTH$ciLow,rev(e060817citronN2PSTH$ciUp)),col=rgb(1,0,0,0.5),border=NA)
#'polygon(c(e060817terpiN2PSTH$mids-0.04,rev(e060817terpiN2PSTH$mids-0.04)),c(e060817terpiN2PSTH$ciLow,rev(e060817terpiN2PSTH$ciUp)),col=rgb(0,0,1,0.5),border=NA)
#'lines(e060817terpiN2PSTH$mids-0.04,e060817terpiN2PSTH$freq,col=rgb(0,0,1),lwd=2)
#'lines(e060817citronN2PSTH$mids,e060817citronN2PSTH$freq,col=rgb(1,0,0),lwd=2)
#'lines(e060817mixN2PSTH$mids-0.02,e060817mixN2PSTH$freq,col=rgb(0,0,0),lwd=2)
#'legend(8,0.9*max(e060817mixN2PSTH$ciUp),c("Terpineol","Citronellal","Mixture"),col=c(4,2,1),lwd=2)
#'
#'## Neuron 3
#'plot(e060817mixN3PSTH$mids-0.02,e060817mixN3PSTH$ciUp,type="n",ylim=c(0,max(e060817mixN3PSTH$ciUp)),xlim=c(5,14),xlab="Time (s)",ylab="Firing rate (Hz)",main="Neuron 3 e060817")
#'rect(5.99,0,6.49,max(e060817mixN3PSTH$ciUp),col="grey80",border=NA)
#'abline(h=0)
#'polygon(c(e060817mixN3PSTH$mids-0.02,rev(e060817mixN3PSTH$mids-0.02)),c(e060817mixN3PSTH$ciLow,rev(e060817mixN3PSTH$ciUp)),col=rgb(1,0,1,0.5),border=NA)
#'polygon(c(e060817citronN3PSTH$mids,rev(e060817citronN3PSTH$mids)),c(e060817citronN3PSTH$ciLow,rev(e060817citronN3PSTH$ciUp)),col=rgb(1,0,0,0.5),border=NA)
#'polygon(c(e060817terpiN3PSTH$mids-0.04,rev(e060817terpiN3PSTH$mids-0.04)),c(e060817terpiN3PSTH$ciLow,rev(e060817terpiN3PSTH$ciUp)),col=rgb(0,0,1,0.5),border=NA)
#'lines(e060817terpiN3PSTH$mids-0.04,e060817terpiN3PSTH$freq,col=rgb(0,0,1),lwd=2)
#'lines(e060817citronN3PSTH$mids,e060817citronN3PSTH$freq,col=rgb(1,0,0),lwd=2)
#'lines(e060817mixN3PSTH$mids-0.02,e060817mixN3PSTH$freq,col=rgb(0,0,0),lwd=2)
#'legend(8,0.9*max(e060817mixN3PSTH$ciUp),c("Terpineol","Citronellal","Mixture"),col=c(4,2,1),lwd=2)
#'}
#'
#'
NULL





#'Spike Trains of a Purkinje Cells (PC) Recorded in Control Conditions and With
#'Bath Applied Bicuculline
#'
#'An object of class \code{"SpikeTrain"}. Spontaneous discharge of a single PC
#'recorded during 300 s in normal saline conditions and during 300 s in the
#'presence of 25 \eqn{\mu}{micro}M bath applied bicuculline.
#'
#'The recording contained in \code{sPK} was done in cell-attached mode. The one
#'in \code{mPK} was done with a NeuroNexus silicon probe.
#'
#'Bicuculline is a GABAA receptor antagonist. It blocks all GABAA inhibition.
#'
#'@name purkinjeCellData
#'@aliases sPK mPK
#'@docType data
#'@format \code{sPK} is a named list with 2 components (\code{"ctl"},
#'\code{"bicu"}. Each component contains the spike train (ie, action potentials
#'occurrence times) of one Purkinje cell recorded during 300 s of spontaneous
#'activity in control (\code{"ctl"}) condition and with bath applied
#'bicuculline (\code{"bicu"}). \emph{Times are expressed in seconds}.
#'
#'\code{mPK} is a named list with 8 components (\code{"neuron 1"},
#'\code{"neuron 2"}, ..., \code{"neuron 8"}. Each component is itself a list
#'with the spike train (ie, action potentials occurrence times) of one Purkinje
#'cell recorded during 300 s of spontaneous activity in control (\code{"ctl"})
#'condition and with bath applied bicuculline (\code{"bicu"}). \emph{Times are
#'expressed in seconds}.
#'@source Recording and spike sorting performed by Matthieu Delescluse at the
#'Cerebral Physiology Lab, CNRS UMR 8118:
#'\url{http://www.biomedicale.univ-paris5.fr/physcerv/physiologie_cerebrale.htm}.
#'@keywords datasets
#'@examples
#'
#'\dontrun{
#'## load spontaneous data of 1 Purkinje cell
#'## recorded in cell attached mode from a cerebellar
#'## slice in control and bath applied bicuculline conditions
#'data(sPK)
#'## coerce data to spikeTrain objects
#'sPK <- lapply(sPK,as.spikeTrain)
#'## Get a summary of the ctl data
#'summary(sPK[["ctl"]])
#'## Look at the control train
#'## Don't show the rug plot for clarity
#'plot(sPK[["ctl"]],addRug=FALSE)
#'## Generate the renewal test plot taking into account
#'## the size of the data set (a lot of spikes!).
#'renewalTestPlot(sPK[["ctl"]],d=10,orderPlotPch=".",lag.max=250)
#'## Get a summary of the bicu data
#'summary(sPK[["bicu"]])
#'## Look at the control train
#'## Don't show the rug plot for clarity
#'plot(sPK[["bicu"]],addRug=FALSE)
#'## Generate the renewal test plot taking into account
#'## the size of the data set (a lot of spikes!).
#'renewalTestPlot(sPK[["bicu"]],d=10,orderPlotPch=".",lag.max=250);par(oldpar)
#'## This time the data are NOT stationary. This is seen clearly on a acf
#'## plot with very large lag.max
#'acf.spikeTrain(sPK[["bicu"]],lag.max=2000)
#'
#'}
#'
NULL



