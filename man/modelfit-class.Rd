\name{modelfit-class}
\docType{class}
\alias{modelfit-class}
\alias{modelfit}
\alias{show,modelfit-method}
\alias{visualize,modelfit,missing,missing-method}
\title{"Modelfit"}

\description{Object returned from a call to \link{fitModelParameters}.}

\section{Slots}{
  \describe{
    \item{model}{The chosen model, either \code{"Gaussian"} or \code{"EMG"}.}
    \item{fitting}{The mode of fitting which generated the object. Either
      \code{"most_intense"} or \code{"model"}.}
    \item{alphafunction}{A function of m/z that computes an estimate for the
  \code{EMG} parameter \code{alpha}. See \link{fitModelParameters} for
  more detailed information on this parameter.}
    \item{sigmafunction}{A function of m/z that computes an estimate for the
  \code{EMG} parameter \code{sigma}. See \link{fitModelParameters} for
  more detailed information on this parameter.}
    \item{mufunction}{A function of m/z that computes an estimate for the
  \code{EMG} parameter \code{mu}. See \link{fitModelParameters} for
  more detailed information on this parameter.}
    \item{peakfitresults}{A \code{matrix} of five columns if \code{model =
	"Gaussian"} and of six columns if \code{model = "EMG"},
  respectively. The matrix contains basic information on the extracted
  peaks and the corresponding parameters. The first column contains the
  number of (mz, intensity)-pairs used for fitting, the second column
  contains the residual sums of squares of the model fit, the middle
  columns the parameter estimates and the last column the m/z position.}
    \item{bestpeak}{A \code{list} containing detailed information on the
  'best' peak, where 'best' is equivalent to 'minimum residual sums of
  squares' after fitting. Mainly used indirectly by calling the function \code{visualize}.}
}
}


\section{Methods}{
  \describe{    
    \item{show}{Use \code{show(object)} for brief information about the object.}
    \item{visualize}{A function to display graphically the result of
  parameter estimation. One can either visualize a single peak
  (\code{slot(object, "bestpeak")}) or the fit of the linear models
  postulated for the parameters \code{alpha}, \code{sigma} and
  \code{mu}. The function \code{visualize} is called with the following
  set of parameters.
  \describe{
    \item{\code{object}}{An object of class \code{modelfit}.}
    \item{\code{type}}{A \code{character} specifying the object to be
      visualized. If \code{type = "peak"},
     the fit of a single peak stored in \code{object@bestpeak} is
  displayed. If \code{type = "model"}, one obtains scatterplots of the
  form \code{parameter} vs. \code{mz} for each parameter in
  \code{parameters}, s. below.}
   \item{\code{parameters}}{Needed if and only if \code{type = "model"} in
  order to choose the y-variable of the scatterplot. Several parameters
  may be specified at a time, in which case one obtains a multi-panel
  plot.}
\item{\code{modelfit}}{A \code{logical} indicating whether the estimated
  regression functions should be added to the
  scatterplots (if \code{type = "model"}).
 }
}
}
}
} 

\keyword{models}

 
 
