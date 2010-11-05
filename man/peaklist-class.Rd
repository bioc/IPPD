\name{peaklist-class}
\docType{class}
\alias{peaklist-class}
\alias{peaklist}
\alias{show,peaklist-method}
\alias{threshold}
\alias{threshold,peaklist-method}
\alias{visualize}
\alias{visualize,peaklist,numeric,numeric-method}
\title{"Peaklist"}

\description{Object returned from a call to \link{getPeaklist}.}

\section{Slots}{
  \describe{
    \item{peaklist}{A \code{matrix} with rows equal to the number of
  peak patterns used in fitting and four
  columns. For each row (pattern), the first column
  contains the initial position of the peak pattern (\code{loc_init}, in
  Thomson), the second one
  contains the position of the most intense peak within a peak pattern
  (\code{loc_most_intense}, in Thomson), the third one contains the
  charge state (\code{charge}), and the last 
  one contains the intensity (\code{amplitude}).}
\item{peaklistprocessed}{The result after applying the postprocessing
  procedure to the raw list \code{peaklist}. It is a matrix of the same structure as
  \code{peaklist}, with two up to four additional columns: the fifth column contains the local noise level at the
  position, and the last one contains the ratio of intensity to local
  noise level (signal-to-noise ratio), hence quantitatively describing the significance of the
  peak pattern. If \code{goodnessoffit} has been set to \code{TRUE},
  there are two additional columns labelled \code{goodness_of_fit} and
  \code{ratio_adj}. The column \code{goodness_of_fit} contains a local assessment of the goodness-of-fit
  at the respective positions in the spectrum, while the column
  \code{ratio_adj} contains a goodness-of-fit adjusted signal-to-noise ratio. 
  If \code{postprocessed} has been set to \code{FALSE}
  when calling \link{getPeaklist}, this slot is an empty matrix.}
\item{model}{One of "Gaussian" or "EMG".}
\item{loss}{One of "L2" or "L1".}
\item{alpha}{Function for the parameter "alpha" (\code{model = "EMG"}).}
\item{sigma}{Function for the parameter "sigma" (both models).}
\item{mu}{Function for the parameter "mu" (\code{model = "EMG"}).}
\item{charges}{A \code{numeric} vector containing the set of charge states used.}
\item{basis}{The matrix of basis functions if \code{returnbasis} has
  been set to \code{TRUE} when calling \link{getPeaklist}. Otherwise,
  this slot is an empty matrix.}
\item{book}{A \code{matrix} annotating each column of \code{basis}. Each
  row corresponds to a template/peak pattern characterized by the
  initial position of the peak pattern (column \code{initial}),
  the position of the most intense peak (column \code{most_intense}) and the charge
  state (column \code{charge}).}
\item{beta}{A \code{numeric} vector storing the intensities for each
  peak pattern in \code{book}/each column of \code{basis}, obtained by
  non-negative least squares estimation (if \code{loss = "L2"}) or by
  non-negative least absolute deviation estimation (if \code{loss = "L1"}). Hence, together with \code{book},
  \code{cbind(book, beta)} constitutes a predecessor of \code{peaklist}.}
\item{locnoise}{A \code{matrix} of local quantiles of the intensities
  computed along the spectrum. Each row corresponds to a m/z value of
  the spectrum (ordered increasingly), and each column to a quantile as
  indicated by the column names. For instance, \code{slot(object,
    "locnoise")[,"0.5"]} would give the vector of local medians.}
\item{noiselevel}{A \code{numeric} vector storing the noise level used for computation of the
  signal-to-noise ratios appearing in \code{peaklistprocessed}.}
\item{goodnessoffit}{A \code{numeric} vector storing the local goodness-of-fit criterion. An empty
  vector unless \code{control.postprocessing$goodnessfittrue} was
  \code{TRUE} when \code{getPeaklist} was called.}
\item{data}{A \code{list} storing the spectrum used to generate the object.}
}
}


\section{Methods}{
  \describe{    
    \item{show}{Use \code{show(object)} for brief information about the
      object.}
    \item{threshold}{Only applicable if \code{postprocessed} has been
  set to \code{TRUE} when calling \code{getPeaklist}. A function to threshold \code{peaklistprocessed},
  yielding a peaklist maintaining only those peaks exceeding a specified
  threshold. Optionally, re-fitting of the spectrum using only those
  templates that have passed the thresholding step can be performed. The
  argument list is given by
  \describe{\item{\code{object}}{An object of class \code{peaklist}.}
            \item{\code{threshold}}{The value to be used as threshold.}
            \item{\code{ratio}}{Whether to use the signal-to-noise ratio
  \code{"ratio"} or the goodness-of-fit adjusted signal-to-noise ratio
  \code{"ratio_adj"}}
 \item{\code{refit}}{Whether re-fitting of the spectrum should be done,
   using only the templates above the threshold.}
 \item{\code{trace}}{A \code{logical} indicating whether information tracing the different steps of the
   fitting process should be displayed. Only used if \code{refit = TRUE}.}
  \item{\code{eps}}{Function values below \code{eps} are set equal to
      precisely zero in order to make the basis function matrix sparse. Only used if \code{refit = TRUE}.}.}}

    \item{visualize}{A function to displaying the results of the fitting
  procedure, in particular significant peak patterns, in selected m/z
  ranges. The description of the arguments is as follows:
  \describe{
         \item{\code{object}}{An object of class \code{peaklist}.}
	 \item{\code{mz}}{The vector of m/z values used to obtain \code{object}.}
         \item{\code{intensities}}{The vector of intensities used to obtain \code{object}.}
	 \item{\code{lower, upper}}{The m/z range for which the result of the
	   fitting procedure should be visualized (\code{lower <
	     upper}). \emph{Hint}: the difference \code{upper - lower} should be
           between 0-30. If it is significantly larger, computation
  times increase considerably. Moreover, this will lead to a reduced
  quality of the plot.}
	 \item{\code{truth}}{May be set to \code{TRUE} if one has a precise
  knowledge of the true underlying signal (normally only occuring for
  simulated data) and one wants to compare the results with the gold
  standard. If \code{truth = TRUE}, \code{booktrue} (s. below) has to be
  specified. In this case, the true signal is displayed as upper panel
  in the plot.}
\item{\code{signal}}{A \code{logical} indicating whether the raw spectrum in the range \code{[lower,upper]}
  should be plotted. If \code{TRUE}, it is displayed as upper panel or as one
  of the middle panels.}
\item{\code{fitted}}{A \code{logical} indicating whether basically all (with restrictions,
  s. \code{cutoff.functions} and \code{cutoff.eps} below) fitted
  templates should be plotted. If \code{TRUE}, the result is displayed as a middle
  panel or as lower panel.}
\item{\code{postprocessed}}{A \code{logical} indicating whether fitted
  templates after postprocessing should be
  plotted. If \code{TRUE}, the result is displayed as lower panel.}
\item{\code{fittedfunction}}{A \code{logical} indicating whether the
  fitted function should be drawn as lines in
  the panel displaying the raw spectrum (requires \code{signal =
  TRUE}). Note that the fitted function is the result of summing up all
fitted templates multiplied by their intensity. Setting
  \code{fittedfunction = TRUE} requires that \code{object} contains the
  matrix of basis functions.}
\item{\code{fittedfunction.cut}}{A \code{logical} indicating whether the fitted function after cutting peak
  patterns falling below a prespecified cutoff should be drawn as lines in the panel displaying the raw spectrum (requires \code{signal =
  TRUE}). Setting \code{fittedfunction.cut = TRUE} requires that \code{object} contains the
matrix of basis functions.}
\item{\code{quantile}}{Optional argument. If \code{quantile} is
  specified, it has to be one of the values \code{0.1,...,0.9}, and the
  corresponding local quantile as contained in the slot \code{locnoise}
  (s. above) will be added to the plot as dotted line.}
\item{\code{booktrue}}{Required if and only if \code{truth = TRUE}. In this
  case, \code{booktrue} has to be a matrix having the same structure as
  \code{cbind(slot(object, "book"), slot(object, "beta"))}, i.e. the first
  two columns of \code{booktrue} contain initial- and most intense peak
  locations within a peak pattern, the third column contains the charge
  state and the fourth one the intensity. \emph{Note:} the 'truth' is computed
  according to \code{slot(object, model)}. This excludes the following
  scenario: assuming that the true peak shapes is \code{"EMG"}, but the
  chosen  model is \code{"Gaussian"}. Then the true model in this function is
  computed according to the latter model.}
\item{\code{cutoff.eps, cutoff.functions}}{Control arguments reducing the
  number of templates to be plotted if \code{fitted = TRUE}. Specifying
  \code{cutoff.eps} removes all templates with intensity smaller than
  \code{cutoff.eps}. Specifying \code{cutoff.functions} as a positive
  integer only plots the \code{cutoff} 'most intense' templates.}
 \item{\code{...}}{Additional options passed to \code{plot}.}	
      }
      Note that the colours in the plot have no meaning. They
  are only used to distinguish between different patterns within a
  panel. Further note that the colour of a pattern appearing in several panels
  may vary from panel to panel. 
   }
 }
} 

\keyword{models}

 
