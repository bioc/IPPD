\name{getPeaklist}
\alias{getPeaklist,numeric,numeric-method}
\alias{getPeaklist}
\alias{getPeaklist-methods}
\title{Peak pattern extraction by non-negative ls/lad template matching}
\description{
  Generates a candidate list of isotopic peak patterns
  present in a protein mass spectrum. This is achieved by matching
  templates calculated according to the so-called Averagine model to the raw
  spectrum using either non-negative least squares (ls) or non-negative
  least absolute deviation (lad) estimation. The presence of multiple charge
  states is supported. In particular, the approach is capable of
  deconvolving overlapping patterns. The method can be applied with two
  different kind of peak shapes, Gaussians and Exponentially Modified Gaussians (EMG).
}

\usage{
getPeaklist(mz, intensities, model = c("Gaussian", "EMG"),
model.parameters = list(alpha = function(){},
                        sigma = function(){},
                        mu = function(){}),
loss = c("L2", "L1"), binning = FALSE,
postprocessing = TRUE, trace = TRUE, returnbasis = TRUE,
control.basis = list(charges = c(1,2,3,4),
                     eps = 1e-05),
control.localnoise = list(quantile = 0.5,
                          factor.place = 1.5,
                          factor.post = 0,
                          window = NULL,
                          subtract = FALSE),
control.postprocessing = list(mzfilter = FALSE,
                              prune = FALSE,
                              factor.prune = NULL,
                              ppm = NULL,
                              goodnessoffit = FALSE),
control.binning = list(tol = 0.01))
}
\arguments{
   \item{mz}{A \code{numeric} vector of m/z (mass/charge) values (in
     Thomson), ordered increasingly.}
   \item{intensities}{A \code{numeric} vector of intensities
    corresponding to \code{mz}.}
  \item{model}{Basic model for the shape of a single peak. Must be
    \code{"Gaussian"} or \code{"EMG"} (exponentially modified
    Gaussian). See \link{fitModelParameters} for further information on
    these models.}
  \item{loss}{The loss function to be used. The choice \code{loss =
      "L2"} yield a nonnegative least squares fit, \code{loss = "L1"} a
    nonnegative least absolute deviation fit. The second choice is
    more robust when deviations from model assumptions (peak model,
    Averagine model,\code{...}) frequently occur in the data. Note,
    however, that computation time is much higher for the second choice
    (at least by a factor two).}
  \item{model.parameters}{A list of functions with precisely one
    argument representing \code{mz}. The parameters of a single peak are
    typically modeled as a function of m/z. If \code{model = "Gaussian"},
    the peak shape depends on the parameter \code{sigma} (a function
    \code{sigma(mz)}). If \code{model = "EMG"}, the peak shape
    additionally depends on two parameters \code{alpha} and \code{mu}
    (two functions \code{alpha(mz)} and \code{mu(mz)}). Note that
    constant functions are usually specified by using
    a construction of the form \code{parameter(mz) <- function(mz)
      rep(constant, length(mz))}. Moreover, note that a valid function 
    has to be vectorized. For the automatic generation
    of those functions from a raw spectrum and further details on the
    meaning of the parameters, see \link{fitModelParameters}. The output
    of a call to \link{fitModelParameters} can directly be plugged into
    \code{getPeaklist} via the argument \code{model.parameters}.} 
\item{binning}{A \code{logical} indicating whether the fitting process
  should be done sequentially in 'bins'. If \code{TRUE}, the spectrum is cut into pieces
  (bins). Each bin is then fitted separately, and the results of all
  bins are combined in the end. Division into bins may be configured using \code{control.binning}. See also
  the 'Details' section below.}
\item{postprocessing}{A \code{logical} indicating whether a
  post-processing correction should be applied to the raw
  peaklist. See also the argument \code{control.postprocessing} and the 'Details' section below.}
\item{trace}{A \code{logical} indicating whether information tracing the different steps of the
  fitting process should be displayed.}
\item{returnbasis}{A \code{logical} indicating whether the matrix of basis functions (template
  functions evaluated at \code{mz}) should be returned. Note that this may be expensive in terms of storage.}
\item{control.basis}{A list of arguments controlling the computation of
  the matrix of basis functions:
  \describe{
    \item{\code{charges}}{The set of charge states present in the
      spectrum.}
    \item{\code{eps}}{Function values below \code{eps} are set equal to
      precisely zero in order to make the basis function matrix sparse.}
  }
}
\item{control.localnoise}{A list of arguments controlling the placement
  and selection of basis functions on the basis of a 'local noise
  level':
  \describe{ 
    \item{\code{quantile}}{A value from \code{0.1, 0.2, ..., 0.9}, specifying
      the quantile of the intensities residing in a sliding
      m/z \code{window} (s. below) to be used as 'local noise level'.}
    \item{\code{factor.place}}{Controls the placement of basis
      functions. A basis function is placed at an element of \code{mz}
      if and only if the intensity at that position exceeds the 'local
      noise level' by a factor at least equal to \code{factor.place}.}
    \item{\code{factor.post}}{Controls which basis functions
      \emph{enter} the postprocessing step. A basis function is
      discarded before the postprocessing step if its estimated
      amplitude does not exceed the \code{factor.post} times the 'local
      noise level'. By default \code{factor.post = 0}. The pre-filtering
    step before postprocessing is mainly done for computational
    speed-up, and \code{factor.post = 0} is supposed to yield the
    qualitatively best solution, though it may take additional time.}
    \item{\code{window}}{The length of the sliding window used to compute
      \code{quantile}, to be specified in Thomson. By default,
      \code{window} is chosen in such a way that it equals the length of
      the support of an 'average' basis function for charge state one.}
    \item{\code{subtract}}{A \code{logical} indicating whether the
      'local noise level' should be subtracted
      from the observed intensities. Setting \code{subtract = TRUE} is
      typically beneficial in the sense that fitting of noise is reduced. 
      }
    }
  }
  \item{control.postprocessing}{A list of arguments controlling the
    postprocessing step (provided \code{postprocessing = TRUE}):
  \describe{ 
    \item{\code{mzfilter}}{Setting \code{mzfilter = TRUE} removes basis
    functions at positions where peak patterns are highly
    improbable to occur, thereby removing peaks from the list which are
    likely to be noise peaks. This filter is sometimes called 'peptide
    mass rule', see \emph{Zubarev et al. (1996): Accuracy Requirements
      for Peptide Characterization by Monoisotopic Molecular
      Measurements, Anal. Chem.,88,4060-4063}.}
    \item{\code{prune}, \code{factor.prune}}{Setting \code{prune = TRUE} activates a crude
      scheme that removes low-intensity peaks (likely to be noise
      peaks), as frequently occurring in regions with extremely intense
      peaks. According to this scheme, a peak is removed from the peak
      list if its amplitude is less than \code{factor.prune} times the
      locally most intense amplitude, where \code{factor.prune}
      typically ranges from \code{0.01} to \code{0.1}.}
    \item{\code{ppm}}{A ppm (= parts per million) tolerance value within
    which basis functions at different m/z positions are considered to
    be merged, s. 'Details' below. By default, that value is computed
    from the spacing of the first two m/z positions.}
    \item{\code{goodnessoffit}}{A \code{logical} indicating whether a local goodness-of-fit
      adjustment of the signa-to-noise ratio should be computed. Yields
      usually more reliable evaluation of the detected patterns, but is
      computationally more demanding.}
   }}
 \item{control.binning}{Controls the division of the spectrum into bins
  (if \code{binning = TRUE}). Based on the 'local noise level'
  described in \code{control.localnoise}, if within
  a range of \code{(1+tol)} Thomson no further significant position occurs, a bin
  is closed, and a new one is not opened  as long as a new significant
  position occurs..
  }
}



\details{\itemize{
    \item While setting \code{binning = TRUE} yields a procedure which
    is less memory consuming than fitting the whole
  spectrum simultaneously (\code{binning = FALSE}), it may be inferior from a quality aspect, since division
  into bins has to be done with care. Otherwise, peak patterns might be
  split up into different bins, which would result into erroneous
  fitting.
\item Postprocessing of the raw list usually yields a
  considerably improved result by counteracting the 'peak-splitting
  phenomenon': due to a limited sampling rate and discrepancies  
  between theoretical and observed peak patterns, several templates at
  adjacent positions are used to fit the same peak pattern.}}
 
\section{Warning}{Although we have tried to choose default values
  expected to produce sensible results, the user should carefully
  examine all options.}
\section{Warning}{Depending on the length and the resolution of the
  raw spectrum, fitting the whole spectrum simultaneously as recommended
  is expensive from a computational point of view, and may take up to several
  minutes per spectrum.}





\value{An object of class \link{peaklist}.}

\seealso{\link{fitModelParameters}, \link{peaklist}}
\keyword{models}
\examples{
### load data
data(toyspectrum)
data(toyspectrumsolution)
mz <- toyspectrum[,"x"]
intensities <- toyspectrum[,"yyy"]
### select mz range
filter <- mz >= 2800 & mz <= 3200
### Extract peak patterns with model = "Gaussian"
sigmafun <- function (mz)  -8.5e-07 * mz + 6.09e-10 * mz^2 + 0.00076
gausslist <- getPeaklist(mz = mz[filter], intensities = intensities[filter],
                   model = "Gaussian",
                   model.parameters = list(sigma = sigmafun,
                                           alpha = function(mz){},
                                           mu = function(mz){}),
                    control.localnoise = list(quantile = 0.5, factor.place = 3))

show(gausslist)
### threshold list at signal-to-noise ratio = 2
peaklist <- threshold(gausslist, threshold = 2)

### Extract peak patterns with model = "EMG" and loss = "L1"
alpha0 <- function(mz) 0.00001875 * 0.5 * 4/3 * mz
sigma0 <- function(mz) 0.00001875 * 0.5 * mz
mu0 <- function(mz) return(rep(-0.06162891, length(mz)))
EMGlist <- getPeaklist(mz = mz[filter], intensities = intensities[filter],
                   model = "EMG", loss = "L1",
                   model.parameters = list(sigma = sigma0,
                                           alpha = alpha0,
                                           mu = mu0),
                   control.localnoise = list(quantile = 0.5, factor.place = 3))
show(EMGlist)
peaklist2 <- threshold(EMGlist, threshold = 2)

### plot results of the 1st list and compare vs. 'truth' 

### 'ground truth'
solution <- toyspectrumsolution[toyspectrumsolution[,1] >= 2800 & toyspectrumsolution[,1] <= 3200,]

visualize(gausslist, mz[filter], intensities[filter], lower = 3150, upper = 3170,
          truth = TRUE,
          signal = TRUE,
          fitted = TRUE,
          postprocessed = TRUE,
          booktrue = as.matrix(toyspectrumsolution),
          cutoff.eps = 0.2)
}
