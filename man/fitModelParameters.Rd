\name{fitModelParameters}
\alias{fitModelParameters,numeric,numeric-method}
\alias{fitModelParameters}
\alias{fitModelParameters-methods}
\title{Peak parameter estimation}
\description{In the template-based approach of this package, each
  template/peak pattern is composed of several single basic peaks.
  The shape of such a basic peak may be modeled either as Gaussian
  or as Exponentially Modified Gaussian (EMG). The second model assumes
  that the shape of each peak equals the shape of the function one
  obtains when convolving the probability density function of a Gaussian
  distribution with the probability density function of the exponential
distribution. This is a more complex as well as more
flexible model, since it allows one to account for skewness. Both peak models depend on a set of
parameters which are usually unknown a priori. Moreover, these parameters
tend to vary over the spectrum. The documented method provides the
following functionality: given a raw spectrum and a linear model that describes
how the set of parameters varies as a function of m/z, we detect
well-resolved peaks within the spectrum. For each detected peak, we
determine the parameters of our model in such a way that the resulting
peak shape matches that of the detected peak as good as possible in a
least-squares sense. In this manner, we obtain a set of estimated
parameters for different m/z positions. These estimates are used as response,
the m/z positions as explanatory variable in a linear regression model (not
necessarily linear in m/z). To be resistant to outliers (e.g. as
occurring due to overlapping peaks), we use least absolute deviation
regression to infer the model parameters of these linear models. The
result can directly be used for the argument \code{model.parameters} in
the main method \code{getPeakList}.}

\usage{fitModelParameters(mz, intensities, model = c("Gaussian", "EMG"),
                                           fitting = c("most_intense", "model"), 
                                           formula.alpha =  formula(~1),
                                           formula.sigma = formula(~1),
                                           formula.mu = formula(~1),
                                           control = list(window = 6,
                                           threshold = NULL))
}
\arguments{
   \item{mz}{A \code{numeric} vector of m/z (mass/charge) values (in Thomson).}
   \item{intensities}{A \code{numeric} vector of intensities
    corresponding to \code{mz}.}
  \item{model}{Basic model for the shape of a single peak. Must be
    \code{"Gaussian"} or \code{"EMG"} (exponentially modified
    Gaussian). See details below.}
  \item{fitting}{A character specifying the mode of peak extraction and
    -fitting. If \code{fitting = "most_intense"}, then the most intense
    peak is detected and the parameter(s) is (are) fitted using only this
    peak. Then, the resulting functions of the form
    \code{parameter(mz)} are all constant functions. If \code{fitting =
      "model"}, as many peaks as possible satisfying the criteria of
    \code{control} are used to estimate linear models of the form
    \code{parameter(mz) = beta_0 + beta_1 g_1(mz) + ... + beta_p g_p(mz)},
    where the \code{g}'s are fixed functions. Their specification is
    performed by specifying one-sided formulae according to the
    Wilkinson-Roger notation for linear models as in the function
    \code{lm}. The model formulae have to satisfy the following criteria:
    \itemize{\item The formula is one sided, i.e. \emph{no} term
	appears on the left hand side of \code{~}.
      \item The right hand side consists only of functions in
	\code{mz}, and \code{mz} is the only variable that may be used.
	Product terms involving \code{*} are \emph{not} admitted.
	\item \bold{Important}: Note that, for example \code{~ 1 +
	  mz + sqrt(mz)} is a valid formula in the sense that no error
	will occur, but it does \emph{not} correspond to the linear model 
	\code{parameter(mz) = beta_0 +  beta_1 mz + beta_2
	  sqrt(mz)}. The correct model formula instead
	reads \code{~ 1 + mz + I(sqrt(mz))}, i.e. each function has to
	be bracketed by \code{I()}.}
      }
  \item{formula.alpha}{A one-sided formula describing the dependence
                       of the EMG-parameter \code{alpha} as a function
		       of m/z, s. \code{fitting}. The default assumes
		       that the parameter is independent of m/z, in
		       which case one obtains a function returning the
		       median of the estimates obtained for different
		       detected peaks. \code{formula.alpha},
		       \code{formula.sigma} and \code{formula.mu} are
		       needed if and only if \code{fitting = "model"}.}
  \item{formula.sigma}{Parameter used for both peak models. For further
    information, see \code{formula.alpha}}
  \item{formula.mu}{See \code{formula.alpha}. Used only if \code{model =
      "EMG"}.}
  \item{control}{A \code{list} controlling peak detection. The parameter
  \code{window} refers to the minimal resolution of a peak. According to
  \code{window}, a sequence of intensities at adjacent m/z positions is
  considered as peak if and only if there are at least \code{window}
  m/z positions with \emph{increasing} intensity followed by a second sequence
  of \code{window} m/z
  positions with \emph{decreasing} intensity. If in addition \code{threshold}
  is specified, only peaks whose maximum intensity is equal to or
  greater than \code{threshold} is considered. \bold{Note}: Usually, \code{threshold}
is specified, since otherwise the maximum intensity among the complete
spectrum minus some epsilon is taken as \code{threshold}.}
  }



\details{Let the variable \code{x} represent m/z. Then \code{model =
    "Gaussian"} assumes that a single peak can be described as \cr
  \cr
  \code{gaussfun(x;sigma,mu) = exp(-(x - mu)^2/sigma)} 
  \cr
  \cr
    The parameter \code{mu} is \emph{not} considered as model parameter:
    in the computation of the resulting basis function matrix, \code{mu}
    is always set to a known m/z position where the leading peak of a peak
    pattern might be present.\cr
    \code{Model = "EMG"} assumes that a single peak can be described as
    \cr
    \cr
  \code{EMG(x;alpha,sigma,mu) = exp(sigma^2/(2 * alpha^2) + (mu -
    x)/alpha) (1 - Phi(sigma/alpha + (mu - x)/(sigma)))/alpha}, \cr
  \cr
   where \code{Phi} represents the cumulative density function of the
    standard Gaussian distribution. Alternatively,
    \code{EMG(.;alpha,sigma,mu)} can be expressed as \cr
    \code{EMG(x;alpha,sigma,mu) = (phi ** gamma)(x)},
    \cr
    where \code{**} denotes convolution, \code{phi} is the density function of
    the Gaussian distribution with mean \code{mu} and standard deviation \code{sigma}
    and \code{gamma} is the density function of an exponential distribution
    with expectation \code{alpha}.\cr
    The parameters of \code{EMG} can be interpreted as follows.
    \describe{
      \item{\code{alpha}}{The lower \code{alpha}, the more the shape of the peak
    resembles that of a Gaussian. Conversely, large values of
    \code{alpha} lead to long right tails.}   
      \item{\code{sigma}}{Controls the width of the peak (together with
	\code{alpha}).}
      \item{\code{mu}}{A location parameter. Note that in general \code{mu}
    does \emph{not} coincide with the mode of \code{EMG}. Therefore, if
    \code{model = "EMG"}, all three parameters are estimated from
    detected peaks.}
    }  
    Moreover, the skewness of \code{EMG} is characterized by the ratio \code{alpha/sigma}.
  }
  
  \section{Warning}{
     Parameter estimation by fitting detected peaks is possible only if
     single peaks are sufficiently well-resolved. A peak composed of, say,
     five (m/z, intensity) pairs, is inappropriate to infer three parameters.
   }

   \section{Warning}{
     The choice \code{model = "EMG"} in conjunction with \code{fitting =
     "model"} can be extremely slow (taking up to
     several minutes of computing time) if many peaks are detected and
     fitted. This is caused by a grid search over a grid of \code{10^6}
     different combinations of \code{alpha}, \code{sigma} and \code{mu}
     performed prior to nonlinear least squares estimation in order to
     find suitable starting values.
    }
  

\value{An object of class \link{modelfit}.}

\seealso{\link{getPeaklist}, \link{modelfit}}
\keyword{models}
\examples{
### load data
data(toyspectrum)
### estimate parameter sigma of a Gaussian model,
### assumed to be independent of m/z

simplegauss <- fitModelParameters(toyspectrum[,1],
             toyspectrum[,2],
             model = "Gaussian",
             fitting = c("model"),
             formula.sigma = formula(~1),
             control = list(window = 6, threshold = 1))

show(simplegauss)
visualize(simplegauss, type = "peak", xlab = "m/z", ylab = "intensity",
          main = "Gaussian fit")

### fit the model sigma(m/z) = beta_0 + beta_1 m/z + beta_2 m/z^2

gaussquadratic <- fitModelParameters(toyspectrum[,1],
             toyspectrum[,2],
             model = "Gaussian",
             fitting = "model",
             formula.sigma = formula(~mz + I(mz^2) ),
             control = list(window = 6, threshold = 1))

show(gaussquadratic)
visualize(gaussquadratic, type = "model", modelfit = TRUE)

### estimate parameters for EMG-shaped peaks

EMGlinear <- fitModelParameters(toyspectrum[,1],
             toyspectrum[,2],
             model = "EMG",
             fitting = "model",
             formula.alpha = formula(~mz),
             formula.sigma = formula(~mz),
             formula.mu = formula(~1),
             control = list(window = 6, threshold = 1))

show(EMGlinear)

visualize(EMGlinear, type = "peak", xlab = "m/z", ylab = "intensities",
          main = "EMG fit")

visualize(EMGlinear, type = "model", parameters = c("alpha", "sigma"), modelfit = TRUE)
} 
 
