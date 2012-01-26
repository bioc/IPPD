\name{IPPD-package}
\alias{IPPD-package}
\docType{package}
\title{Peak pattern deconvolution for Protein Mass Spectrometry by
  non-negative ls/lad template matching}
\description{The package provides functionality to extract isotopic peak
  patterns from raw mass spectra. This is done by fitting a large set of
  template basis functions to the raw spectrum using nonnegative least
  squares (ls) or nonnegative least absolute deviation (lad). Ideally, the nonnegativity constraint in combination with
  nonnegativity of the template basis functions effects
  that templates not matching the data are assigned an extremely low
  weight such that one can easily identify isotopic patterns present and
  not present in the spectrum. In practice, the templates only
  approximate the peak patterns, where the quality the approximation
  crucially depends on how well the shapes of the templates fit the
  isotopic patterns contained in a spectrum. For this reason, the
  package offers the flexible function \link{fitModelParameters} which
  tries to estimate model parameters, e.g. the width of a
  Gaussian bump, in a way tailored to the peak shapes in the data. As
  second peak model in addition to the standard Gaussian, the package
  offers full support for the Exponential Modified Gaussian.\cr 
  The function \link{getPeaklist} predicts the set of isotopic peak
  patterns present in the spectrum in a fully automatic, yet
  customizable way. The main benefits of our approach are that
  \enumerate{
    \item Overlapping peak patterns can be resolved.
    \item The complete spectrum can be processed as a whole or in
      large sections by exploiting the sparse nature of the problem.
     \item The set of parameters in \link{getPeaklist} are easy to
       interpret and require only very basic knowledge of statistics.
     \item A theoretically well-founded post-processing procedure is used.
    \item The result can be analyzed visually in a detailed way using the function \code{visualize}.
    }
\tabular{ll}{
Package: \tab IPPD\cr
Type: \tab Package\cr
Version: \tab 0.99.2\cr
Date: \tab 2012-01-17\cr
License: \tab GPL (version 2 or later) \cr
}
}
\author{
  Martin Slawski \email{ms@cs.uni-saarland.de}, \cr
  Rene Hussong \email{rene.hussong@uni.lu}, \cr
  Matthias Hein \email{hein@cs.uni-saarland.de}

  Maintainer: Martin Slawski \email{ms@cs.uni-saarland.de}.
}
\keyword{models}
 
