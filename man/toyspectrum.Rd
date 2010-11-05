\name{toyspectrum}
\alias{toyspectrum}
\alias{toyspectrumsolution}
\docType{data}
\title{Simulated protein mass spectrum}
\description{
A simulated dataset for testing purposes. m/z ranges from 800 to about 4000 
Thomson, sampled at approximately 118,000 points. The isotopic patterns present 
in this artificial dataset exactly obey the Averagine model in
combination with Exponentially Modified Gaussians peak shapes, s. 
\link{fitModelParameters} if no noise were present. The noise has been generated
from two sources: the first is a baseline pattern of Gaussian bumps with low 
intensity and the second one is additive noise drawn from a truncated Gaussian 
distribution. The dataset contains about 300, partially heavily overlapping, peak
patterns whose data are stored in \code{data(toyspectrumsolution)}.
}
\usage{data(toyspectrum)}

\examples{
data(toyspectrum)
plot(toyspectrum)
}

\keyword{models}
