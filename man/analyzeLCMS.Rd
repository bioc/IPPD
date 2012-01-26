\name{analyzeLCMS}
\alias{analyzeLCMS}
\title{Functionality for LC-MS data}
\description{The function \code{analyzeLCMS} applies the function
\code{getPeaklist} to all scans of the dataset. The sweep line scheme of
Schulz-Trieglaff et al. (2008) is then used to find retention time
intervals of consecutive signals.}

\usage{analyzeLCMS(data, arglist.fitModelParameters = list(),
                         arglist.getPeaklist = list(),
                         arglist.threshold = list(),
                         arglist.sweepline = list(),
                         trace = TRUE)}

\arguments{
  \item{data}{\code{data} can be one of the following:
    \itemize{\item{A \code{character} for the path to the file that is
	to be processed. The file has to be either an \code{mzXML} file
      or a three-column matrix whose columns contains retention
  times, m/z positions and intensities (exactly in that order).}
\item{An object of class \code{mzXML}.}
\item{A \code{matrix} or \code{data.frame} containing three columns
  containing retention
  times, m/z positions and intensities (exactly in that order).}}
}
  \item{arglist.fitModelParameters}{Optional arguments to be passed to
  \link{fitModelParameters}.}
  \item{arglist.getPeaklist}{Optional arguments to be passed to
  \link{getPeaklist}.} 
\item{arglist.threshold}{Optional arguments to be passed to the method
\code{threshold}, see \link{peaklist} for details.}
  \item{arglist.sweepline}{Optional arguments to be passed to the
    function \link{sweepline}.}
  \item{trace}{A \code{logical} controlling whether status information
  is displayed.}
}

\value{A \code{list} composed of
\describe{
     \item{peaklists}{A list of matrices as in the slot
       \code{peaklistprocessed} of an object of class \code{peaklist}. These
       matrices are obtained from running \code{getPeaklist} for
       all scans contained in \code{data} and applying hard thresholding
     as implemented in the method \code{threshold}.}
     \item{boxes}{A matrix of retention time intervals as
       returned as output by a call to the function \link{sweepline}.}
  }
}

\examples{}

\references{0. Schulz-Trieglaff and R. Hussong and C. Groepl and
  A. Leinenbach and A. Hildebrandt and C. Huber and K. Reinert. (2008)
  \emph{Journal of Computational Biology, 15, 685-704}}

\keyword{models}
