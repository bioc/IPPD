\name{sweepline}
\alias{sweepline}
\title{Sweep line scheme to process LC-MS data}
\description{The function implements a version of the sweep line
  scheme suggested in Schulz-Trieglaff et al. (2008) to aggregate
  the results of multiple scans of an LC-MS dataset.
}

\usage{sweepline(peaklists, rt, tol = 100, gap = 2, minboxlength = 5)}

\arguments{ 
  \item{peaklists}{A \code{list} of matrices as in the slot
  \code{peaklistprocessed} of an object of class \code{peaklist}. It is
  assumed that the entries of \code{peaklists} correspond to the
  sequence retention times when ordered in ascending order.}
  \item{rt}{A \code{numeric} vector of retention times.}
  \item{tol}{A \code{m/z} tolerance (in parts per million, ppm) below
  which signals of different scans (retention times) are considered for merging.} 
\item{gap}{Retention time intervals are formed by merging signals of
adjacent (with respect to retention time) scans. An intervals is closed
(i.e. no further signals are added) once no additional signal has been added for more than \code{gap} consecutive scans.}
  \item{minboxlength}{Minimum number of adjacent signals found such that
  the corresponding retention time interval is kept in the output.} 
}

\value{A \code{matrix} having the following columns. Each row
       corresponds to a retention time interval of consecutive signals.
       \describe{  
	 \item{\code{loc}}{\code{m/z}-position corresponding to the interval.}
	 \item{\code{charge}}{charge of the signals found.}
	 \item{\code{quant}}{Cumulative intensities of the signals found.}
	 \item{\code{rt_begin},\code{rt_begin}}{Boundaries of the retention time interval.}
	 \item{\code{npeaks}}{Total number of scans containing signals
       assigned to the interval.}
	 \item{\code{gapcount}}{Total number of scans corresponding to
       the retention time interval not containing a signal.}
	 }
}

\examples{}

\references{0. Schulz-Trieglaff and R. Hussong and C. Groepl and
  A. Leinenbach and A. Hildebrandt and C. Huber and K. Reinert. (2008)
  \emph{Journal of Computational Biology, 15, 685-704}}

\seealso{\link{analyzeLCMS}}
\keyword{models}
