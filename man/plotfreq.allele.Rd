% plotfreq.allele

\name{plotfreq.allele}

\alias{plotfreq.allele}

\title{
   Plot a histogram of the allele frequences.
}

\description{
  Makes a histogram of the allele frequencies of one generation of a certain gene.
}

\usage{
  plotfreq.allele(struc.ea,
                  gene=1,
                  gen=NULL,
                  breaks=NULL,
                  ...)
}

\arguments{
  \item{struc.ea}{An object of type evoutionary algorithms. This
    contains the data to be plotted}
  \item{gene}{The number of the gene to be plotted.}
  \item{gen}{A number pointing to which generation to plot. The default
    value is NULL, resulting in the last generation present in the optimization.}
  \item{breaks}{Where the breaks on each plot should be set. See
    function 'hist' for more information.}
  \item{...}{Arguments to set for making the histogram. See \link{hist}.}
}

\details{
}


\references{
}

\author{
   Kurt Sys (\email{kurt.sys@rug.ac.be})
}

\seealso{
}

\keyword{
  plot, evolution, gene evolution
}

\examples{
}
