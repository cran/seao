% terminology

\name{seao.terminology}

\alias{seao.terminology}

\title{
   Short overview of terminology
}

\description{
   This package uses a terminology which may seem somewhat strange and
   not obvious. The function does actually notting special, it just shows
   a short overview of the terms used in this package)
}

\usage{
  seao.terminology()
}

\arguments{
}

\details{
   This package uses a terminology which may seem somewhat strange and
   not obvious. The reason why this terminology has been used is however
   quite simple: it are genetic (or evolutionary) algoritms so the
   evolutionary terms are used:
   \describe{
      \item{population}{A number of experimental set-ups. This can be one
        generation or a list of several generations.}
      \item{generation}{All experimental set-ups of one experiment. One
        experiment is defined as a number of parallel experimental
        set-ups.}
      \item{individual}{One experimental set-up with own parameter-values.}
      \item{gene}{Parameter that can be changed for each experimental
        setup. It are these parameters that are optimized.}
      \item{allele}{Value of a parameter in an experimental setup.}
      \item{fitness}{The output of each experimental setup, i.e. the
        property you are using for optimizing the parameters.}
   }
}


\references{
}

\author{
   Kurt Sys (\email{kurt.sys@rug.ac.be})
}

\seealso{
}

\keyword{
   terminology
}

\examples{
}
