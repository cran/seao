% plotevolution.fit

\name{plotevolution.fit}

\alias{plotevolution.fit}

\title{
   Plot the evolution of the fitness of an optimization.
}

\description{
  Makes a plot of the evolution of the fitness. One can choose whether
  the minimum, maximum, average or all fitness is plotted, or if
  boxplots must be used, of if the standard deviation must be
  plotted. Also, the generations to plot can be chosen.
}

\usage{
  plotevolution.fit(struc.ea, gens=NULL, show="average", plot=TRUE, ...)
}

\arguments{
  \item{struc.ea}{An object of type evoutionary algorithms. This
    contains the data to be plotted}
  \item{gens}{A vector containing the generations to plot. The default
    value is NULL, giving all generations present in the optimization.}
  \item{show}{Defines what should be plotted: 'min' for minimum
    fitnesses, 'max' for maximum fitnesses, 'average' for average
    fitnesses, 'all' for all fitnesses, 'box' for boxplots of fitnesses
    per generations and 'stdev' for the standard deviations.}
  \item{plot}{If set to FALSE, the result is not plotted but assigned
    to an object or written to the screen.}
  \item{...}{Arguments to set for plotting. See \link{plot}}
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
  plot, evolution, fitness evolution
}

\examples{
}
