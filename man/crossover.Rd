% crossover

\name{crossover}

\alias{crossover}

\title{
  Perform crossover
}

\description{
  The parameter values are recombined to make new batches. The rate can
  be set, as well as the generation on which this recombination should happen.}
}

\usage{
  crossover(struc.ea, gen=NULL, rate=90)
}

\arguments{
  \item{struc.ea}{An evolutionary algorithm structure}
  \item{gen}{The number of the generation (batch) for which the parents
    have to be chosen. The standard value is the last generation.}
  \item{rate}{The parameter value recombination rate, i.e. the
    percentage individuals which should be recombined.}
}

\details{
}

\value{
  The function returns an evolutionary algorithm object, which is
  actually the same as the input argument \code{struc.ea}, but in which
  the a generation is transformed:
  \item{genes}{a list containing the genome (parameter set)
    stucture. This is the same as the genome (parameter set) structure
    of in the input (\code{struc.ea$genes}).}
  \item{generations}{a list of structures with information of each
    generation: 
    \describe{
      \item{fit:}{a vector containing the fitnesses of all
	individuals. These values are all set to NA, since the fitnesses
	of the new individuals are not known yet.}
      \item{allele:}{a matrix containing the allele-values of all
	individuals (experimental set-ups) of the generation (batch);
	these values are recombined paramter values.} 
      \item{parents:}{a string vector with the numbers of the parent
	generations.}
      \item{selection:}{a list containing the selection base and
	rescaling.}
      \item{crossover:}{the crossover rate; this value is added.}
	
  }}
}

\note{
  As for all functions in this package, evolutionary terminology is used
  as described in \code{\link{seao.terminology}}. 
}

\author{
  Kurt Sys (\email{kurt.sys@rug.ac.be})
}

\seealso{
  \link{genomestruc}, \link{newgen}, \link{chooseparents}, \link{selectparents}
}

\keyword{
  child
  crossover
}

\examples{
}
