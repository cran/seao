% mutation

\name{mutation}

\alias{mutation}

\title{
  Perform mutation
}

\description{
  Ad random parameter values are selected and changed, giving a uniform
  or normal probability distribution.
}

\usage{
  mutation(struc.ea, gen=NULL, method=list(base="unif", spread=1, rate=15))
}

\arguments{
  \item{struc.ea}{An evolutionary algorithm structure}
  \item{gen}{The number of the generation (batch) for which the parents
    have to be chosen. The standard value is the last generation.}
  \item{method}{Some evolutionary parameters to set. This contains
    whether a uniform ('unif') or normal ('norm') distribution should be
    used, whether the spread value (if this is set to 0, a value of 1 is
    used), and to set the paramater value randomization rate.}
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
	a portion of these values are randomized paramater values.} 
      \item{parents:}{a string vector with the numbers of the parent
	generations.}
      \item{selection:}{a list containing the selection base and
	rescaling.}
      \item{crossover:}{the crossover rate.}
      \item{mutation:}{a list containing the mutation parameters. This
	value is added.}
	
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
  \link{genomestruc}, \link{newgen}, \link{chooseparents},
  \link{selectparents}, \link{crossover}
}

\keyword{
  child
  mutation
}

\examples{
}
