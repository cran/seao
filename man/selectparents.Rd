% selectparents

\name{selectparents}

\alias{selectparents}

\title{
  Select the parents to mate
}

\description{
  Select the parents which are the basis for the next (child)
  generation. The individuals (set-ups) with highest fitness have the
  highest chance to be chosen, so some parents may be chosen a few
  times. 
}

\usage{
  selectparents(struc.ea, gen=NULL,
                method=list(base="fit", rescale=0))
}

\arguments{
  \item{struc.ea}{An evolutionary algorithm structure}
  \item{gen}{The number of the generation (batch) for which the parents
    have to be chosen. The standard value is the last generation.}
  \item{method}{A list containing the base for selection, and possible
    selection. The base can be either "fit" (standard value) or "rank";
    rescaling a number, standard value is 0, but 1.5 is a value which
    gives often nice results.}
}

\details{
  The method is a list of \code{base} and \code{rescale}. \code{base}
  gives which basis should be used for calculating the weights of the
  different individuals:
  \describe{
    \item{"fitness"}{The fitness values are the weights given, before
      rescaling}
    \item{"rank"}{The individuals are ranked from 1 to the number
      of individuals in the generation. The higher the number, the higher
      the fitness of the individual. If two individuals have exactly the
      same fitness, the average of the next and previous ranking is used,
      e.g. imagene that two fitnesses should be at ranking 7 and
      8. However, they have the same value, so they are both set to 7.5
      (in this case, the sum is still 15).}
  }
  In many cases, the difference between the lowest and highest
  weights are very high. This would result in choosing the same
  parent(s) every time. In order to decrease this extreme selection,
  one can rescale the weights. Rescaling is done so that the highest
  fitness is \code{recale} times the average value of the fitnesses.}

  For the selection itself, a roulette-like system is used. The chance
  that a certain numer is chosen is linear correlated with the weight
  (which is either the fitness or ranking value).
}

\value{
  The function returns an evolutionary algorithm object, which is
  actually the same as the input argument \code{struc.ea}, but with new
  individuals (experimental set-ups) and the selection method added in
  the generation (batch) of selection: 
  \item{genes}{a list containing the genome (parameter set)
    stucture. This is the same as the genome (parameter set) structure
    of in the input (\code{struc.ea$genes}).}
  \item{generations}{a list of structures with information of each
    generation: 
    \describe{
      \item{fit:}{a vector containing the fitnesses of all individuals
	(experimental set-ups) of the generation (batch); this are just
	copies of the fitnesses of the selected parents.}
      \item{allele:}{a matrix containing the allele-values of all
	individuals (experimental set-ups) of the generation (batch);
	this are just copies of the alleles (parameter values) of the
	selected parents.} 
      \item{parents:}{a string vector with the numbers of the parent
	generations.}
      \item{selection:}{a list containing the selection base and
	rescaling.}
	
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
  \link{genomestruc}, \link{newgen}, \link{chooseparents}
}

\keyword{
  parent
  selection
}

\examples{
}
