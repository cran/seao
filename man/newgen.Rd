% new generation

\name{newgen}

\alias{newgen}

\title{
  New generation
}

\description{
  Generates a new (random) generation of individuals (experiments) using an
  already defined genomestructure.
}

\usage{
  newgen(struc.ea, n.ind=25, method="random")
}

\arguments{
  \item{struc.ea}{An evolutionary algorithm structure. It must contain
    at least \code{$genes}, which is a list describing the structure of
    the genome (parameter set).}
  \item{n.ind}{The number of individuals one wants to generate. Default
    value is set to 25.}
  \item{method}{This describes the method that must be used the make a
    (random) generation (batch). There are two possibilities,
    i.e. \code{random}, for pure random allele (parameter value)
    generation and \code{maxgenediv}, in which case the diversity per
    gene is maximised. The default value is "random".} 
}

\value{
  The function returns an evolutionary algorithm object:
  \item{genes}{a list containing the genome (parameter set)
    stucture. This is the same as the genome (parameter set) structure
    of in the input (\code{struc.ea$genes}).}
  \item{generations}{a list of structures with information of each
    generation (batch): 
    \describe{
      \item{fit:}{a vector containing the fitnesses of all individuals
	(experimental set-ups).}
      \item{allele:}{a matrix containing the alleles (parameter values of
	all individuals (experimental set-ups).} 
      \item{parents:}{a string vector with the numbers of the parent
	generations (batches). For this function, this is set to a vector with
	values "-1" and the method used ("random" or "maxgenediv")}
  }}
}

\note{
  \code{newgen} is the actual function that calculates the (random)
  alleles of a  new generation. However, if package \code{seao.gui} is
  also loaded, a graphical interface can be called by
  \code{newgen.gui}. The return value of both functions is exactly the
  same, since the former is called by the latter. 

  As for all functions in this package, evolutionary terminology is used
  as described in \code{\link{seao.terminology}}..
}

\references{
}

\author{
  Kurt Sys (\email{kurt.sys@rug.ac.be})
}

\seealso{
  \link{genomestruc}
}

\keyword{
  initial generation
  new generation
}

\examples{
}
