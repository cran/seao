% choose parents

\name{chooseparents}

\alias{chooseparents}

\title{
  Choose parent generations
}

\description{
  Choose the best individuals (set-ups) of the given generations
  (batches) as possible parents. Parent selection is generally the next
  step. 
}

\usage{
  chooseparents(struc.ea, n.ind=NULL, gen.parent=NULL)
}

\arguments{
  \item{struc.ea}{An evolutionary algorithm structure}
  \item{n.ind}{Number of individuals (set-ups) present in the new
    generation (batch); this is equal to the number of parents chosen in
    this function. If this argument is not given, the average of the
    number of individuals (set-ups) in the chosen parent generations
    (batches) is taken as default.}
  \item{gen.parent}{A vector containing the numbers of the parent
    generations (batches). As default, the function includes all
    previous calculated generations (batches) for choosing the parents
    of the next generation (batch).} 
}

\value{
  The function returns an evolutionary algorithm object, which is
  actually the same as the input argument \code{struc.ea}, but with one
  generation (batch) added, i.e. the chosen parents:
  \item{genes}{a list containing the genome (parameter set)
    stucture. This is the same as the genome (parameter set) structure
    of in the input (\code{struc.ea$genes}).}
  \item{generations}{a list of structures with information of each
    generation: 
    \describe{
      \item{fit:}{a vector containing the fitnesses of all individuals
	(experimental set-ups) of the generation (batch); this are just
	copies of the fitnesses of the chosen parents.}
      \item{allele:}{a matrix containing the allele-values of all
	individuals (experimental set-ups) of the generation (batch);
	this are just copies of the alleles (parameter values) of the
	chosen parents.} 
      \item{parents:}{a string vector with the numbers of the parent
	generations.}
  }}
}

\note{
  As for all functions in this package, evolutionary terminology is used
  as described in \code{\link{seao.terminology}}. 
}

\references{
}

\author{
  Kurt Sys (\email{kurt.sys@rug.ac.be})
}

\seealso{
  \link{genomestruc}, \link{newgen}, \link{selectparents}

  \link{seao.terminology}
}

\keyword{
  parent
}

\examples{
}
