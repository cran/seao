% next generation

\name{nextgen}

\alias{nextgen}

\title{
  Next generation
}

\description{
  Calculates a next generation, i.e. chooses and selects parents and
  performs crossover and mutation.
}

\usage{
  nextgen(struc.ea,
          n.ind=NULL,
          gen.parent=NULL,
          selection=list(base="fit", rescale=0),
          corate=90,
          mutation=list(base="unif", spread=1, rate=15))
}

\arguments{
  \item{struc.ea}{An evolutionary algorithm structure. It must contain
    at least \code{$genes}, which is a list describing the structure of
    the genome (parameter set).}
  \item{n.ind}{The number of individuals one wants to generate. Default
    value is set to 25.}
  \item{gen.parent}{A vector containing the numbers of the parent generations
    (batches). As default, the function includes all previous
    calculated generations (batches) for choosing the parents of
    the next generation (batch).}
  \item{selection}{A list containing the base for selection, and possible
    selection. The base can be either "fit" (standard value) or
    "rank"; rescaling a number, standard value is 0, but 1.5 is a
    value which gives often nice results.}
  \item{corate}{The parameter value recombination rate, i.e. the percentage
    individuals which should be recombined.}
  \item{mutation}{Some evolutionary parameters to set. This contains whether a
    uniform ('unif') or normal ('norm') distribution should be
    used, whether the spread value (if this is set to 0, a value
    of 1 is used), and to set the paramater value randomization
    rate.}  
}

\value{
  The function returns an evolutionary algorithm object:
  \item{genes}{a list containing the genome (parameter set)
    stucture. This is the same as the genome (parameter set) structure
    of in the input (\code{struc.ea$genes}).}
  \item{generations}{a list of structures with information of each
    generation (batch): 
    \describe{
      \item{fit:}{a vector containing the fitnesses of all individuals.
        These values are all set to NA, since the fitnesses of
        the new individuals are not known yet.}
      \item{allele:}{a matrix containing the allele-values of all
        individuals (experimental set-ups) of the new generation
        (batch)}
      \item{parents:}{a string vector with the numbers of the parent
        generations.}
      \item{selection:}{a list containing the selection base and
	rescaling.}
      \item{crossover:}{the crossover rate.}
      \item{mutation:}{a list containing the mutation parameters.}
  }}
}

\note{
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
