% genomestructure

\name{genomestruc}

\alias{genomestruc}

\title{
  Genome structure definition
}

\description{
  Defines the experiment genome (or parameter set) structure.
}

\usage{
  genomestruc(n.genes=NULL)
}

\arguments{
  \item{n.genes}{The number of genes (parameters) of each individual
  (experimental set-up). The number of genes (parameters) will be asked
  if no argument is given.}
}

\value{
  A structure which is the first part of an evolutionary algorithm
  structure.
  \item{genes}{An list of genes (parameters); each gene (parameter)
    itself is a list of: 
    \describe{
      \item{name:}{the name of the gene (parameter)}
      \item{min:}{the minimal value that the gene (parameter) can have}
      \item{max:}{the maximum value that the gene (parameter) can have}
      \item{n.alleles:}{the number of the alleles (parameter values) the
	gene (parameter) can have} 
      \item{names:}{a vector containing the names of all the different
      alleles (parameter values); this is optional and set to NULL if the
      number of names given is not equal to the number of alleles
      (parameter values)}  
    }
  }
  \item{generations}{The second part of an evolutionary algorithm
    structure which contains a list of generations (fitnesses, allele
    values etc); the value is initially NULL.}
}

\note{
  \code{genomestruc} is the actual genome (parameter set) definition
  function. However, if package \code{seao.gui} is also loaded, a graphical
  interface can be called by \code{genomestruc.gui}. The return value of
  both functions is exactly the same, since the former is called by the
  latter.
}

\references{
}

\author{
  Kurt Sys (\email{kurt.sys@rug.ac.be})
}

\seealso{
  \link{newgen}
}

\keyword{
  gene, parameter,
  genome, parameter set,
  genome structure, parameter set structure
}

\examples{
}
