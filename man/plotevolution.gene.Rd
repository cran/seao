% plotevolution.gene

\name{plotevolution.gene}

\alias{plotevolution.gene}

\title{
   Plot the evolution of the alleles of one gene for an optimization.
}

\description{
  Makes a histogram of the evolution of the alleles of one gene. The
  generations to plot can be chosen.
}

\usage{
  plotevolution.gene(struc.ea,
                     gene=1,
                     gens=NULL,
                     breaks=NULL,
                     main="Gen",
                     cols=NULL,
                     rows=NULL,
                     ...)
}

\arguments{
  \item{struc.ea}{An object of type evoutionary algorithms. This
    contains the data to be plotted}
  \item{gene}{The number of the gene to be plotted.}
  \item{gens}{A vector containing the generations to plot. The default
    value is NULL, giving all generations present in the optimization.}
  \item{breaks}{Where the breaks on each plot should be set. See
    function 'hist' for more information.}
  \item{main}{The main title of the plot.}
  \item{cols}{The number of columns, i.e. the number of histograms on
    one row. If set to NULL, argument 'rows' is checked.}
  \item{rows}{The number of rows, i.e. the number of histogram in one
    column. If 'cols' is set, this argument is omitted. Is both 'cols'
    and 'rows' are set to NULL, there's an automatic calculation of the
    number of rows and columns.}
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
  plot, evolution, gene evolution
}

\examples{
}
