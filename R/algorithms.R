#    sgao: Simple Genetic Algorithms for Optimization
#    R library that contains functions for solving non-linear optimization
#    using (simple) genetic algorithms and some appropriate plots
#    Copyright (C) 2003 Kurt Sys -- LabMET, University of Ghent
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


                                        # Small internal function, nothing special
.add.gen<-function(list.gen, gen.new) {
    i.new<-length(list.gen) + 1
    list.gen[[i.new]]<-gen.new

    return(list.gen)
}





genomestruc<-function(n.genes=NULL) {
    lstruc.genes <- NULL
    if (!length(n.genes)) {      
        cat ("Give the number of genes\n")
        n.genes<-as.integer(readline())
    }


    for (li.gene in 1:n.genes) {
        cat(paste ("\n---- GENE ", li.gene, "\n"))
        cat("Name of the gene\n")
        str.name<-readline()
        cat("Give the minimum value\n")
        d.min<-as.double(readline())
        cat("Give the maximum value\n")
        d.max<-as.double(readline())
        cat("Give the step (precision)\n")
        d.step<-as.double(readline())
        d.max<-d.min + round((d.max-d.min)/d.step)*d.step
        cat(paste("Maximum set to:", d.max, "\n"))
        i.alleles<-(d.max-d.min)/d.step + 1
        cat(paste("There are", i.alleles, "alleles\n"))
        cat("Optional: Give a list of names for the different alleles (seperated by a ';')\n")
        vstr.names<-readline()                           
        vstr.names<-strsplit(vstr.names, ";")[[1]]
                                        # todo: remove spaces at beginning (perl?)
        if (length(vstr.names) != (d.max-d.min)/d.step + 1) {
            vstr.names<-NULL
            cat("WARNING: Number of allelenames not equal to number of alleles; namevector set to NULL\n\n")
        } else {
            cat("\n\n")
        }
        lstruc.genes[[li.gene]]<-list(name=str.name,
                                      min=d.min,
                                      max=d.max,
                                      n.alleles=i.alleles,
                                      names=vstr.names)
    }
   
    return(list(genes=lstruc.genes))
}





newgen<-function(struc.ea,
                 n.ind=25,
                 method="random") {



    makemaxgenedivgen<-function() {
        vi.allgeneval <- NULL
        for (li.gene in 1:length(struc.ea$genes)) {
            d.f <- n.ind/struc.ea$genes[[li.gene]]$n.alleles
            vi.geneval<-NULL
                                        # calculate gene values with maximum spread per gene
            for (li.ind in 1:n.ind) {
                vi.geneval<-c(vi.geneval, floor((li.ind-0.5)/d.f)*
                              (struc.ea$genes[[li.gene]]$max -
                               struc.ea$genes[[li.gene]]$min)/(struc.ea$genes[[li.gene]]$n.alleles-1) +
                              struc.ea$genes[[li.gene]]$min)
            }
                                        # scramble gene values
            for (li.ind in 1:n.ind) {
                d.random <- floor(runif(1)*n.ind) + 1
                buffer <- vi.geneval[d.random]
                vi.geneval[d.random] <- vi.geneval[n.ind]
                vi.geneval[n.ind] <- buffer
            }
            vi.allgeneval <- cbind(vi.allgeneval, vi.geneval)
        }

        return(matrix(vi.allgeneval, nrow=n.ind))
    }



    makerandomgen<-function() {
        vi.allgeneval <- NULL
        for (li.gene in 1:i.genes) {
            for (li.ind in 1:n.ind) {
                vi.allgeneval <- c(vi.allgeneval, floor(struc.ea$genes[[li.gene]]$n.alleles*runif(1)) *
                                   (struc.ea$genes[[li.gene]]$max -
                                    struc.ea$genes[[li.gene]]$min)/(struc.ea$genes[[li.gene]]$n.alleles-1) +
                                   struc.ea$genes[[li.gene]]$min)
            }
        }

        return(matrix(vi.allgeneval, nrow=n.ind))
    }



                                        # start function newgen
    i.genes <- length(struc.ea$genes)
    list.gen <- list()

    list.gen$fit <- c(rep(NA, n.ind))
    switch(method,
           random = list.gen$allele <- makerandomgen(),
           maxgenediv = list.gen$allele <- makemaxgenedivgen()
           )

    return(list(genes=struc.ea$genes,
                generations=.add.gen(struc.ea$generation,
                                     list(fit=list.gen$fit,
                                          allele=list.gen$allele,
                                          parents=c(-1, method) ))
                )
           )
}





chooseparents<-function(struc.ea,
                        n.ind=NULL,
                        gen.parent=NULL) {
                                        # use all previous generations if none is given
    if (!length(gen.parent))
        gen.parent <- c(1:length(struc.ea$generations))
    i.parentgens <- length(gen.parent)
    
    list.parentpop <- NULL
    
                                        # Combine all parent generations
    for (li.gen in 1:i.parentgens) {
        list.parentpop <- list(fit=c(list.parentpop$fit,struc.ea$generations[[gen.parent[li.gen]]]$fit),
                               allele=rbind(list.parentpop$allele, struc.ea$generations[[gen.parent[li.gen]]]$allele))
    }

                                        # Make order vector
    vi.ordering <- order(list.parentpop$fit)
    i.parents <- length(list.parentpop$fit)


                                        # If n.ind is not given, choose average of all parent generations
    if (!length(n.ind))
        n.ind<-floor(i.parents/i.parentgens)


    list.parentgen <- NULL
    
                                        # Choose first 'nind' individiuals
    for(li.parent in 1:n.ind) {
        list.parentgen <- list(fit=c(list.parentgen$fit, list.parentpop$fit[vi.ordering[i.parents-li.parent+1]]),
                               allele=rbind(list.parentgen$allele, list.parentpop$allele[vi.ordering[i.parents-li.parent+1],]))
    }


                                        # Add new generation to object
    return(list(genes=struc.ea$genes,
                generations=.add.gen(struc.ea$generation,
                list(fit=list.parentgen$fit,
                     allele=list.parentgen$allele,
                     parents=as.character(gen.parent)))))
}



selectparents<-function(struc.ea,
                        gen=NULL,
                        method=list(base="fit", rescale=0)) {

    if (!length(gen))
        gen<-length(struc.ea$generations)

    list.parents <- NULL
    i.parents <- length(struc.ea$generations[[gen]]$fit)


                                        # Ranking base selection procedure
    rankbased<-function() {
        vd.ranking <- rank(struc.ea$generations[[gen]]$fit)
        vd.roulette <- ranking[1]

        if (method$rescale > 1) {
            d.maxfit <- max(d.ranking)
            d.avfit <- sum(d.ranking)/length(d.ranking)
            if (d.maxfit-d.avfit > 1E-12) {
                i.parent <- 1
                while (vd.ranking[i.parent] > d.avfit) {
                    vd.ranking[i.parent] <- d.avfit +
                        (vd.ranking[i.parent]-d.avfit) * (method$rescale*d.maxfit-d.avfit)/(d.maxfit-d.avfit)
                    i.parent <- i.parent+1
                }
            }
        }
    
        for (li.parent in 2:i.parents) {
            vd.roulette <- c(vd.roulette, vd.roulette[li.parent-1] + vd.ranking[li.parent])
        }
        return(vd.roulette)
    }


                                        # Fitness based selection procedure, optionally using rescaling
    fitbased<-function() {
        vd.roulette<-struc.ea$generations[[gen]]$fit[1]
        if (method$rescale > 1) {
            d.maxfit <- max(struc.ea$generations[[gen]]$fit)
            d.avfit <- sum(struc.ea$generations[[gen]]$fit)/length(struc.ea$generations[[gen]]$fit)
            if (d.maxfit-d.avfit > 1E-12) {
                i.parent <- 1
                while (struc.ea$generations[[gen]]$fit[i.parent] > d.avfit) {
                    struc.ea$generations[[gen]]$fit[i.parent] <- d.avfit +
                        (struc.ea$generations[[gen]]$fit[i.parent]-d.avfit) * (method$rescale*d.maxfit-d.avfit)/(d.maxfit-d.avfit)
                    i.parent <- i.parent+1
                }
            }
        }

        for (li.parent in 2:i.parents) {
            vd.roulette<-c(vd.roulette, vd.roulette[li.parent-1] + struc.ea$generations[[gen]]$fit[li.parent])
        }
        return(vd.roulette)
    }


    d.maxfit <- max(struc.ea$generations[[gen]]$fit)
    d.minfit <- min(struc.ea$generations[[gen]]$fit)
    
    if (d.maxfit-d.minfit < 1E-12) {
        vd.weights <- c(1:length(struc.ea$generations[[gen]]$fit))
    } else {
        switch(method$base,
               rank=vd.weights <- rankbased(),
               fit=vd.weights <- fitbased()
               )
    }

        
                                        # Select parents using roulette system
    for (li.parent in 1:i.parents) {
        d.random <- runif(1)*max(vd.weights)
        i.testparent <- 1
        while (vd.weights[i.testparent] <= d.random)
            i.testparent <- i.testparent + 1
        list.parents<-list(fit=c(list.parents$fit, struc.ea$generations[[gen]]$fit[i.testparent]),
                           allele=rbind(list.parents$allele, struc.ea$generations[[gen]]$allele[i.testparent,]))
    }
 
    struc.ea$generations[[gen]]$fit <- list.parents$fit
    struc.ea$generations[[gen]]$allele <- list.parents$allele
    struc.ea$generations[[gen]]$selection <- method
    
    return(struc.ea)
}



crossover<-function(struc.ea,
                    gen=NULL,
                    rate=90) {

    if (!length(gen))
        gen <- length(struc.ea$generations)
    
    i.parents <- length(struc.ea$generations[[gen]]$fit)

    while (rate >= 1)
        rate <- rate/100

    i.crossover <- 0

    for (li.parent in 1:floor(i.parents/2)) {
        if (runif(1) <= rate) {
            i.position <- ceiling(runif(1)*(length(struc.ea$genes)-1))
            vi.buffer <- struc.ea$generations[[gen]]$allele[2*li.parent-1, 1:i.position]
            struc.ea$generations[[gen]]$allele[2*li.parent-1, 1:i.position] <-
                struc.ea$generations[[gen]]$allele[2*li.parent, 1:i.position]
            struc.ea$generations[[gen]]$allele[2*li.parent, 1:i.position]<-vi.buffer
            i.crossover<-i.crossover+1
        }
    }

    struc.ea$generations[[gen]]$fit<-c(rep(NA, i.parents))
    d.rate <- i.crossover/(i.parents/2)

    struc.ea$generations[[gen]]$crossover <- d.rate

    return(struc.ea)
}




mutation<-function(struc.ea,
                   gen=NULL,
                   method=list(base="unif", spread=1, rate=15)) {

    if (!length(gen))
        gen<-length(struc.ea$generations)

    
    normalmutation<-function(spread=1) {

        if (!spread )
            calcallele<-function(spread=NULL) {
                spread<-1/6
                d.step<-(struc.ea$genes[[gene]]$max - struc.ea$genes[[gene]]$min)/
                    (struc.ea$genes[[gene]]$n.alleles - 1)
                floor(rnorm(1,
                            (struc.ea$generations[[gen]]$allele[parent, gene] - struc.ea$genes[[gene]]$min)/d.step,
                            (spread*struc.ea$genes[[gene]]$n.alleles)
                            ) * d.step) +
                     struc.ea$genes[[gene]]$min
            }
        else
            calcallele<-function(spread=NULL) {
                spread<-spread/6
                d.step<-(struc.ea$genes[[gene]]$max - struc.ea$genes[[gene]]$min)/
                    (struc.ea$genes[[gene]]$n.alleles - 1)
                floor(rnorm(1,
                            (struc.ea$generations[[gen]]$allele[parent, gene] - struc.ea$genes[[gene]]$min)/d.step,
                            (spread*struc.ea$genes[[gene]]$n.alleles)
                            ) * d.step) +
                     struc.ea$genes[[gene]]$min
            }
        return(calcallele)
    }

        
    uniformmutation<-function(spread=1) {
        if (! spread )
            calcallele<-function(spread=NULL) {
                floor(struc.ea$genes[[gene]]$n.alleles * runif(1)) *
                     (struc.ea$genes[[gene]]$max - struc.ea$genes[[gene]]$min)/
                     (struc.ea$genes[[gene]]$n.alleles - 1) +
                     (struc.ea$genes[[gene]]$min)
            }
        else
            calcallele<-function(spread=NULL) {
                floor(struc.ea$genes[[gene]]$n.alleles * (runif(1)-0.5) * 2 * spread) *
                     (struc.ea$genes[[gene]]$max - struc.ea$genes[[gene]]$min)/
                     (struc.ea$genes[[gene]]$n.alleles - 1) +
                     (struc.ea$generations[[gen]]$allele[parent, gene])
            }
        return(calcallele)
    }



    i.parents<-length(struc.ea$generations[[gen]]$fit)

    while (method$rate >= 1)
        method$rate <- method$rate/100

    
    if (length(method$spread))
        while (method$spread > 1)
            method$spread<-method$spread/100
    

    switch(method$base,
           unif = calcallele<-uniformmutation(method$spread),
           norm = calcallele<-normalmutation(method$spread)
           )


    i.mutation<-0
    
    for (parent in 1:i.parents) 
    for (gene in 1:length(struc.ea$genes))
    if (runif(1) <= method$rate) {
        d.randomallele<-calcallele(method$spread)
        while (
               (d.randomallele > struc.ea$genes[[gene]]$max) ||
               (d.randomallele < struc.ea$genes[[gene]]$min) ||
               (abs(d.randomallele - struc.ea$generations[[gen]]$allele[parent, gene]) < 1E-12)
               ) 
            d.randomallele<-calcallele(method$spread)
                        
        struc.ea$generations[[gen]]$allele[parent, gene]<-d.randomallele
        i.mutation<-i.mutation+1
    }

    d.rate<-i.mutation/(i.parents*length(struc.ea$genes))

    struc.ea$generations[[gen]]$mutation <- list(base=method$base, spread=method$spread, rate=d.rate)

    return(struc.ea)
}


nextgen <- function(struc.ea,
                    n.ind=25,
                    gen.parent=NULL,
                    selection=list(base="fit", rescale=0),
                    corate=90,
                    mutation=list(base="unif", spread=1, rate=15)) {
    
    struc.ea <- chooseparents(struc.ea, n.ind, gen.parent)
    
    struc.ea <- selectparents(struc.ea, length(struc.ea$generations), selection)
    
    struc.ea <- crossover(struc.ea, length(struc.ea$generations), corate)
    
    struc.ea <- mutation(struc.ea, length(struc.ea$generations), mutation)
    
    return(struc.ea)
}




plotevolution.fit<-function(struc.ea, gens=NULL, show="average", plot=TRUE, ...) {

    x<-NULL
    y<-NULL

    if (!length(gens))
        gens<-c(1:length(struc.ea$generations))
    
    for (i in gens) {
        x<-c(x, rep(i, length(struc.ea$generations[[i]]$fit)))
        y<-c(y, struc.ea$generations[[i]]$fit)
    }

    if (!plot)
        return (list(gens=gens,
                     min=tapply(y, x, min),
                     average=tapply(y, x, mean),
                     max=tapply(y, x, max),
                     stdev=tapply(y, x, function(z) {sqrt(var(z))})))
    else
        switch(show,
               average = plot(gens,
                         tapply(y, x, mean), ...),
               max = plot(gens,
                          tapply(y, x, max), ...),
               min = plot(gens,
                          tapply(y, x, min), ...),
               all = plot(x, y, ...),
               stdev = plot(gens,
                            tapply(y, x, function(z) {sqrt(var(z))} ), ... ),
               box = boxplot(split(y, x), ...),
               )
}


plotfreq.allele<-function(struc.ea, gene=1, gen=NULL, breaks=NULL, ...) {
    
    if (!length(gen))
        gen<-length(struc.ea$generations)


    allelevalues<-struc.ea$generations[[gen]]$allele[,gene]
   

    step<-(struc.ea$genes[[gene]]$max - struc.ea$genes[[gene]]$min)/(struc.ea$genes[[gene]]$n.alleles - 1)

    
    if (!length(breaks))
        hist(allelevalues,
             breaks=c(struc.ea$genes[[gene]]$min - 1.5*step + 1:(struc.ea$genes[[gene]]$n.alleles+1)*step),
             axes=F, ...)
    else
        hist(allelevalues,
             axes=F, ...)
    
    if (length(struc.ea$genes[[gene]]$namevector) > 0 )
        axis(1, at=c(struc.ea$genes[[gene]]$min + 0:(struc.ea$genes[[gene]]$n.alleles-1)*step),
             labels=struc.ea$genes[[gene]]$namevector)
    else
        axis(1)
    axis(2)
}




plotevolution.gene<-function(struc.ea, gene=1, gens=NULL, breaks=NULL, main="Gen", cols=NULL, rows=NULL, ...) {

    if (!length(gens))
        gens <- c(1:length(struc.ea$generations))

    i.depth<-length(gens)

    if (!length(rows)) {
        if ((i.depth %% 4 == 0) && (i.depth > 14)) {
            r <- 4
        } else if (i.depth %% 3 == 0) {
            r <- 3
        } else if (i.depth %% 2 == 0) {
            r <- 2
        } else {
            r <- floor(log(i.depth, base=2))
        }
    } else {
        r <- rows
    }
    
    if (!length(cols))
        c <- ceiling(i.depth/r)
    else {
        c <- cols
        r <- ceiling(i.depth/c)
    }


    l<-layout(matrix(c(1:i.depth, rep(0, c*r - i.depth)), r, c, byrow=T))

    

    for (d in gens) {
        plotfreq.allele(struc.ea,
                       gene,
                       d,
                       breaks,
                       main=paste(main, d), ...)
    }
}
