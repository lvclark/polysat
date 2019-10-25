# Non-iterative function for calculating allele frequencies
# under polysomic inheritance.  For allele copy number ambiguity,
# assumes that all alleles have an equal chance of being present in more
# than one copy.
simpleFreq <- function(object, samples=Samples(object), loci=Loci(object)){
  # subset the object to avoid a lot of indexing later in the function
  object <- object[samples,loci]

    # we need PopInfo and Ploidies; check for these
    if(!all(!is.na(PopInfo(object)))){
        cat("PopInfo needed for samples:",
            samples[is.na(PopInfo(object))], sep="\n")
        stop("PopInfo needed for this function.")
    }
    if(!all(!is.na(Ploidies(object)))){
        stop("Ploidies needed for this function.")
    }

    # we need a genbinary object for this function; make one if necessary
    if(is(object, "genambig")){
        object <- genambig.to.genbinary(object)
    }
    Present(object) <- as.integer(1)
    Absent(object) <- as.integer(0)

    # get the populations that will be evaluated
    pops <- PopNames(object)[unique(PopInfo(object))]

    # Get ploidies into a format where you can total over the populations
    # (and don't index by locus if you don't need to)
    if(is(object@Ploidies, "ploidymatrix")){
      object <- reformatPloidies(object, output="collapse", erase=FALSE)
    }
    if(!is(object@Ploidies, "ploidysample") &&
       length(unique(Ploidies(object, samples, loci)))==1){
      object <- reformatPloidies(object, output="sample", erase=FALSE)
    }
    if(is(object@Ploidies, "ploidylocus")){
      object <- reformatPloidies(object, output="matrix", erase=FALSE)
    }

    # Get total genomes per population if loci are uniform ploidy
    if(is(object@Ploidies, "ploidysample")){
      # get the total number of genomes per population
      totgenomes <- rep(0, length(pops))
      names(totgenomes) <- pops
      for(p in pops){
          totgenomes[p] <- sum(Ploidies(object)[Samples(object,populations=p)])
      }

          # set up data frame to contain allele frequencies
      freqtable <- data.frame(Genomes=totgenomes,row.names=pops)
    } else { # or if loci have different ploidies
      freqtable <- data.frame(row.names=pops)
    }

    # loop to get frequency data
    for(L in loci){
        # get all samples without missing data at this locus
        xsamples <- samples[!isMissing(object, samples, L)]
        # get the total number of genomes per population
        totgenomes <- rep(0, length(pops))
        names(totgenomes) <- pops
        for(p in pops){
            totgenomes[p] <- sum(Ploidies(object, xsamples[xsamples %in%
                                                           Samples(object,
                                                           populations=p)], L))
        }
        # write genomes to the table if necessary
        if(is(object@Ploidies, "ploidymatrix")){
          totg <- data.frame(totgenomes)
          names(totg) <- paste(L, "Genomes", sep=".")
          freqtable <- cbind(freqtable, totg)
        }
        # make a conversion factor to weight allele presence based on ploidy
        # of each individual and number of alleles at this locus
        numalleles <- rep(0, length(xsamples))
        names(numalleles) <- xsamples
        for(s in xsamples){
            numalleles[s] <- sum(Genotype(object, s , L))
        }
        convf <- as.vector(Ploidies(object,xsamples,L))/numalleles
        # make table of weighted allele presence
        loctable <- Genotypes(object, xsamples, L) * convf
        # loop through all alleles at this locus
        for(al in dimnames(loctable)[[2]]){
            theseallelefreqs <- rep(0, length(pops))
            names(theseallelefreqs) <- pops
            # loop through populations
            for(p in pops){
                theseallelefreqs[p]<-sum(loctable[xsamples[xsamples %in%
                                                Samples(object, populations=p)],
                                                  al])/totgenomes[p]
            }
            freqtable <- cbind(freqtable,theseallelefreqs)
            names(freqtable)[length(freqtable)] <- al
        }
    }

    # return data frame
    return(freqtable)
}

# Iterative estimation of allele frequencies under polysomic inheritance,
# with a uniform, even-numbered ploidy and a known selfing rate.
# Much of this code is translated directly from:
# De Silva HN, Hall AJ, Rikkerink E, McNeilage MA, and LG Fraser (2005).
# Estimation of allele frequencies in polyploids under certain patterns
# of inheritance.  Heredity 95:327-334.
deSilvaFreq <- function(object, self,
                        samples=Samples(object),
                        loci=Loci(object), initNull = 0.15,
                        initFreq=simpleFreq(object[samples,loci]),
                        tol = 0.00000001, maxiter = 1e4){
  # make sure self argument has been given
  if(missing(self)) stop("Selfing rate required.")
  
  # make sure initFreq is in right format
  if(!"Genomes" %in% names(initFreq))
    stop("initFreq must have single Genomes column.")
  
  # convert object to genambig if necessary
  if("genbinary" %in% class(object)) object <- genbinary.to.genambig(object,
                                                                     samples,
                                                                     loci)
  
  # get the ploidy (m2), and check that there is only one and that it is even
  m2 <- unique(as.vector(Ploidies(object,samples,loci)))
  if(length(m2) != 1)
    stop("Only one ploidy allowed.  Try running subsets of data one ploidy at a time.")
  if(is.na(m2)) stop("Function requires information in Ploidies slot.")
  if(m2 %% 2 != 0) stop("Ploidy must be even.")
  
  # check and set up initNull
  if(!length(initNull) %in% c(1,length(loci)))
    stop("Need single value for initNull or one value per locus.")
  if(length(initNull) == 1)
    initNull <- rep(initNull, times=length(loci))
  if(is.null(names(initNull)))
    names(initNull) <- loci
  
  # calculate the number of alleles from one gamete
  # (from de Silva et al's INIT subroutine)
  m <- m2 %/% 2L
  
  # get the populations that will be used
  pops <- PopNames(object)[unique(PopInfo(object)[samples])]
  if(!identical(pops, row.names(initFreq)))
    stop("Population names in initFreq don't match those in object.")
  
  # set up the data frame where final allele frequencies will be stored
  # has columns from initFreq, plus columns for nulls
  finalfreq <- data.frame(row.names=pops, Genomes=initFreq$Genomes,
                          matrix(0, nrow=length(pops),
                                 ncol=dim(initFreq)[2]-1+length(loci),
                                 dimnames=list(NULL, sort(c(
                                   paste(loci, ".null", sep=""),
                                   names(initFreq)[2:dim(initFreq)[2]])))))
  
  # get the divisor for the selfing matrices
  smatdiv <- (G(m-1,m+1))^2
  
  # INDEXF function from de Silva et al
  # af1 = vector of alleles in phenotype
  # m = number of observable alleles in phenotype
  # na = number of alleles, calculated for each locus
  indexf <- function(m, af1, na){
    x <- 1L
    if(m == 1){
      x <- x + af1[1]
    } else {
      if(m > 1){
        for(q in 1:(m-1)){
          x <- x + G(q-1,na-q+1)
        }
        x <- x + G(m-1,na+1-m) - G(m-1,na+2-m-af1[1])
        if(m > 2){
          for(q in 1:(m-2)){
            x <- x + G(q,na-q-af1[m-q-1]) -
              G(q,na+1-q-af1[m-q])
          }
        }
        x <- x + af1[m] - af1[m-1]
      }
    }
    return(x)
  }
  
  # FENLIST function
  fenlist <- function(na){
    # set up temporary holding vector for phenotpes (FENLIST)
    af1 <- rep(0L, m2)
    # set up array to contain all phenotypes (FENLIST)
    af <- array(0L, dim=c(1, m2))
    # set up vector for number of alleles in each phenotype (FENLIST)
    naf <- integer(0)
    
    # fill in af and naf (FENLIST)
    # This is done with rbind rather than setting up whole array,
    # because the method for calculating the number of phenotypes
    # doesn't seem to work for certain numbers of alleles.
    f <- 1L
    naf <- 0L
    for(m in 1:min(m2, na)){
      af1[1] <- 1L
      if(m > 1){
        for(j in 2:m){
          af1[j] <- af1[j-1] + 1L
        }
      }
      f <- f + 1L
      naf[f] <-  m
      af <- rbind(af, af1)
      a <- m
      while(a > 0){
        if(af1[a] == (na+a-m)){
          a <- a - 1L
        } else {
          if(a > 0){
            af1[a] <- af1[a] + 1L
            if(a < m){
              for(a1 in (a+1):m) af1[a1] <- af1[a1-1] + 1L
            }
            f <- f + 1L
            naf[f] <- m
            af <- rbind(af, af1)
            a <- m
          }
        }
      }
    }
    
    return(list(af, naf))
  }
  
  # CONVMAT function
  convmat <- function(ng, nf, na1, ag){
    na <- na1 - 1L
    af1 <- rep(0L, m2)
    # CONVMAT subroutine - make a matrix for conversion from
    # genotypes to phenotypes
    cmat <- matrix(0L, nrow=nf, ncol=ng)
    for(g in 1:ng){
      ag1 <- ag[g,]
      # find the phenotype (af1) that matches this genotype (ag1)
      if(ag1[1] == na1){
        naf1 <- 0L # this is the homozygous null genotype
      } else {
        naf1 <- 1L
        af1[naf1] <- ag1[1]
        for(a in 2:m2){
          if(ag1[a] == na1) break # exit loop if null allele
          if(ag1[a] > af1[naf1]){
            naf1 <- naf1 + 1L
            af1[naf1] <- ag1[a]
          }
        }
      }
      # fill in the extra alleles with zeros
      if(naf1 < m2){
        for(j in (naf1 + 1):m2){
          af1[j] <- 0L
        }
      }
      f <- indexf(naf1, af1, na) # This is the one phenotype to match
      # this genotype.
      cmat[f,g] <- 1L
    }
    return(cmat)
  }
  
  # Loop to do calculations one locus and population at a time
  for(L in loci){
    cat(paste("Starting", L), sep="\n")
    
    ## Begin looping through populations
    for(pop in pops){
      cat(paste("Starting", L, pop), sep="\n")
      # get a list of samples in this pop being analyzed
      psamples <- Samples(object, populations=pop)[!isMissing(object,
                                                              Samples(object, populations=pop),
                                                              L)]
      psamples <- psamples[psamples %in% samples]
      
      # extract the initial allele frequencies and add a null
      subInitFreq <- initFreq[pop, grep(paste("^", L,"\\.",sep=""),
                                        names(initFreq))]
      
      # set up matrices if this has not already been done
      templist <- names(subInitFreq)[subInitFreq !=0]
      templist <- strsplit(templist, split=".", fixed=TRUE)
      alleles <- rep(0L, length(templist))
      for(i in 1:length(alleles)){
        alleles[i] <- as.integer(templist[[i]][2])
      }
      alleles <- sort(alleles)
      na <- length(alleles)
      
      # get the number of alleles with a null
      na1 <- na + 1L
      # get the number of genotypes (from the GENLIST function)
      ng <- na1
      for(j in 2:m2){
        ng <- ng * (na1 + j - 1L) / j
      }
      
      ag <- GENLIST(ng, na1, m2)
      temp <- fenlist(na)
      af <- temp[[1]]
      naf <- temp[[2]]
      nf <- length(naf)
      temp <- RANMUL(ng, na1, ag, m2)
      rmul <- temp[[1]]
      arep <- temp[[2]]
      smatt <- SELFMAT(ng, na1, ag, m2)/smatdiv
      cmat <- convmat(ng, nf, na1, ag)
      rm(temp)
      
      # calculate pp, the frequency of each phenotype in this population
      pp <- rep(0, nf)
      for(s in psamples){
        phenotype <- sort(unique(Genotype(object, s, L)))
        phenotype <- match(phenotype, alleles)
        f <- indexf(length(phenotype), phenotype, na)
        pp[f] <- pp[f] + 1
      }
      pp <- pp/sum(pp)
      
      # Get initial allele frequencies
      p1 <- rep(0, na1) # vector to hold frequencies
      p1[na1] <- initNull[L] # add null freq to the last position
      # make sure everything will sum to 1
      subInitFreq <- subInitFreq * (1 - initNull[L])/sum(subInitFreq)
      # get each allele frequency
      for(a in alleles){
        p1[match(a, alleles)] <- subInitFreq[1,paste(L, a, sep = ".")]
      }
      
      ## Begin the EM algorithm
      converge <- 0
      niter <- 1L
      oneg <- rep(1, ng)
      while(converge == 0){
        # Expectation step
        # GPROBS subroutine
        pa <- rep(0, na1)
        pa[na1] <- 1
        for(j in 1:na){
          pa[j] <- p1[j]
          pa[na1] <- pa[na1]-p1[j]
          # it seems like this is the same as pa <- p1
          # unless p1 is not supposed to have the null allele
          # (conflicting information in orignal code comments)
        }
        rvec <- rep(0, ng)
        for(g in 1:ng){
          rvec[g] <- rmul[g]
          for(j in 1:m2){
            rvec[g] <- rvec[g]*pa[ag[g,j]]
          }
          # This gives prob of genotype by multiplying
          # probs of alleles and coefficients for a multinomial
          # distribution.
        }
        
        # make an identity matrix
        id <- diag(nrow=ng)
        # calculate gprob using selfing and outcrossing matrices
        s3 <- id - self * smatt
        gprob <- (1-self) * solve(s3, rvec)
        # end GPROBS
        
        # equation (12) from the paper
        xx1 <- matrix(0, nrow=nf, ncol=ng)
        for(i in 1:nf){
          xx1[i,] <- cmat[i,] * gprob
        }
        xx2 <- xx1 %*% oneg
        xx3 <- matrix(0, nrow=nf, ncol=ng)
        for(i in 1:ng){
          xx3[,i] <- xx1[,i] * (xx2^(-1))
        }
        EP <- t(xx3) %*% pp
        
        # Maximization step
        p2 <- t(arep) %*% EP/m2
        
        # check for convergence
        pB <- p1 + p2
        pT <- p1 - p2
        pT <- pT[abs(pB) > 1e-14]
        pB <- pB[abs(pB) > 1e-14]
        
        if(length(pB) == 0){
          converge <- 1
        } else {
          if(sum(abs(pT)/pB) <= tol){
            converge <- 1
          }
        }
        if(niter >= maxiter){
          converge <- 1
        }
        
        niter <- niter + 1L
        p1 <- p2
      }
      
      # write frequencies in p2 to finalfreqs
      for(a in alleles){
        finalfreq[pop, match(paste(L, a, sep="."), names(finalfreq))]<-
          p2[match(a, alleles)]
      }
      finalfreq[pop, match(paste(L, "null", sep="."), names(finalfreq))]<-
        p2[na1]
      # print the number of reps
      cat(paste(niter-1, "repetitions for", L, pop),
          sep="\n")
    }
  }
  return(finalfreq)
}

calcPopDiff<-function(freqs, metric, pops=row.names(freqs),
                  loci=unique(as.matrix(as.data.frame(strsplit(names(freqs),
                  split=".",
                                                      fixed=TRUE),
                                             stringsAsFactors=FALSE))[1,]), 
                     global = FALSE, bootstrap = FALSE, n.bootstraps = 1000,
                  object = NULL){
  # check metric
  if(!metric %in% c("Fst", "Gst", "Jost's D", "Rst")){
    stop("metric must be Fst, Gst, Rst, or Jost's D")
  }
  
  # check pop names
  if(!all(pops %in% row.names(freqs))){
    stop("pops must all be in row names of freqs.")
  }
  freqs <- freqs[pops,]
  # Clean up loci
  loci<-loci[loci!="Genomes"]
  
  if(metric == "Rst" && (is.null(object) || any(is.na(Usatnts(object)[loci])))){
    stop("gendata object required with Usatnts for Rst metric")
  }
  
  # internal function to get differentiation statistic from any set of two or more pops
  cpd <- function(freqs, metric, loci, bootstrap, n.boostraps){
    # Get genome number from the table
    if("Genomes" %in% names(freqs)){
      genomes<-freqs$Genomes
      names(genomes)<-row.names(freqs)
      GbyL <- FALSE
    } else {
      GbyL <- TRUE
    }
    # set up array for HT and HS values
    hets<-array(0,dim=c(length(loci),4),
                dimnames=list(loci,c("HT","HS","HTest","HSest")))
    for(L in loci){
      # population sizes for this locus
      if(!GbyL){
        thesegenomes <- genomes
      } else {
        thesegenomes <- freqs[,paste(L, "Genomes", sep=".")]
      }
      nsubpop <- length(thesegenomes)
      # allele frequencies for this locus
      thesefreqs<-freqs[,grep(paste("^", L,"\\.",sep=""), names(freqs)), drop = FALSE]
      thesefreqs <- thesefreqs[,names(thesefreqs)!=paste(L,"Genomes",sep="."), drop = FALSE]
      # expected heterozygosity by pop
      hsByPop <- apply(as.matrix(thesefreqs), 1, function(x) 1 - sum(x^2))
      if(metric == "Fst"){
        # weighted mean frequencies for each allele
        avgfreq <- unlist(lapply(thesefreqs, weighted.mean, w = thesegenomes))
        # weighted mean Hs
        hets[L, "HS"] <- weighted.mean(hsByPop, thesegenomes)
      }
      if(metric %in% c("Jost's D", "Gst")){
        # unweighted mean frequencies for each allele
        avgfreq <- colMeans(thesefreqs)
        # unweighted mean Hs
        hets[L, "HS"] <- mean(hsByPop)
      }
      if(metric == "Rst"){
        replen <- Usatnts(object)[L] # microsatellite repeat length
        alleles <- sapply(strsplit(grep(paste("^", L,"\\.", sep = ""), names(freqs), value = TRUE), ".",
                                   fixed = TRUE),
                          function(x) x[2])
        alleles <- alleles[alleles != "Genomes"]
        alleles[alleles == "null"] <- 0
        alleles <- as.integer(alleles)
        if(0 %in% alleles){ # remove null alleles if they are there
          nullfreqs <- thesefreqs[,match(0, alleles)] # frequency of null allele in each pop
          for(pop in 1:dim(thesefreqs)[1]){
            thesefreqs[pop,] <- thesefreqs[pop,]/(1-nullfreqs[pop]) # scale so non-null allele frequencies sum to 1
          }
          thesegenomes <- round(thesegenomes * (1-nullfreqs)) # remove null allele copies from total genomes
        }
        totgenomes <- sum(thesegenomes)
        avgfreq <- colMeans(thesefreqs)
        SSalleledistS <- numeric(dim(freqs)[1]) # for totaling sums of squares of allele differences for each pop
        SSalleledistT <- 0 # for totaling sums of squares of allele differences across all pops
        for(i in 1:(length(alleles)-1)){ # loop through all pairs of different alleles
          for(j in (i+1):length(alleles)){
            if(alleles[i] == 0 || alleles[j] == 0) next
            sqdiff <- (abs(alleles[i] - alleles[j])/replen)^2 # squared difference in repeat number
            nocc <- thesefreqs[,i] * thesefreqs[,j] * thesegenomes^2 # number of times these alleles would be compared
            SSalleledistS <- SSalleledistS + sqdiff * nocc
            SSalleledistT <- SSalleledistT + sqdiff * totgenomes^2 * avgfreq[i] * avgfreq[j]
          }
        }
        hets[L, "HS"] <- mean(SSalleledistS / (thesegenomes * (thesegenomes - 1)))
        hets[L, "HT"] <- SSalleledistT/(totgenomes * (totgenomes - 1))
      }
      # estimate expected heterozygosity for panmictic pop
      if(metric %in% c("Fst", "Gst", "Jost's D")){
        hets[L,"HT"]<-1-sum(avgfreq^2)
      }
      if(metric %in% c("Jost's D","Gst")){
        # harmonic mean of the sample size
        meanGenomes <- 1/mean(1/thesegenomes)
        # Nei and Chesser's (1983) estimates of expected heterozygosity
        hets[L, "HSest"] <- hets[L, "HS"] * meanGenomes/(meanGenomes - 1)
        hets[L, "HTest"] <- hets[L, "HT"] + hets[L, "HSest"] / (nsubpop * meanGenomes)
      }
    }
    # set up vector to contain results, and bootstrapping
    if(!bootstrap){
      n.bootstraps <- 1
    }
    result <- numeric(n.bootstraps)
    for(b in 1:n.bootstraps){ # loop through bootstraps (just one if not bootstrapping)
      if(bootstrap){
        thishets <- hets[sample(loci, replace = TRUE),] # H matrix with bootstrapped set of loci
      } else {
        thishets <- hets
      }
      if(metric == "Fst"){ # calculate Fst
        HT<-mean(thishets[,"HT"])
        HS<-mean(thishets[,"HS"])
        result[b] <- (HT-HS)/HT
      }
      if(metric == "Rst"){
        R <- (thishets[,"HT"] - thishets[,"HS"])/thishets[,"HT"]
        result[b] <- mean(R)
      }
      if(metric == "Gst"){ # calculate Gst and average across loci
        G <- (thishets[,"HTest"] - thishets[,"HSest"])/thishets[,"HTest"]
        result[b] <- mean(G)
      }
      if(metric == "Jost's D"){ # calculate Jost's D and average across loci
        D <- nsubpop / (nsubpop - 1) * (thishets[,"HTest"] - thishets[,"HSest"]) / (1 - thishets[,"HSest"])
        result[b] <- mean(D)
        if(nsubpop == 1) result[b] <- 0 # prevent NaN
      }
    }
    
    return(result)
  } # end internal function
  
  # wrappers for internal function (global or pairwise)
  if(global){ # estimate single global statistic
    result <- cpd(freqs, metric, loci, bootstrap, n.bootstraps)
  } else { # estimate pairwise statistics
    # Set up matrix for pairwise diffentiation statistics
    if(bootstrap){
      result <- array(list(), dim = c(length(pops), length(pops)), dimnames = list(pops,pops)) # needs to be array-list if each item will be vector of bootstraps
    } else {
      result<-matrix(0,nrow=length(pops),ncol=length(pops),dimnames=list(pops,pops)) # matrix for single non-bootstrapped values
    }
    for(m in 1:length(pops)){
      for(n in m:length(pops)){
        thisres <- cpd(freqs[unique(c(m,n)),], metric, loci, bootstrap, n.bootstraps)
        if(bootstrap){
          result[[m, n]] <- result[[n, m]] <- thisres
        } else {
          result[m, n] <- result[n, m] <- thisres
        }
      }
    }
  }
  # return matrix of differentiation statistics (or single global value)
  return(result)
}

# wrapper function for backwards compatibility
calcFst<-function(freqs, pops=row.names(freqs),
                  loci=unique(as.matrix(as.data.frame(strsplit(names(freqs),
                                                               split=".",
                                                               fixed=TRUE),
                                                      stringsAsFactors=FALSE))[1,]), ...){
  return(calcPopDiff(freqs = freqs, metric = "Fst", pops = pops, loci = loci, ...))
}

# put allele frequencies into a format for SPAGeDi for it to use in estimates
# of kinship and relatedness coefficients
write.freq.SPAGeDi <- function(freqs, usatnts, file="", digits=2,
                               pops=row.names(freqs),
                  loci=unique(as.matrix(as.data.frame(strsplit(names(freqs),
                  split=".", fixed=TRUE), stringsAsFactors=FALSE))[1,])){
    if(is.null(names(usatnts))) names(usatnts) <- loci
    loci <- loci[loci != "Genomes"]
    # subset the data frame
    freqs <- freqs[pops,]
    # make a list to contain the columns to write
    datalist <- list(0)
    length(datalist) <- ( length(loci) * 2 )
    item <- 1
    genbyloc <- !"Genomes" %in% names(freqs)
    if(!genbyloc) genomes <- freqs$Genomes

    for(L in loci){
        # find the alleles
        alleles <- as.matrix(as.data.frame(strsplit(names(freqs)[grep(paste("^", L,"\\.",sep=""),
                                                            names(freqs))],
                                          split=".", fixed=TRUE))[2,])
        alleles <- as.vector(alleles)
        if(genbyloc){
          alleles <- alleles[alleles!="Genomes"]
        }
        # convert alleles to numbers used in SPAGeDi file
        alleles[alleles == "null"] <- 0
        alleles <- as.integer(alleles)
        alleles <- floor(alleles/usatnts[L])
        suballele<-function(value){if(value[1]!=0){
            value-(10^(digits-1))} else{
            0}}
        while(max(alleles) >= 10^digits){
            alleles<-mapply(suballele, alleles)
        }
        # add alleles to datalist
        datalist[[item]] <- alleles
        names(datalist)[item] <- L
        item <- item + 1
        # make a weighted average of allele frequencies across populations
        if(genbyloc){
          genomes <- freqs[[paste(L, "Genomes", sep=".")]]
          locindex <- grep(paste("^", L, "\\.", sep=""), names(freqs))
          locindex <- locindex[!locindex %in% grep("Genomes", names(freqs))]
        } else {
          locindex <- grep(paste("^", L, "\\.", sep=""), names(freqs))
        }
        avgfreq <- (genomes %*% as.matrix(freqs[,locindex])) / sum(genomes)
        # add frequencies to list
        datalist[[item]] <- as.vector(avgfreq)
        names(datalist)[item] <- length(avgfreq)
        item <- item + 1
    }

    # find the maximum number of alleles
    maxal <- max(mapply(length, datalist))
    # set up data frame to write
    fr <- data.frame(row.names=1:maxal)
    # put list elements into the data frame
    for(i in datalist){
        fr <- data.frame(fr, c(i, rep("", maxal-length(i))))
    }
    names(fr) <- names(datalist)

    # write the file
    write.table(fr, file=file, sep="\t", col.names=TRUE, row.names=FALSE,
                quote=FALSE)
}

freq.to.genpop <- function(freqs, pops=row.names(freqs),
            loci=unique(as.matrix(as.data.frame(strsplit(names(freqs),
                  split=".",fixed=TRUE),stringsAsFactors=FALSE))[1,])){
    # clean up loci
    loci <- loci[loci!="Genomes"]

    # errors
    if(!"Genomes" %in% names(freqs))
      stop("Only one Genomes column allowed.")

    # get population sizes
    popsizes <- freqs[pops, "Genomes"]

    # get an index of columns containing these loci
    loccol <- integer(0)
    for(L in loci){
        loccol <- c(loccol, grep(paste("^", L, "\\.",sep=""), names(freqs)))
    }
    # get a table just for these populations and loci
    subfreq <- freqs[pops, loccol]

    # convert allele frequencies to allele counts
    for(i in 1:length(subfreq)){
        subfreq[[i]] <- round(popsizes * subfreq[[i]])
    }

    return(subfreq)
}

## Genotype diversity functions
# p is a vector of genotype counts
# base = exp(1) for natural log, or base = 2 for log base 2
Shannon <- function(p, base=exp(1)){
    N <- sum(p)
    return(-sum(p/N * log(p/N, base=base)))
}

Simpson <- function(p){
    N <- sum(p)
    return(sum(p*(p-1)/(N*(N-1))))
}

Simpson.var <- function(p){
    N <- sum(p)
    p <- p/N
    v <- (4*N*(N-1)*(N-2)*sum(p^3) + 2*N*(N-1)*sum(p^2) -
            2*N*(N-1)*(2*N-3)*(sum(p^2))^2)/((N*(N-1))^2)
    return(v)
}

# ... is passed to index. (So you can adjust base of Shannon index.)
# threshold is highest distance between two individuals that can be considered
# to be one clone.
genotypeDiversity <- function(genobject,
                   samples = Samples(genobject), loci = Loci(genobject),
                    d = meandistance.matrix(genobject,samples,loci,
                    all.distances=TRUE,distmetric=Lynch.distance),
                    threshold = 0, index = Shannon,
                              ...){
    # get populations
    if(all(is.na(PopInfo(genobject)[samples]))){
        PopInfo(genobject)[samples] <- rep(1, length(samples))
    }
    if(!all(!is.na(PopInfo(genobject)[samples]))){
        stop("PopInfo needed.")
    }
    pops <- PopNames(genobject)[unique(PopInfo(genobject)[samples])]

    # set up results matrix
    results <- matrix(NA, ncol=length(loci)+1, nrow=length(pops),
                      dimnames=list(pops, c(loci,"overall")))

    # statistics for individual loci
    for(L in loci){
        for(p in pops){
            # get all samples for this pop with data for this locus
            psamples <- samples[samples %in% Samples(genobject, populations=p)]
            psamples <- psamples[!isMissing(genobject, psamples, L)]

            # assign to genotype groups
            dsub <- d[[1]][L, psamples, psamples]
            if(is.vector(dsub)){ # fix for if there was only one sample
                dsub <- array(dsub, dim=c(1,1),
                              dimnames=list(psamples, psamples))
            }
            clones <- assignClones(dsub,
                                   threshold=threshold)
            # get genotype frequencies
#            n <- length(clones) # number of individuals
            cl <- length(unique(clones)) # number of groups
            counts <- rep(NA, cl) # vector to hold counts of individuals
            for(i in 1:cl){
                counts[i] <- length(clones[clones == i])
            }
            # get diversity index
            results[p,L] <- index(counts, ...)
        }
    }

    # get a mean array for only the loci being used here
    if(length(loci) > dim(d[[1]])[1]){
        d2 <- meandist.from.array(d[[1]], samples=samples, loci=loci)
    } else {
        d2 <- d[[2]]
    }
    # statistics for multilocus genotypes
    for(p in pops){
        psamples <- samples[samples %in% Samples(genobject, populations=p)]
        clones <- assignClones(d2, psamples, threshold)
#        n <- length(clones) # number of individuals
        cl <- length(unique(clones)) # number of groups
        counts <- rep(NA, cl) # vector to hold counts of individuals
        for(i in 1:cl){
            counts[i] <- length(clones[clones == i])
        }
        # get diversity index
        results[p,"overall"] <- index(counts, ...)
    }

    return(results)
}

# function to get unique alleles and counts
alleleDiversity <- function(genobject, samples=Samples(genobject),
                            loci=Loci(genobject),
                            alleles=TRUE, counts=TRUE){
  if(!alleles && !counts)
    stop("At least one of alleles and counts must be set to TRUE")
  # get populations
    if(all(is.na(PopInfo(genobject)[samples]))){
        PopInfo(genobject)[samples] <- rep(1, length(samples))
    }
    if(!all(!is.na(PopInfo(genobject)[samples]))){
        stop("PopInfo needed.")
    }
    pops <- PopNames(genobject)[unique(PopInfo(genobject)[samples])]

    # set up results tables
    rcounts <- matrix(NA, nrow=length(pops)+1, ncol=length(loci),
                     dimnames=list(c(pops,"overall"),loci))
    ralleles <- array(list(NA), dim=c(length(pops)+1,length(loci)),
                      dimnames=dimnames(rcounts))
    # get unique alleles
    for(L in loci){
      for(p in pops){
        psamples <- samples[samples %in% Samples(genobject, populations=p)]
        xalleles <- .unal1loc(genobject, psamples, L)
        ralleles[[p,L]] <- xalleles
        rcounts[p,L] <- length(xalleles)
      }
      xalleles <- .unal1loc(genobject, samples, L)
      ralleles[["overall",L]] <- xalleles
      rcounts["overall",L] <- length(xalleles)
    }

    if(alleles && counts)
      return(list(alleles=ralleles, counts=rcounts))
  if(alleles && !counts)
    return(ralleles)
  if(!alleles && counts)
    return(rcounts)
}

PIC <- function(freqs, pops=row.names(freqs),
                loci=unique(as.matrix(as.data.frame(strsplit(names(freqs),
                                                             split=".",
                                                             fixed=TRUE),
                                                    stringsAsFactors=FALSE))[1,]),
                bypop = TRUE, overall = TRUE){
  if(!bypop && !overall){
    stop("At least one of bypop and overall must be TRUE.")
  }
  # check pop names
  if(!all(pops %in% row.names(freqs))){
    stop("pops must all be in row names of freqs.")
  }
  freqs <- freqs[pops,]
  # Clean up loci
  loci<-loci[loci!="Genomes"]
  # is number of genomes (pop size) specified by locus?
  GbL <- !"Genomes" %in% names(freqs)
  # convert freqs to matrix for math
  freqs <- as.matrix(freqs)
  
  # set up results matrix
  if(bypop){
    results <- matrix(NA, nrow = length(pops), ncol = length(loci), dimnames = list(pops, loci))
  } else {
    results <- matrix(NA, nrow = 0, ncol = length(loci), dimnames = list(character(0), loci))
  }
  if(overall){
    results <- rbind(results, matrix(NA, nrow = 1, ncol = length(loci), dimnames = list("Overall", loci)))
  }
  
  # function to get PIC for one locus and pop
  pic <- function(fr){
    if(any(is.na(fr))) return(NA)
    if(!isTRUE(all.equal(sum(fr), 1))) stop("Allele frequencies don't sum to 1.")
    sq <- fr^2 # square each frequncy to get p_i^2 and p_j^2
    nAl <- length(fr) # number of alleles
    # matrix of the squares multiplied by each other
    mt <- matrix(sq, nrow = nAl, ncol = 1) %*% matrix(sq, nrow = 1, ncol = nAl)
    # sum of the squared allele freqs
    sum1 <- sum(sq)
    # twice the sum of the product of allele freqs (only between different alleles)
    sum2 <- sum(mt) - sum(diag(mt))
    return(1 - sum1 - sum2)
  }
  
  # loop through loci 
  for(L in loci){
    lcol <- grep(paste("^", L, "\\.", sep = ""), dimnames(freqs)[[2]]) # columns for this locus
    if(length(lcol) == 0) stop(paste("Locus", L, "not found in freqs."))
    if(GbL){
      lgencol <- grep(paste("^", L, "\\.Genomes", sep = ""), dimnames(freqs)[[2]])
      lcol <- lcol[lcol != lgencol]
    }
    if(bypop){ # get PIC for each population for this locus
      for(p in pops){
        results[p,L] <- pic(freqs[p,lcol])
      }
    }
    if(overall){ # get PIC for mean allele frequency
      if(GbL){
        genomes <- freqs[,lgencol]
      } else {
        genomes <- freqs[,"Genomes"]
      }
      meanfreq <- apply(freqs[,lcol], 2, function(x) weighted.mean(x, w = genomes))
      results["Overall",L] <- pic(meanfreq)
    }
  }
  
  return(results)
} # end of PIc function

