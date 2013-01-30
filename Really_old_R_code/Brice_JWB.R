# load dependencies
require(phyclust)
require(ape)
require(geiger)
require(TreeSim)

#simulate trees
holder<-0
nTaxa<-25
numbsim<-1
lambda<-1
mu<-0.5
age<-0.9

# These don't seem to be used
lengths<-0
lengths<-as.list(lengths)

# I bundled this into a function; makes it easier to loop over diversification parameters, or whatever
getTreesSeqs <- function(nTaxa=25, numRep=10, numSites=1000, lambda=1, mu=0.5, numbsim=1, age=0.9, model=c("HKY", "GTR", "HKY+G", "GTR+G"), quiet=TRUE)
{
	trees <-list()
	sequen <-list()

	counter<-0
	
	while (counter<numRep){
	counter<-counter+1
	trees[counter]<-sim.bd.taxa.age(n=nTaxa, numbsim=numbsim, frac=1.0, lambda=lambda, mu=mu, age=age, mrca=TRUE)
#	holder<-sort(trees[[counter]]$edge.length, decreasing=TRUE)
	
	for (i in 1:length(model))
	{
		seed <- runif(1, max=10000000)
		freqs <- getNucFreq()
		
		if (model[i] == "GTR+G")
		{
			alpha <- rlnorm(1, meanlog = 1, sdlog = 1)
			rates <- toString(runif(6, max=5)) # You might want to think about acceptable ranges
			settings <- paste("-m", model[i], " -l", numSites, " -r", rates, " -f", freqs, " -on", " -z", seed, sep="")
		} else if (model[i] == "HKY+G") {
			alpha <- rlnorm(1, meanlog = 1, sdlog = 1)
			tratio <- runif(1, max=5) # You might want to think about acceptable ranges
			settings <- paste("-m", model[i], " -l", numSites, " -t", tratio, " -f", freqs, " -on", " -z", seed, sep="")
		} else if (model[i] == "GTR") {
			rates <- toString(runif(6, max=5)) # You might want to think about acceptable ranges
			settings <- paste("-m", model[i], " -l", numSites, " -r", rates, " -f", freqs, " -on", " -z", seed, sep="")
		} else if (model[i] == "HKY") {
			tratio <- runif(1, max=5) # You might want to think about acceptable ranges
			settings <- paste("-m", model[i], " -l", numSites, " -t", tratio, " -f", freqs, " -on", " -z", seed, sep="")
		}
		
		if (quiet) settings <- paste(settings, " -q", sep="")
		sequen[[counter]]<-seqgen(opts=settings, rooted.tree=trees[[counter]])
		write(sequen[[counter]], file=paste("b",lambda,"d",mu,"a",age,"n",nTaxa,".",model[i],".", counter, ".NEX", sep=""))
	}
	write.tree(trees[[counter]], file=paste("b",lambda,"d",mu,"a",age,"n",nTaxa,".", counter, ".phy", sep=""))
	
	cat(counter, "\n")
	}
}
#stores ONLY trees with a single-taxon root; useful for later
#	if (holder[1]==age){
#		trees[[counter+1]]<-tree
#		counter<-counter+1
#		print(counter)
#		}	
#	}

getNucFreq <- function() {
# NOTE: using a uniform range of 0.1 to 0.4; seems reasonable, but you may want to tweak it.
	cool <- FALSE
	while (!cool)
	{
		sum <- 0
		freq <- numeric()
		for (i in 1:3)
		{
			freq[i] <- runif(1, max=min(0.4, 1-sum))[1]
			sum <- sum + freq[i]
		}
		freq[4] <- 1 - sum
		if (any(freq > 0.5) || any(freq < 0.05)) cool <- FALSE # disregard any where a single nucleotide is > 50% or < 5%
		else cool <- TRUE
	}
# Don't want to bias results (say, always having a high A and low T content), so randomize frequencies
	rand <- runif(4)
	freq <- freq[order(rand)]
	
	return (toString(freq))
}

