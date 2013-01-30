# load dependencies
require(ape)
require(geiger)
require(TreeSim)

#define parameters
nTaxa<-c(25,50,75,100)
numbsim<-1
lambda<-1
mu<-c(0,0.5)
age<-c(0.1,0.5,0.9)
numReps<-5
model<-c("JC","HKY","GTR","JC+G","HKY+G","GTR+G")
mastergrid<-expand.grid(nTaxa,mu,age)
names(mastergrid)<-c("nTaxa","mu","age")

#simulation parameters defined following Sullivan and Swofford (2001)
numSites<-"-l 1000"
freq<-("-f 0.3,0.25,0.15,0.3")
titv<-("-t 3.0")
alpha=("-a 1.0")
ratematrix<-("-r 0.5,10,3,1,15,1")

seqFilenames <- list()

getTreesSeqs<-function(nTaxa,numbsim,lambda,mu,age,numReps,numSites,freq,titv,alpha,ratematrix,model,seqFilenames){

	quiet=TRUE
	print(nTaxa)
	nextrees<-as.list(1:numReps)
	counter<-0
	while (counter<numReps){
		counter<-counter+1
		nextrees[[counter]]<-sim.bd.taxa.age(n=nTaxa,numbsim=1,mu=mu,frac=1,lambda=lambda,age=age,mrca=TRUE)
		
#		treeFileName=paste("b_",lambda,"_d_",mu,"_a_",age,"_n_",nTaxa,"_rep_",counter,".phy", sep="")
		
		write.tree(nextrees[[counter]][[1]],append=FALSE,file=paste("b_",lambda,"_d_",mu,"_a_",age,"_n_",nTaxa,"_rep_",counter,".phy", sep=""))
		for (x in 1:length(model)){
			seed<-abs(round(runif(1,max=100000)))
			print(counter)
			print(seed)
			print(model[x])
			if (model[x]=="GTR+G")
		{
			settings <- paste("-mGTR ", numSites," ", ratematrix," ", freq, " -on", " -z ", seed," ", alpha, sep="")
		} 
		
		else if (model[x] == "HKY+G")
		{
			settings <- paste("-mHKY ", numSites," ", titv," ", freq, " -on", " -z ", seed," ", alpha, sep="")
		} 
		
		else if (model[x] == "GTR")
		{
			settings<-paste("-mGTR ", numSites," ", ratematrix," ", freq, " -on", " -z ", seed, sep="")
		} 
		
		else if (model[x] == "HKY") 
		{
			settings<-paste("-mHKY ", numSites," ", titv," ", freq, " -on", " -z ", seed, sep="")
		} 
		
		else if (model[x]=="JC")
		{
			settings<-paste("-mHKY ", numSites," -t 0.5 -f 0.25,0.25,0.25,0.25 -on", " -z ", seed, sep="")
		} 
		
		else if (model[x]=="JC+G")
		{
			settings<-paste("-mHKY ", numSites," -t 0.5 -f 0.25,0.25,0.25,0.25 -on", " -z ", seed," ",alpha, sep="")
		}
		
		else {print("NOT A MODEL IMPLEMENTED")}
		
		#append quiet to eliminate verbose logging
		if (quiet) {settings <- paste(settings, " -q", sep="")}
		
		filename<-paste("b_",lambda,"_d_",mu,"_a_",age,"_n_",nTaxa,"_sim_",model[x],"_rep_", counter, ".NEX", sep="")
		cat(filename,"\n")
		seqFilenames[length(seqFilenames)+1] <- filename
		
		system(paste("./seq-gen",settings,"<",paste("b_",lambda,"_d_",mu,"_a_",age,"_n_",nTaxa,"_rep_", counter,".phy", sep=""),">",filename,sep=" "))
		}
	}
	return(seqFilenames)
}

# I've added the bit below so that filenames are logged; this can be passed into my code
# Use as:
# files <- runIt()
# write(files, file="files.txt")
runIt <- function()
{
	seqFilenames <- list()
	for (xx in 1:nrow(mastergrid))
	{
		nTaxanew<-mastergrid$nTaxa[xx]
		print(nTaxanew)
		munew<-mastergrid$mu[xx]
		print(munew)
		agenew<-mastergrid$age[xx]
		print(agenew)
		seqFilenames <- getTreesSeqs(nTaxa=nTaxanew,numbsim,lambda,mu=munew,age=agenew,numReps,numSites,freq,titv,alpha,ratematrix,model,seqFilenames)
	}
	return(unlist(seqFilenames))
}