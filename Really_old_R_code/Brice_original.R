# load dependencies
require(phyclust)
require(ape)
require(geiger)
require(TreeSim)

#simulate trees
counter<-0
holder<-0
trees<-0
tree<-0
trees<-as.list(trees)
n<-25
numbsim<-1
lambda<-1
mu<-0.5
age<-0.9
lengths<-0
lengths<-as.list(lengths)
while (counter<10){
	tree<-sim.bd.taxa.age(n=n, numbsim=numbsim, frac=1.0, lambda=lambda, mu=mu, age=age, mrca=TRUE)
	holder<-sort(tree[[1]]$edge.length, decreasing=TRUE)
	trees[[counter+1]]<-tree
	counter<-counter+1
	print(counter)
	}
#stores ONLY trees with a single-taxon root; useful for later
#	if (holder[1]==age){
#		trees[[counter+1]]<-tree
#		counter<-counter+1
#		print(counter)
#		}	
#	}

#write trees
#for (i in 1:counter){
#write.tree(trees[[i]][[1]], append=TRUE, file=paste("b",lambda,"d",mu,"a",age,"n",n,".phy", sep=""))
#	}

#generate trees given topology
sequen<-seqgen(opts="-mHKY -l1000 -on", rooted.tree=trees[[1]][[1]])
write<-(sequen, file="test.nex")
