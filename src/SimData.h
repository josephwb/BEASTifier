#ifndef _SIMDATA_H_
#define _SIMDATA_H_

class SimData {
	
	int numTaxa, numChar;
	vector < vector <string> > taxaAlignment;
	string seqFileName, root, simModel, treeFileName, starterTree;
	bool interleavedData, starterTreePresent;
	
public:
	void setNumTaxaChar (string&, int&, int&, bool&);
	vector < vector <string> > collectTaxaAlignment (string&, int const&, int const&, bool const&);
	string setRootName (string const&);
	string getTreeName (string const&);
	string collectStartingTreePhylip (string&, bool &);
	
	// allow easy access to data
	friend class BEASTXML;
	
	SimData (string const& fileName);
	~SimData () {};
};

#endif /* _SIMDATA_H_ */
