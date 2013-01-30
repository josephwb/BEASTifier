#ifndef _USER_INTERFACE_H_
#define _USER_INTERFACE_H_

void printProgramInfo();
void processCommandLineArguments (int argc, char *argv[], vector <string> & listFileNames,
	AnalysisSettings & ASet);
void usage ();
//vector <string> readFileList (string const& fileName);
// void printMenu (string const& prefix, vector <string> const& elementNames,
// 	vector <string> const& currentParameterValues, int const& numElements);

void readConfigFile (string const& fileName, AnalysisSettings & ASet,
	vector <string> & listFileNames);
bool checkComment (string const& val);
// void processCommandLineArgumentsOLD(int argc, char *argv[], vector <string> & listFileNames,
// 	vector <string> & models, string & treePrior, string & clockFlavour,
// 	bool & manipulateTreeTopology, bool & logPhylograms, int & mcmcLength, int & screenSampling,
// 	int & parameterSampling, int & treeSampling, vector <string> & rootPrior, bool & overwrite,
// 	AnalysisSettings & ASet);
// void getNumTaxaChar (string fileName, int & numTaxa, int & numChar, bool & interleavedData);
// string collectStartingTreePhylip (string fileName, bool & starterTreePresent);
// vector < vector <string> > collectTaxaAlignment (string fileName, int const& numTaxa,
// 	int const& numChar, bool const& interleavedData);
// string getTreeName(string const& rootName);
// bool checkClockFlavour (string & clockString);
// bool checkTreePrior (string & treePriorString);
// bool checkPriorFlavour (string & priorString);
// void parseModel (string & analyzeModel, vector <string> & partitionSubstitutionModels,
// 	vector <string> & partitionSiteModels);

#endif /* _USER_INTERFACE_H_ */