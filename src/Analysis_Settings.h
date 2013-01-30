#ifndef _ANALYSIS_SETTINGS_H_
#define _ANALYSIS_SETTINGS_H_

class AnalysisSettings {
	
	string treePrior, clockFlavour;
	bool manipulateTreeTopology, logPhylograms, overwrite;
	int mcmcLength, screenSampling, parameterSampling, treeSampling;
	vector <string> rootPrior, models, clockFlavours; // allow looping over clock flavours ***
		
public:
	
	vector <string> readListFromFile(string const& fileName);
	void initializeSubModels ();
	
	void setSubModels (vector <string> const& subModels);
	//void setSubModelsOLD (string & filename);
	vector <string> getSubModels ();
	string getSubModel (int const& modelIndex);
	
	
	void setMcmcLength (string val);
	int getMcmcLength ();
	void setScreenSampling (string val);
	int getScreenSampling ();
	void setParameterSampling (string val);
	int getParameterSampling ();
	void setTreeSampling (string val);
	int getTreeSampling ();
	
	void setClockFlavours (vector <string> const&  clockVals);
	void setClockFlavour (string val);
	bool checkClockFlavour (string const& clockString);
	string getClockFlavour (int const& clockIndex);
	vector <string> getClockFlavours ();
	
	void setTreeManipulationFalse ();
	bool getTreeManipulation ();
	void setLogPhylogramsTrue ();
	//bool getLogPhylograms ();
	
	void setTreePrior (string val);
	bool checkTreePrior (string const& treePriorString);
	//string getTreePrior ();
	
	void setRootPrior (vector <string> const& rootPriorVals);
	//void setRootPriorOLD (string val1, string val2, string val3);
	
	
	bool checkPriorFlavour (string & priorString);
	//vector <string> getRootPrior ();
	void setOverwriteTrue ();
	//bool getOverwrite ();
	
	void readConfigFile (string const& fileName);
	
	// allow easy access to settings
	friend class BEASTXML;
	
	AnalysisSettings ();
	~AnalysisSettings () {};
};

#endif /* _ANALYSIS_SETTINGS_H_ */