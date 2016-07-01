#ifndef _ANALYSIS_SETTINGS_H_
#define _ANALYSIS_SETTINGS_H_

class AnalysisSettings {
    
    string treePrior;
    bool manipulateTreeTopology, logPhylograms, overwrite;
    int mcmcLength, screenSampling, parameterSampling, treeSampling;
    vector <string> rootPrior, models, clockFlavours, treePriors;
        
public:
    
    vector <string> readListFromFile(string const& fileName);
    void intializeDefaults ();
    
// mcmc parameters
    void setMcmcLength (string val);
    void setScreenSampling (string val);
    void setParameterSampling (string val);
    void setTreeSampling (string val);
    
    void setTreeManipulationFalse ();
    void setLogPhylogramsTrue ();
    
    void setRootPrior (vector <string> const& rootPriorVals);
    bool checkPriorFlavour (string & priorString);
    
    void setOverwriteTrue ();
    
// values to loop over
    void setSubModels (vector <string> const& subModels);
    int getNumSubModels ();
    string getSubModel (int const& modelIndex);
    
    void setClockFlavours (vector <string> const&  clockVals);
    bool checkClockFlavour (string const& clockString);
    string getClockFlavour (int const& clockIndex);
    int getNumClockFlavours ();
    
    
    void setTreePriors (vector <string> const& treePriorVals);
//    void setTreePrior (string val);
    bool checkTreePrior (string const& treePriorString);
    //string getTreePrior ();
    string getTreePrior (int const& treePriorIndex);
    int getNumTreePriors ();
    
    void readConfigFile (string const& fileName);
    
    // allow easy access to settings
    friend class BEASTXML;
    
    AnalysisSettings ();
    ~AnalysisSettings () {};
};

#endif /* _ANALYSIS_SETTINGS_H_ */
