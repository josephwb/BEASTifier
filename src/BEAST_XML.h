#ifndef _BEAST_XML_H_
#define _BEAST_XML_H_

class BEASTXML {
    
    ofstream BEAST_xml_code;
    string root, analyzeModel, clockFlavour, treePrior, XMLOutFileName, starterTree,
        partitionSubstitutionModel, partitionSiteModel;
    int numTaxa, numChar, mcmcLength, screenSampling, parameterSampling, treeSampling;
    bool manipulateTreeTopology, overwrite, logPhylograms, starterTreePresent;
    vector <string> rootPrior;
    vector <string> partitionSubstitutionModels;
    vector <string> partitionSiteModels;
    vector < vector <string> > taxaAlignment;
    
public:
    
    string setXMLOutFileName (bool const& overwrite);
    string getXMLOutFileName ();
    void setDNASubModel (string const& analyzeModel, string & partitionSubstitutionModel,
        string & partitionSiteModel);
    void setDNAModel (string const& analyzeModel, vector <string> & partitionSubstitutionModels,
        vector <string> & partitionSiteModels);
    
    void writeFile ();
    
    // Writing functions
    void writeXMLHeader ();
    void writeXMLTail ();
    void writeTaxonList (ofstream & BEAST_xml_code, int const& numTaxa,
        vector < vector <string> > const& taxaAlignment);
    void writeAlignment (ofstream & BEAST_xml_code, int const& numTaxa, int const& numChar,
        vector < vector <string> > const& taxaAlignment);
    void writePartitionInformation (ofstream & BEAST_xml_code);
    void writeTreePrior (ofstream & BEAST_xml_code, string const& treePrior,
        bool const& starterTreePresent, string const& starterTree);
    void writeTreeModel (ofstream & BEAST_xml_code, string const& treePrior);
    void writeClockModel (ofstream & BEAST_xml_code, int const& numTaxa, string const& clockFlavour);
    
    void writeSubstitutionModel (ofstream & BEAST_xml_code, string const& partitionSubstitutionModel);
    
//    void writeSubstitutionModels (ofstream & BEAST_xml_code, vector <string> const& partitionSubstitutionModels);
    
    void writeSiteModel (ofstream & BEAST_xml_code, string const& partitionSubstitutionModel,
        string const& partitionSiteModel);

//     void writeSiteModels (ofstream & BEAST_xml_code, vector <string> const& partitionSubstitutionModels,
//         vector <string> const& partitionSiteModels);
    
    void writeTreeLikelihoods (ofstream & BEAST_xml_code, string const& clockFlavour);
    
    void writeOperators (ofstream & BEAST_xml_code, string const& treePrior, bool const& manipulateTreeTopology,
        string const& partitionSubstitutionModel, string const& partitionSiteModel,
        string const& clockFlavour, int const& numTaxa);
    
//     void writeOperators (ofstream & BEAST_xml_code, string const& treePrior,
//         bool const& manipulateTreeTopology, vector <string> const& partitionSubstitutionModels,
//         vector <string> const& partitionSiteModels, string const& clockFlavour, int const& numTaxa);
    
    void writeMCMCParameters (ofstream & BEAST_xml_code, int const& mcmcLength, string const& clockFlavour, 
        string const& partitionSubstitutionModel, vector <string> const& rootPrior, string const& treePrior);
    
//     void writeMCMCParameters (ofstream & BEAST_xml_code, int const& mcmcLength, string const& clockFlavour, 
//         vector <string> const& partitionSubstitutionModels, vector <string> const& rootPrior, string const& treePrior);
    
    void writeScreenLog (ofstream & BEAST_xml_code, int const& screenSampling, string const& clockFlavour);
    
    void writeParameterLog (ofstream & BEAST_xml_code, int const& parameterSampling,
        string const& treePrior, string & clockFlavour, string const& partitionSubstitutionModel,
        string const& partitionSiteModel);
    
//     void writeParameterLog (ofstream & BEAST_xml_code, int const& parameterSampling,
//         string const& treePrior, string & clockFlavour, vector <string> const& partitionSubstitutionModels,
//         vector <string> const& partitionSiteModels);
    
    void writeTreeLogs (ofstream & BEAST_xml_code, string const& clockFlavour,
        bool const& logPhylograms, int const& treeSampling);

    BEASTXML (SimData & data, int const& modelIndex, int const& clockIndex, int const& treePriorIndex, AnalysisSettings ASet);
    ~BEASTXML () {};
};

#endif /* _BEAST_XML_H_ */
