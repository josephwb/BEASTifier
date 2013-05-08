#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>

using namespace std;

#include "General.h";
#include "Analysis_Settings.h";


AnalysisSettings::AnalysisSettings ()
:
	treePrior("bd"), manipulateTreeTopology(true), logPhylograms(false), overwrite(false),
	mcmcLength(20000000), screenSampling(500), parameterSampling(1000), treeSampling(5000)
{
	intializeDefaults();
}

vector <string> AnalysisSettings::readListFromFile(string const& fileName) {
	vector <string> listElements;
	ifstream inputNames;
	inputNames.open(fileName.c_str());
	
	string line;
	
	while (getline(inputNames,line)) {
		if (!checkWhiteSpaceOnly(line)) {
			listElements.push_back(line);
		}
	}
	return (listElements);
}

void AnalysisSettings::setMcmcLength (string val) {
	mcmcLength = convertStringtoInt(val);
}

void AnalysisSettings::setScreenSampling (string val) {
	screenSampling = convertStringtoInt(val);
}

void AnalysisSettings::setParameterSampling (string val) {
	parameterSampling = convertStringtoInt(val);
}

void AnalysisSettings::setTreeSampling (string val) {
	treeSampling = convertStringtoInt(val);
}

void AnalysisSettings::setTreeManipulationFalse () {
	manipulateTreeTopology = false;
}

void AnalysisSettings::setLogPhylogramsTrue () {
	logPhylograms = true;
}

// default settings for substitution models, clock flavour, and tree priors
void AnalysisSettings::intializeDefaults () {
	string init[] = {"JC", "HKY", "GTR", "JC+G", "HKY+G", "GTR+G"}; // initialized set; may be superseded by user input
	vector <string> defaultModels(init, init + sizeof(init) / sizeof(string));
	models = defaultModels;
	
	string defClock = "ucln";
	string defTreePr = "bd";
	vector <string> defaultClock;
	defaultClock.push_back(defClock);
	vector <string> defaultTreePr;
	defaultTreePr.push_back(defTreePr);
	
	clockFlavours = defaultClock;
	treePriors = defaultTreePr;
}

void AnalysisSettings::setSubModels (vector <string> const& subModels) {
	models = subModels;
}

int AnalysisSettings::getNumSubModels () {
	return models.size();
}

string AnalysisSettings::getSubModel (int const& modelIndex) {
	return models[modelIndex];
}

void AnalysisSettings::setClockFlavours (vector <string> const& clockVals) {
	clockFlavours = clockVals;
	for (int i = 0; i < (int)clockFlavours.size(); i++) {
		checkClockFlavour(clockFlavours[i]);
	}
}

bool AnalysisSettings::checkClockFlavour (string const& clockString) {
	bool cool = true;
	if (clockString != "strict" && clockString != "ucln" && clockString != "uced" && clockString != "randlocal") {
		cool = false;
		ofstream errorReport("Error.BEASTifier.txt");
		errorReport << "BEASTifier  failed." << endl << "Error: clock flavour '";
		errorReport << clockString << "' not recognized." << endl;
		errorReport.close();
		cerr << endl << "BEASTifier failed." << endl << "Error: clock flavour '";
		cerr << clockString << "' not recognized. You fucked up, yo. Exiting." << endl << endl;
		exit(1);
	}
	return cool;
}

string AnalysisSettings::getClockFlavour (int const& clockIndex) {
	return clockFlavours[clockIndex];
}

int AnalysisSettings::getNumClockFlavours () {
	return clockFlavours.size();
}

void AnalysisSettings::setTreePriors (vector <string> const& treePriorVals) {
	treePriors = treePriorVals;
	for (int i = 0; i < (int)treePriors.size(); i++) {
		checkTreePrior(treePriors[i]);
	}
}

bool AnalysisSettings::checkTreePrior (string const& treePriorString) {
	bool cool = true;
	if (treePriorString != "bd" && treePriorString != "yule" && treePriorString != "concoal"
		&& treePriorString != "expcoal" && treePriorString != "logcoal") {
		cool = false;
		ofstream errorReport("Error.BEASTifier.txt");
		errorReport << "BEASTifier  failed." << endl << "Error: tree prior flavour '";
		errorReport << treePriorString << "' not recognized." << endl;
		errorReport.close();
		cerr << endl << "BEASTifier failed." << endl << "Error: tree prior flavour '";
		cerr << treePriorString << "' not recognized. You fucked up, yo. Exiting." << endl << endl;
		exit(1);
	}
	return cool;
}

string AnalysisSettings::getTreePrior (int const& treePriorIndex) {
	return treePriors[treePriorIndex];
}

int AnalysisSettings::getNumTreePriors () {
	return treePriors.size();
}

void AnalysisSettings::setRootPrior (vector <string> const& rootPriorVals) {
	rootPrior = rootPriorVals;
	if (rootPrior.size() != 3) {
		ofstream errorReport("Error.BEASTifier.txt");
		errorReport << "BEASTifier  failed." << endl;
		errorReport << "Error: Root prior must have three elements:" << endl;
		errorReport << "   '-rprior unif min_value max_value'" << endl;
		errorReport << "or" << endl;
		errorReport << "   '-rprior norm mean_value stdev_value'" << endl;
		errorReport.close();
		cerr << "BEASTifier  failed." << endl;
		cerr << "Error: Root prior must have three elements:" << endl;
		cerr << "   '-rprior unif min_value max_value'" << endl;
		cerr << "or" << endl;
		cerr << "   '-rprior norm mean_value stdev_value'" << endl << endl;
		exit(1);
	}
	
	checkPriorFlavour(rootPrior[0]); // flavour; unif or norm
	checkValidFloat(rootPrior[1]); // min (unif) or mean (norm)
	checkValidFloat(rootPrior[2]); // max (unif) or stdev (norm)
}

bool AnalysisSettings::checkPriorFlavour (string & priorString) {
	bool cool = true;
	if (priorString != "unif" && priorString != "norm") {
		cool = false;
		ofstream errorReport("Error.BEASTifier.txt");
		errorReport << "BEASTifier  failed." << endl << "Error: prior flavour '";
		errorReport << priorString << "' not recognized." << endl;
		errorReport.close();
		cerr << endl << "BEASTifier failed." << endl << "Error: prior flavour '";
		cerr << priorString << "' not recognized. You fucked up, yo. Exiting." << endl << endl;
		exit(1);
	}
	return cool;
}

void AnalysisSettings::setOverwriteTrue () {
	overwrite = true;
}
