/*

************************************************
            BEASTifier version 0.3
                Joseph W. Brown
              University of Idaho
              josephwb@umich.edu
                 January, 2013
************************************************

Program description: Generate a heap of BEAST xml file(s), looping over 1) alignment files,
and 2) substitution models.

To compile, type the following in a unix prompt:
g++ -Wall -m64 -O3 BEASTifier.cpp -o BEASTifier

Usage:
./BEASTifier -alist listOfAlignmenteFileNames [-mods listOfModels] [-clock clockFlavour] [-mcmc mcmcLength] [-ts treeSamplingInterval]
[-ps parameterSamplingInterval] [-ss screenSamplingInterval] [-logphy logPhylograms] [-tprior] [-fixtree] [-rprior args]

where:
	-alist: name of file listing alignment file names. Required; all other arguments are optional.
	-mods: name of file listing substitution model(s) to analyze data; default = JC,HKY,GTR,JC+G,HKY+G,GTR+G.
	   - list one model per line.
	   - supported models: K80, HKY, TrNef, TrN, K3P, K3Puf, TIMef, TIM, TVMef, TVM, SYM, GTR.
	   - model names must contain no spaces and at most one '+'.
	   - 'GTR+IG' = good; 'GTR+I+G' no es bueno.
	-clock: flavour of clock model. Options = 'strict' or 'ucln' or 'uced'; default = 'ucln'
	-mcmc: the number of mcmc generations to run analysis; default = 20000000 generations.
	-ts: the interval for sampling trees; default = every 5000 generations.
	-ps: the interval for sampling parameter values; default = every 1000 generations.
	-ss: the interval for printing results to standard output; default = every 500 generations.
	-logphy: turn on logging of phylograms (in addition to chronograms); default = don't.
	-tprior: specify tree prior. Options = 'bd', 'yule', 'concoal', 'expcoal', 'logcoal'; default = 'bd'.
	-fixtree: turn off topology manipulation operators (i.e. fix to input topology); default = estimate topology.
	-rprior: specify prior for root age. One of 'unif'' or 'norm'.
	   - if 'unif', expecting '-rprior unif min_value max_value'.
	   - if 'norm', expecting '-rprior norm mean_value stdev_value'.

*/

#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>

#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_num_procs() 1
#endif

using namespace std;

#include "General.h"
#include "Analysis_Settings.h"
#include "User_Interface.h"
#include "SimData.h"
#include "BEAST_XML.h"

// version information
double version = 0.37;
string month = "January";
int year = 2013;

bool DEBUG = false;

int main (int argc, char *argv[]) {

	vector <string> listFileNames;
 	int fileCounter = 0;
	
// default options are now in a AnalysisSettings object.
	AnalysisSettings ASet;
	
	printProgramInfo();
	processCommandLineArguments(argc, argv, listFileNames, ASet);
	
// 	cout << "listFileNames.size() = " << listFileNames.size() << endl;
// 	cout << "ASet.getSubModels().size() = " << ASet.getSubModels().size() << endl;
//	cout << "ASet.getTreeManipulation() = " << ASet.getTreeManipulation() << endl;
	
	for (int i = 0; i < int(listFileNames.size()); i++) { // loop over file names
	
// File-specific parameters are now stored in SimData object
		SimData Data(listFileNames[i]);
		
 		cout << endl << "Processing alignment '" << listFileNames[i] << "'..." << endl;
		
		for (int j = 0; j < int(ASet.getSubModels().size()); j++) { // loop over model flavours
			
			for (int k = 0; k < int(ASet.getClockFlavours().size()); k++) { // loop over clock flavours
			
				BEASTXML BXML(Data, j, k, ASet); // add looping over clock flavours as well ***
			
				cout << "	- creating BEAST file using substitution model '" << ASet.getSubModel(j)
					<< "' and clock flavour '" << ASet.getClockFlavour(k) << "'." << endl;
			
				fileCounter++;
				if (DEBUG) {cout << "Successfully created file '" << BXML.getXMLOutFileName() << "'." << endl;}
 			}
		}
	}
	
	cout << endl << endl << "Successfully created " << fileCounter << " BEAST input files. Hazzah!" << endl;
	cout << endl << "Fin." << endl;
	return 0;
}


// Not used
string collectStartingTreeNexus (string fileName, bool & starterTreePresent) {
// For starters, looking for simple pattern i.e. 'tree tree_name = [&R] tree_description;'
// Allow user to choose between multiple trees, if present
// Later - add ability to remove commented-out sections of tree e.g. rate estimates from a previous analysis

	string treeString;
	bool commentLine = false;
	bool whiteSpaceOnly = false;
	
	cout << endl << endl
	<< "*******************************" << endl
	<< "*** INITIALIZING CHRONOGRAM ***" << endl
	<< "*******************************" << endl << endl;
	
	checkValidInputFile(fileName);
	ifstream treeInput;
	treeInput.open(fileName.c_str());
	vector <string> tempTreeVector;
	vector <string> tempTreeNameVector;
	int treeCounter = 0;
	string line;
	
// Read in every non-empty (or non-whitespace), non-commented-out line
	while (getline(treeInput,line)) {
		int stringPosition = 0;
		commentLine = checkCommentLineNexus(line);
		whiteSpaceOnly = checkWhiteSpaceOnly(line);
		if (line.empty() || commentLine || whiteSpaceOnly) {
			continue;
		} else {
			if (checkStringValue(line, "tree", stringPosition)) {
				stringPosition++;
				string treeName = parseString(line, stringPosition);
				bool treeRead = false;
				while (!treeRead) {
					stringPosition++;
					if (checkStringValue(line, "=", stringPosition)) {
						stringPosition++;
					}
// Test for possible rooting option e.g. '[&R]'
					if (checkCommentLineNexus(parseString(line,stringPosition))) {
						stringPosition++;
					}
					tempTreeVector.push_back(parseString(line, stringPosition));
					tempTreeNameVector.push_back(treeName);
					treeCounter++;
					line.clear();
					treeRead = true;
				}
			} else {
				line.clear();
			}
		}
	}

	int userSelection = 0;
	bool validChoice = false;
	
	while (!validChoice) {
		cout << "Available trees:" << endl << endl;
		cout << " (1) Coalescent starting tree" << endl;

		for (int treeIter = 0; treeIter < treeCounter; treeIter++) {
			if (treeIter < 9) {
				cout << " ";
			}
			cout << "(" << treeIter + 2 << ") Input tree name: " << tempTreeNameVector[treeIter] << endl;
		}
		cout << endl << "Please select tree from list above: ";
		cin >> userSelection;
		
		if (cin.fail() || userSelection < 1 || userSelection > treeCounter + 1) {
			cout << endl << "*** Invalid input. Integer must be between 0 and " << treeCounter << " ***" << endl << endl;
			cin.clear(); 
			cin.ignore(200, '\n');
			continue;
		} else if (userSelection == 1) {
			starterTreePresent = false;
			validChoice = true;
			cout << endl << "Initializing chronogram will be generated under the coalescent process." << endl;
		} else {
			treeString = tempTreeVector[userSelection - 2];
			cout << endl << "Initializing chronogram set to tree '" << tempTreeNameVector[userSelection - 2] << "'." << endl;
			starterTreePresent = true;
			validChoice = true;
		}
	}
	
	cin.ignore(200, '\n');
	
	if (starterTreePresent) {
		cout << endl 
		<< "*** PLEASE NOTE ***" << endl
		<< "NO error-checking is currently performed to determine if this starting tree is consistent with monophyly/temporal constraints." << endl;
	}
	
	return treeString;
}