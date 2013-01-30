#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <iomanip>

using namespace std;

#include "General.h"
#include "Analysis_Settings.h";
#include "User_Interface.h"

extern bool DEBUG;
extern double version;
extern string month;
extern int year;

void printProgramInfo () {
	cout << endl << 
	"************************************************" << endl <<
	"             BEASTifier version " << version <<      endl <<
	"                Joseph W. Brown"                  << endl <<
	"             University of Michigan"              << endl <<
	"  Department of Ecology & Evolutionary Biology"   << endl <<
	"        Complaints: josephwb@umich.edu"           << endl <<
	"                 " << month <<", " << year <<        endl << 
	"************************************************" << endl << endl;
}

// *** add option to process config file instead - DONE
// need to error-check: parameters (maybe conflicting) in config vs. commandline - low priority
void processCommandLineArguments (int argc, char *argv[], vector <string> & listFileNames,
	AnalysisSettings & ASet)
{
	if (argc == 1) {
		cout << "No arguments given." << endl << endl;
		usage();
		exit(0);
	} else {
		for (int i = 1; i < argc; i++) {
			string temp = argv[i];
			
			if (temp == "-h" || temp == "-help") { // rewrite this; low priority - DONE
				cout
				<< "Program description: Generate a heap of BEAST xml files, looping over:" << endl
				<< " 1) alignment files, 2) substitution models, and 3) clock flavours." << endl
				<< endl
				<< "To compile, type the following in a unix prompt:" << endl << endl
				<< "   make" << endl << endl;
				usage();
				exit(0);
			} else if (temp == "-config") {
				i++;
				string temp = argv[i];
				readConfigFile (temp, ASet, listFileNames);
				continue;
			} else {
				cout
				<< "*** Unknown command-line argument '" << argv[i] << "' encountered. ***" << endl << endl;
				usage();
				cout << endl << "Exiting because of improper command-line argument '" << argv[i] << "' encountered." << endl << endl;
				exit(0);
			}
		}
	}
}

void usage() {
	cout
	<< "Usage:" << endl
	<< endl
	<< "BEASTifier utilizes a configuration file for all analysis parameters. Call as:" << endl
	<< endl
	<< "   ./BEASTifier -config config_filename" << endl
	<< endl
	<< "where 'config_filename' contains all analysis settings." << endl
	<< "Parameters are listed one per line, in any order. The character '#' is used for comments." << endl
	<< endl
	<< "Arguments:" << endl
	<< endl
	<< "   -alist: filename" << endl
	<< "      - name of text file listing alignment filenames." << endl
	<< "      - one alignment filename per line." << endl
	<< "      - Required; all other arguments are optional." << endl
	<< "   -mods: list substitution model(s) to analyze data" << endl
	<< "      - supported models: JC, K80, HKY, TrNef, TrN, K3P, K3Puf, TIMef, TIM, TVMef, TVM, SYM, GTR." << endl
	<< "      - if more than one model, separate by spaces." << endl
	<< "      - models themselves must contain no spaces and at most one '+'." << endl
	<< "         - e.g. 'GTR+IG' = good; 'GTR+I+G' no es bueno." << endl
	<< "      - default: -mods JC HKY GTR JC+G HKY+G GTR+G" << endl
	<< "   -clock: list of flavour(s) of clock model to implement." << endl
	<< "      - supported models: 'strict' or 'ucln' of 'uced'" << endl
	<< "      - default = -clock ucln" << endl
	<< "   -mcmc: the number of mcmc generations to run analysis." << endl
	<< "      - default: -mcmc 20000000" << endl
	<< "   -tsamp: the interval (in generations) for sampling trees." << endl
	<< "      - default: -ts 5000" << endl
	<< "   -psamp: the interval (in generations) for sampling parameter values." << endl
	<< "      - default: -ps 1000" << endl
	<< "   -ssamp: the interval (in generations) for printing results to standard output." << endl
	<< "      - default: -ss 500" << endl
	<< "   -logphy: turn on logging of phylograms (in addition to chronograms)." << endl
	<< "      -  default: don't log." << endl
	<< "   -tprior: specify tree prior." << endl
	<< "      - supported: 'bd', 'yule', 'concoal', 'expcoal', 'logcoal'." << endl
	<< "      - default = -tprior bd" << endl
	<< "   -fixtree: turn off topology manipulation operators (i.e. fix to input topology)." << endl
	<< "      - default = estimate topology." << endl
	<< "   -rprior: specify prior for root age." << endl
	<< "      - supported: 'unif'' or 'norm'." << endl
	<< "      - if 'unif', expecting '-rprior unif min_value max_value'." << endl
	<< "      - if 'norm', expecting -rprior norm mean_value stdev_value'." << endl
	<< "   -overwrite: overwrite existing files." << endl
	<< "      - default = don't overwrite; warn instead." << endl
	<< endl
	<< "Consult 'config.example' as a, well, example." << endl << endl;
}

// no longer used
/*
void printMenu (string const& prefix, vector <string> const& elementNames,
	vector <string> const& currentParameterValues, int const& numElements)
{
// Format assumed:
// '('index_value')' element_type(prefix) element_name ':' current_value
	
	string maxLengthName;
	string maxLengthInteger;
	string maxLengthParameter;
	
// Get longest integer
	maxLengthInteger = convertIntToString(numElements);
// Get longest individual element name and parameter value - used for formatting
	maxLengthName = getLongestName(elementNames);
	maxLengthParameter = getLongestName(currentParameterValues);
	
// Print!
	for (int elementIter = 0; elementIter < numElements; elementIter++) {
// Print out leading spaces for index
		string tempIndexString = convertIntToString(elementIter + 1);
		printFormattingSpaces(maxLengthInteger, tempIndexString);
		cout << "(" << tempIndexString <<") " << prefix << " " << elementNames[elementIter];

// Print out trailing spaces for element name
		string tempNameString = elementNames[elementIter];
		printFormattingSpaces(maxLengthName, tempNameString);
		cout << ": ";
		
// Print out leading spaces for parameter value
		string tempParameterString = currentParameterValues[elementIter];
		printFormattingSpaces(maxLengthParameter, tempParameterString);
		cout << currentParameterValues[elementIter] << endl;
	}
}
*/

void readConfigFile (string const& fileName, AnalysisSettings & ASet,
	vector <string> & listFileNames)
{
	ifstream configInput;
	configInput.open(fileName.c_str());
	
	vector <string> tempVect;
	string line;
	
	while (getline(configInput, line)) {
		if (!checkWhiteSpaceOnly(line)) {
			tempVect = tokenizeString(line);
			if (!checkComment(tempVect[0])) {
				if (tempVect[0] == "-alist") {
					string fileName = tempVect[1];
					checkValidInputFile(fileName);
					listFileNames = readFileList(fileName);
					continue;
				} else if (tempVect[0] == "-mods") {
					tempVect.erase(tempVect.begin());
					ASet.setSubModels(tempVect);
					continue;
				} else if (tempVect[0] == "-mcmc") {
					ASet.setMcmcLength(tempVect[1]);
					continue;
				} else if (tempVect[0] == "-tsamp") {
					ASet.setTreeSampling(tempVect[1]);
					continue;
				} else if (tempVect[0] == "-psamp") {
					ASet.setParameterSampling(tempVect[1]);
					continue;
				} else if (tempVect[0] == "-ssamp") {
					ASet.setScreenSampling(tempVect[1]);
					continue;
				} else if (tempVect[0] == "-logphy") {
					ASet.setLogPhylogramsTrue();
					continue;
				} else if (tempVect[0] == "-tprior") {
					ASet.setTreePrior(tempVect[1]);
					continue;
				} else if (tempVect[0] == "-fixtree") {
					ASet.setTreeManipulationFalse();
					continue;
				} else if (tempVect[0] == "-clock") {
					tempVect.erase(tempVect.begin());
					ASet.setClockFlavours(tempVect);
					continue;
				} else if (tempVect[0] == "-rprior") {
					tempVect.erase(tempVect.begin());
					ASet.setRootPrior(tempVect);
					continue;
				} else if (tempVect[0] == "-overwrite") {
					ASet.setOverwriteTrue();
					continue;
				} else {
					cout << endl
					<< "*** Unknown configuration file argument '" << tempVect[0] << "' encountered. ***" << endl << endl;
					usage();
					cout << endl << "Exiting because of improper configuration file argument '" << tempVect[0] << "' encountered." << endl << endl;
					exit(0);
				}
			}
			tempVect.clear();
		}
	}
}

bool checkComment (string const& val) {
	bool comment = false;
	if (val[0] == '#') {
		comment = true;
	}
	return comment;
}

/*
void getNumTaxaChar (string fileName, int & numTaxa, int & numChar, bool & interleavedData) {
	ifstream inputUserFile;
	bool commentLine = false;
	bool whiteSpaceOnly = false;
	bool numTaxaEncountered = false;
	bool numCharEncountered = false;
	bool equalSignEncountered = false;
	bool semicolonEncountered = false;
	bool matrixEncountered = false; // know to stop looking
	
	inputUserFile.open(fileName.c_str());
	string line;
	
// Looking for pattern like 'dimensions ntax=53 nchar=16620;'
// 	- can be in either order, but must be stated on same line (for now)
// 	- no spaces allowed next to equal sign (for now)

// Looking for pattern like 'Format datatype=dna [gap=-] [missing=?] {[interleave=yes] or [interleave]};'
// 	- no spaces allowed next to equal sign (for now)

	while (getline(inputUserFile,line) && !matrixEncountered) {
		int stringPosition = 0;
		commentLine = checkCommentLineNexus(line);
		whiteSpaceOnly = checkWhiteSpaceOnly(line);
		if (line.empty() || commentLine || whiteSpaceOnly) {
			continue;
		} else {
			if (checkStringValue(line, "matrix", stringPosition)) {	// Done - won't find the information further down; really only used for 'interleave'
				matrixEncountered = true;
				continue;
			} else if (checkStringValue(line, "dimensions", stringPosition)) {
				stringPosition = 0;
				while (!numTaxaEncountered || !numCharEncountered) {
					stringPosition++;
					string tempString = removeStringSuffix(parseString(line,stringPosition), '=', equalSignEncountered);
					if (checkStringValue(tempString, "ntax", 0)) {
						tempString = removeStringPrefix(parseString(line,stringPosition), '=');
						tempString = removeStringSuffix(tempString, ';', numTaxaEncountered);
						numTaxa = convertStringtoInt(tempString);
						if (DEBUG) {cout << "NTax = " << numTaxa << endl;}
						numTaxaEncountered = true;
					}
					if (checkStringValue(tempString, "nchar", 0))
					{
						tempString = removeStringPrefix(parseString(line,stringPosition), '=');
						tempString = removeStringSuffix(tempString, ';', numCharEncountered);
						numChar = convertStringtoInt(tempString);
						if (DEBUG) {cout << "NChar = " << numChar << endl;}
						numCharEncountered = true;
					}
				}
			} else if (checkStringValue(line, "format", stringPosition)) {
				stringPosition = 0;
				while (!semicolonEncountered) {
					stringPosition++;
					
					string tempString = removeStringSuffix(parseString(line,stringPosition), '=', equalSignEncountered); // where '=' is used i.e. 'interleave=yes;'
					if (checkStringValue(tempString, "interleave", 0)) {
						tempString = removeStringPrefix(parseString(line,stringPosition), '=');
						tempString = removeStringSuffix(tempString, ';', semicolonEncountered);
						if (checkStringValue(tempString, "yes", 0)) {
							interleavedData = true;
							semicolonEncountered = true;
							cout << "Data are in interleaved format." << endl;
							continue;
						} else if (checkStringValue(tempString, "no", 0)) {
							interleavedData = false;
							semicolonEncountered = true;
							cout << "Data are not in interleaved format." << endl;
							continue;
						}
					}
					tempString = removeStringSuffix(parseString(line,stringPosition), ';', semicolonEncountered); // where '=' is NOT used i.e. 'interleave;'
					if (checkStringValue(tempString, "interleave", 0)) {										  // OR a space follows interleave declaration i.e. 'interleave [=];'
						interleavedData = true;
						semicolonEncountered = true;
						cout << "Data are in interleaved format." << endl;
						continue;
					}
				}
			}
		}
	}
	inputUserFile.close();

	cout << endl;
	inputUserFile.close();
}

string collectStartingTreePhylip (string fileName, bool & starterTreePresent) {
// Rudimentary. if a valid line is found in a file, assume it is a phylip tree.
	string treeString;
	bool commentLine = false;
	bool whiteSpaceOnly = false;
	
	checkValidInputFile(fileName);
	ifstream treeInput;
	treeInput.open(fileName.c_str());
	string line;
	
// Read in every non-empty (or non-whitespace), non-commented-out line
	while (getline(treeInput,line)) {
		commentLine = checkCommentLineNexus(line);
		whiteSpaceOnly = checkWhiteSpaceOnly(line);
		if (line.empty() || commentLine || whiteSpaceOnly) {
			continue;
		} else {
			treeString = line;
			starterTreePresent = true;
		}
	}
	return treeString;
}

vector < vector <string> > collectTaxaAlignment (string fileName, int const& numTaxa,
	int const& numChar, bool const& interleavedData)
{
// PLEASE NOTE: search strategy below uses very strict format assumptions - that which is exported by PAUP*
// 	- do not be surprised if this fucks up - it is probably a simple rearrangement of terms

	vector < vector <string> > taxaAlignment;
	ifstream inputAlignment;
	vector <string> tempStringVector;
	bool commentLine = false;
	bool whiteSpaceOnly = false;
	bool matrixEncountered = false;
	bool allCharacterRead = false;
	int numCharRead = 0;
	
	inputAlignment.open(fileName.c_str());
	string line;
	
	if (!interleavedData) {
// Ignore lines until 'matrix' is encountered
		while (!matrixEncountered) {
			getline(inputAlignment,line);
			commentLine = checkCommentLineNexus(line);
			whiteSpaceOnly = checkWhiteSpaceOnly(line);
			
			if (line.empty() || commentLine || whiteSpaceOnly) {
				continue;
			} else if (checkStringValue(line, "matrix", 0)) {
				matrixEncountered = true;
			} else {
				continue;
			}
		}
		numCharRead = 0;
// Read in every non-empty (or non-whitespace), non-commented-out line
		for (int taxonIter = 0; taxonIter < numTaxa; taxonIter++) {
			getline(inputAlignment,line);
			commentLine = checkCommentLineNexus(line);
			whiteSpaceOnly = checkWhiteSpaceOnly(line);
			
			if (line.empty() || commentLine || whiteSpaceOnly) {
				taxonIter--;
				continue;
			} else {
// First string is taxon name, second is sequence
				if (DEBUG) {cout << "Reading in taxon '" << parseString(line, 0) << "'..." << endl;}
				tempStringVector.push_back(parseString(line, 0));
				tempStringVector.push_back(parseString(line, 1));
				
				taxaAlignment.push_back(tempStringVector);
				
				tempStringVector.clear();
				
// Count sites encountered - for error-checking; only checking first sequence (for now)
				if (taxonIter == 0) {
					int charCounter = 0;
					for (string::size_type iterCharacters = 0; iterCharacters < parseString(line, 1).size(); iterCharacters++) {
						charCounter++;
					}
					numCharRead += charCounter;
				}
			}
		}
		if (numCharRead == numChar) {
			allCharacterRead = true;
			if (DEBUG) {cout << "numCharRead (" << numCharRead << ") == numChar (" << numChar << ") declared in Nexus file. Woo-hoo!" << endl;}
		}
	} else if (interleavedData) {
// Ignore lines until 'matrix' is encountered
		while (!matrixEncountered) {
			getline(inputAlignment,line);
			commentLine = checkCommentLineNexus(line);
			whiteSpaceOnly = checkWhiteSpaceOnly(line);
			
			if (line.empty() || commentLine || whiteSpaceOnly) {
				continue;
			} else if (checkStringValue(line, "matrix", 0)) {
				matrixEncountered = true;
			} else {
				continue;
			}
		}
		bool firstPass = true;
		numCharRead = 0;
		while (!allCharacterRead) {
// Read in every non-empty (or non-whitespace), non-commented-out line
			for (int taxonIter = 0; taxonIter < numTaxa; taxonIter++) {
				getline(inputAlignment,line);
				commentLine = checkCommentLineNexus(line);
				whiteSpaceOnly = checkWhiteSpaceOnly(line);
				
				if (line.empty() || commentLine || whiteSpaceOnly) {
					taxonIter--;
					continue;
				} else {
// First string is taxon name, second is sequence
					string taxonName = (parseString(line, 0));
					string taxonSequence = (parseString(line, 1));
					
					if (firstPass) {
						if (DEBUG) {cout << "Reading in taxon '" << parseString(line, 0) << "'..." << endl;}
						tempStringVector.push_back(taxonName);		// Taxon name
						tempStringVector.push_back(taxonSequence);	// sequence
						taxaAlignment.push_back(tempStringVector);
					}
					if (!firstPass) {
						taxaAlignment[taxonIter][1] += taxonSequence;
					}
					
// Count sites encountered - for error-checking; only checking first sequence (for now)
					if (taxonIter == 0) {
						int charCounter = 0;
						for (string::size_type iterCharacters = 0; iterCharacters < taxonSequence.size(); iterCharacters++) {
							charCounter++;
						}
						numCharRead += charCounter;
					}
					tempStringVector.clear();
				}
			}
			firstPass = false;
			if (numCharRead == numChar) {
				allCharacterRead = true;
				if (DEBUG) {cout << "numCharRead (" << numCharRead << ") == numChar (" << numChar << ") declared in Nexus file. Woo-hoo!" << endl;}
			}
		}
	}
	return taxaAlignment;
}

string getTreeName(string const& rootName) {
// form: b_1_d_0.5_a_0.9_n_100_sim_JC_rep_4 - 12 elements
	bool done = false;
	string treeName = rootName;
	string repNumber = getStringElement(rootName, '_', 12);
	for (int i = 12; i > 8; i--) {
		treeName = removeStringSuffix(treeName, '_', done);
	}
	treeName = treeName + "_rep_" + repNumber + ".phy";
	
	return treeName;
}

bool checkClockFlavour (string & clockString) {
	bool cool = true;
	if (clockString != "strict" && clockString != "ucln" && clockString != "uced") {
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

bool checkTreePrior (string & treePriorString) {
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

bool checkPriorFlavour (string & priorString) {
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

void parseModel (string & analyzeModel, vector <string> & partitionSubstitutionModels,
	vector <string> & partitionSiteModels)
{
	string subModel = getStringElement(analyzeModel, '+', 1);
	partitionSubstitutionModels.push_back(subModel);
	
	if (DEBUG) {cout << "subModel = " << subModel << endl;}
	
	string hetModel = removeStringPrefix(analyzeModel, '+');
	
	if (hetModel == subModel) {
		hetModel = "none";
	}
	partitionSiteModels.push_back(hetModel);
	if (DEBUG) {cout << "hetModel = " << hetModel << endl;}
}


void processCommandLineArgumentsOLD(int argc, char *argv[], vector <string> & listFileNames,
	vector <string> & models, string & treePrior, string & clockFlavour,
	bool & manipulateTreeTopology, bool & logPhylograms, int & mcmcLength, int & screenSampling,
	int & parameterSampling, int & treeSampling, vector <string> & rootPrior, bool & overwrite,
	AnalysisSettings & ASet)
{
	if (argc == 1) {
		cout << "No arguments given." << endl << endl;
		usage();
		exit(0);
	} else {
		for (int i = 1; i < argc; i++) {
			string temp = argv[i];
			
			if (temp == "-h" || temp == "-help") { // rewrite this; low priority
				cout << endl
				<< "Program description: Generate a heap of BEAST xml file(s), looping over:" << endl
				<< " 1) alignment files" << endl
				<< " 2) substitution models." << endl
				<< endl
				<< "To compile, type the following in a unix prompt:" << endl << endl
				<< "g++ -Wall -m64 -O3 BEASTifier.cpp -o BEASTifier" << endl << endl;
				usage();
				exit(0);
			} else if (temp == "-alist") {
				i++;
				string fileName = argv[i];
				checkValidInputFile(fileName);
				listFileNames = readFileList(fileName);
				continue;
			} else if (temp == "-mods") {
				i++;
				string fileName = argv[i];
				checkValidInputFile(fileName);
				models.clear();
				models = readFileList(fileName);
				ASet.setSubModels(fileName);
				continue;
			} else if (temp == "-mcmc") {
				i++;
				mcmcLength = convertStringtoInt(argv[i]);
				ASet.setMcmcLength(argv[i]);
				continue;
			} else if (temp == "-ts") {
				i++;
				treeSampling = convertStringtoInt(argv[i]);
				ASet.setTreeSampling(argv[i]);
				continue;
			} else if (temp == "-ps") {
				i++;
				parameterSampling = convertStringtoInt(argv[i]);
				ASet.setParameterSampling(argv[i]);
				continue;
			} else if (temp == "-ss") {
				i++;
				screenSampling = convertStringtoInt(argv[i]);
				ASet.setScreenSampling(argv[i]);
				continue;
			} else if (temp == "-logphy") {
				logPhylograms = true;
				ASet.setLogPhylogramsTrue();
				continue;
			} else if (temp == "-tprior") {
				i++;
				treePrior = argv[i];
				checkTreePrior(treePrior);
				ASet.setTreePrior(argv[i]);
				continue;
			} else if (temp == "-fixtree") {
				manipulateTreeTopology = false;
				ASet.setTreeManipulationFalse();
				continue;
			} else if (temp == "-clock") {
				i++;
				clockFlavour = argv[i];
				checkClockFlavour(clockFlavour);
				ASet.setClockFlavour(argv[i]);
				continue;
			} else if (temp == "-rprior") {
				i++;
				ASet.setRootPrior(argv[i], argv[i+1], argv[i+2]);
				string prior = argv[i];
				checkPriorFlavour(prior);
				rootPrior.push_back(prior);   // flavour; unif or norm
				i++;
				checkValidFloat(argv[i]);
				rootPrior.push_back(argv[i]); // min (unif) or mean (norm)
				i++;
				checkValidFloat(argv[i]);
				rootPrior.push_back(argv[i]); // max (unif) or stdev (norm)
			} else if (temp == "-overwrite") {
				overwrite = true;
				ASet.setOverwriteTrue();
				continue;
			} else {
				cout
				<< "Unknown command-line argument '" << argv[i] << "' encountered." << endl << endl;
				usage();
				exit(0);
			}
			cout << endl;
		}
	}
}
*/
