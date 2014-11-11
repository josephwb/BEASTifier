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
	"            BEASTifier version " << version <<      endl <<
	"                Joseph W. Brown"                  << endl <<
	"             University of Michigan"              << endl <<
	"  Department of Ecology & Evolutionary Biology"   << endl <<
	"        Complaints: josephwb@umich.edu"           << endl <<
	"                   " << month <<", " << year <<      endl << // fix centering here
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
	<< "      - supported models: 'strict' or 'ucln' or 'uced'" << endl
	<< "      - default = -clock ucln" << endl
	<< "   -mcmc: the number of mcmc generations to run analysis." << endl
	<< "      - default: -mcmc 20000000" << endl
	<< "   -tsamp: the interval (in generations) for sampling trees." << endl
	<< "      - default: -tamps 5000" << endl
	<< "   -psamp: the interval (in generations) for sampling parameter values." << endl
	<< "      - default: -psamp 1000" << endl
	<< "   -ssamp: the interval (in generations) for printing results to standard output." << endl
	<< "      - default: -ssamp 500" << endl
	<< "   -logphy: turn on logging of phylograms (in addition to chronograms)." << endl
	<< "      -  default: don't log." << endl
	<< "   -tprior: list of tree prior(s)" << endl
	<< "      - supported: 'bd', 'yule', 'concoal', 'expcoal', 'logcoal'." << endl
	<< "      - default = -tprior bd" << endl
	<< "   -fixtree: turn off topology manipulation operators (i.e. fix to input topology)." << endl
	<< "      - default = estimate topology." << endl
	<< "   -rprior: specify prior for root age." << endl
	<< "      - supported: 'unif'' or 'norm'." << endl
	<< "      - if 'unif', expecting '-rprior unif min_value max_value'." << endl
	<< "      - if 'norm', expecting '-rprior norm mean_value stdev_value'." << endl
	<< "   -overwrite: overwrite existing files." << endl
	<< "      - default = don't overwrite; warn instead." << endl
	<< endl
	<< "Consult 'config.example' as a, well, example." << endl << endl;
}

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
					tempVect.erase(tempVect.begin());
					ASet.setTreePriors(tempVect);
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
