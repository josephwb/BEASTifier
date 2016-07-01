/*
************************************************
            BEASTifier version 0.41
                Joseph W. Brown
             University of Michigan
  Department of Ecology & Evolutionary Biology
        Complaints: josephwb@umich.edu
                   May, 2013
************************************************

Program description: Generate a heap of BEAST xml file(s), looping over 1) alignment files,
2) substitution models, 3) clock flavours, and 4) tree priors.

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
double version = 0.41;
string month = "May";
int year = 2013;

bool DEBUG = false;

int main (int argc, char *argv[]) {

    vector <string> listFileNames;
     int fileCounter = 0;
    
// default options are now in a AnalysisSettings object.
    AnalysisSettings ASet;
    
    printProgramInfo();
    processCommandLineArguments(argc, argv, listFileNames, ASet);
    
//     cout << "listFileNames.size() = " << listFileNames.size() << endl;
    
    for (int i = 0; i < int(listFileNames.size()); i++) { // loop over file names
    
// File-specific parameters are now stored in SimData object
        SimData Data(listFileNames[i]);
        
         cout << endl << "Processing alignment '" << listFileNames[i] << "'..." << endl;
        
        for (int j = 0; j < ASet.getNumSubModels(); j++) { // loop over model flavours
            for (int k = 0; k < ASet.getNumClockFlavours(); k++) { // loop over clock flavours
                for (int l = 0; l < ASet.getNumTreePriors(); l++) { // loop over tree priors
                    BEASTXML BXML(Data, j, k, l, ASet);
            
                    cout << "    - creating BEAST file using substitution model '" << ASet.getSubModel(j)
                        << "', clock flavour '" << ASet.getClockFlavour(k)
                        << "' and tree prior '" << ASet.getTreePrior(l) << "'." << endl;
            
                    fileCounter++;
                    if (DEBUG) {cout << "Successfully created file '" << BXML.getXMLOutFileName() << "'." << endl;}
                }
             }
        }
    }
    
    cout << endl << endl << "Successfully created " << fileCounter << " BEAST input files. Hazzah!" << endl;
    cout << endl << "Fin." << endl;
    return 0;
}
