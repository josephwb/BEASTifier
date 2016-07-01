#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
#include <cstdlib>

using namespace std;

#include "General.h"
#include "Analysis_Settings.h"
#include "SimData.h";
#include "BEAST_XML.h"

extern bool DEBUG;
extern double version;
extern string month;
extern int year;

// functions for writing BEAST xml files

// functions marked as "not used at the moment" were taken from code for partitioned model analyses

BEASTXML::BEASTXML (SimData & data, int const& modelIndex, int const& clockIndex, int const& treePriorIndex, AnalysisSettings ASet)
{
// extract information from SimData object
    root = data.root;
    numTaxa = data.numTaxa;
    numChar = data.numChar;
    taxaAlignment = data.taxaAlignment;
    starterTree = data.starterTree;
    starterTreePresent = data.starterTreePresent;
    
// analysis parameters. others will include e.g. mcmc parameters
    analyzeModel = ASet.getSubModel(modelIndex);
    setDNASubModel(analyzeModel, partitionSubstitutionModel, partitionSiteModel);
    clockFlavour = ASet.getClockFlavour(clockIndex);
    treePrior = ASet.getTreePrior(treePriorIndex);
    
    rootPrior = ASet.rootPrior;
    mcmcLength = ASet.mcmcLength;
    screenSampling = ASet.screenSampling;
    parameterSampling = ASet.parameterSampling;
    treeSampling = ASet.treeSampling;
    logPhylograms = ASet.logPhylograms;
    overwrite = ASet.overwrite;
    manipulateTreeTopology = ASet.manipulateTreeTopology;
        
    XMLOutFileName = setXMLOutFileName(overwrite);
    BEAST_xml_code.open(XMLOutFileName.c_str());
    
    writeFile();
    
    BEAST_xml_code.close();
}

void BEASTXML::writeFile () {
    writeXMLHeader ();
    writeTaxonList (BEAST_xml_code, numTaxa, taxaAlignment);
    writeAlignment (BEAST_xml_code, numTaxa, numChar, taxaAlignment);
    writePartitionInformation (BEAST_xml_code);
    writeTreePrior (BEAST_xml_code, treePrior, starterTreePresent, starterTree);
    writeTreeModel (BEAST_xml_code, treePrior);
    writeClockModel (BEAST_xml_code, numTaxa, clockFlavour);
    
    writeSubstitutionModel (BEAST_xml_code, partitionSubstitutionModel);
    
//    writeSubstitutionModels (BEAST_xml_code, partitionSubstitutionModels);
    
    writeSiteModel (BEAST_xml_code, partitionSubstitutionModel,    partitionSiteModel);
    
//    writeSiteModels (BEAST_xml_code, partitionSubstitutionModels, partitionSiteModels);
    
    writeTreeLikelihoods (BEAST_xml_code, clockFlavour);
    
    writeOperators (BEAST_xml_code, treePrior, manipulateTreeTopology,
        partitionSubstitutionModel, partitionSiteModel, clockFlavour, numTaxa);
    
//     writeOperators (BEAST_xml_code, treePrior, manipulateTreeTopology,
//         partitionSubstitutionModels, partitionSiteModels, clockFlavour, numTaxa);
    
    writeMCMCParameters (BEAST_xml_code, mcmcLength, clockFlavour,
        partitionSubstitutionModel, rootPrior, treePrior);
    
//     writeMCMCParameters (BEAST_xml_code, mcmcLength, clockFlavour,
//         partitionSubstitutionModels, rootPrior, treePrior);
    
    writeScreenLog (BEAST_xml_code, screenSampling, clockFlavour);
    
    writeParameterLog (BEAST_xml_code, parameterSampling, treePrior, clockFlavour,
        partitionSubstitutionModel, partitionSiteModel);
    
//     writeParameterLog (BEAST_xml_code, parameterSampling, treePrior, clockFlavour,
//         partitionSubstitutionModels, partitionSiteModels);
    
    
    writeTreeLogs (BEAST_xml_code, clockFlavour, logPhylograms,
        treeSampling);
    writeXMLTail();
}

string BEASTXML::setXMLOutFileName (bool const& overwrite) {
    XMLOutFileName = root + "_analyze-" + analyzeModel;
    if (manipulateTreeTopology) {
        XMLOutFileName = XMLOutFileName + "_est-top";
    } else if (!manipulateTreeTopology) {
        XMLOutFileName = XMLOutFileName + "_fix-top";
    }
    
    XMLOutFileName = XMLOutFileName + "_" + treePrior + "-prior";
    
    if (!rootPrior.empty()) {
        XMLOutFileName = XMLOutFileName + "_" + rootPrior[0] + "-root";
    }
    
    XMLOutFileName = XMLOutFileName + '_' + clockFlavour + "-clock";
    XMLOutFileName = XMLOutFileName + ".xml";
    
// Check if file exists/is writable
    bool validFileName = false;
    while (!validFileName) {
        validFileName = checkValidOutputFile(XMLOutFileName, overwrite);
    }
    return XMLOutFileName;
}

string BEASTXML::getXMLOutFileName () {
    return XMLOutFileName;
}

void BEASTXML::setDNASubModel (string const& analyzeModel, string & partitionSubstitutionModel,
    string & partitionSiteModel)
{
    partitionSubstitutionModel = getStringElement(analyzeModel, '+', 1);
    
    if (DEBUG) {cout << "subModel = " << partitionSubstitutionModel << endl;}
    
    partitionSiteModel = removeStringPrefix(analyzeModel, '+');
    
    if (partitionSiteModel == partitionSubstitutionModel) {
        partitionSiteModel = "none";
    }
    
    if (DEBUG) {cout << "partitionSiteModel = " << partitionSiteModel << endl;}
}

void BEASTXML::setDNAModel (string const& analyzeModel, vector <string> & partitionSubstitutionModels,
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

void BEASTXML::writeXMLHeader () {
    BEAST_xml_code 
    << "<?xml version=\"1.0\" standalone=\"yes\"?>" << endl
    << endl
    << " <!-- Generated by BEASTifier version " << version << " -->" << endl
    << " <!--           Joseph W. Brown           -->" << endl
    << " <!--        University of Michigan       -->" << endl
    << " <!--          josephwb@umich.edu         -->" << endl
    << " <!--            " << month <<", " << year << "            -->" <<  endl
    << endl
    << "<beast>" << endl << endl;
}

void BEASTXML::writeXMLTail () {
    BEAST_xml_code << endl
    << "    <report>" << endl
    << "        <property name=\"timer\">" << endl
    << "            <object idref=\"mcmc\"/>" << endl
    << "        </property>" << endl
    << "    </report>" << endl
    << endl
    << "</beast>" << endl;
    
    BEAST_xml_code.close();
}

void BEASTXML::writeAlignment (ofstream & BEAST_xml_code, int const& numTaxa, int const& numChar,
    vector < vector <string> > const& taxaAlignment)
{
    BEAST_xml_code
    << "<!-- *** NUCLEOTIDE ALIGNMENT (refers to taxa above) *** -->" << endl
    << "    <!-- numTaxa = " << numTaxa << " numChar = " << numChar << " -->" << endl
    << "    <alignment id=\"alignment\" dataType=\"nucleotide\">" << endl;
    
    for (int taxonIter = 0; taxonIter < numTaxa; taxonIter++) {
        BEAST_xml_code
        << "        <sequence>" << endl
        << "            <taxon idref=\"" << taxaAlignment[taxonIter][0] << "\"/>" << endl
        << "            " << taxaAlignment[taxonIter][1] << endl
        << "        </sequence>" << endl;
    }
    BEAST_xml_code
    << "    </alignment>" << endl << endl;
}

void BEASTXML::writeTaxonList (ofstream & BEAST_xml_code, int const& numTaxa, vector < vector <string> > const& taxaAlignment) {
    BEAST_xml_code
    << "<!-- *** TAXON LIST *** -->" << endl
    << "    <!-- numTaxa = " << numTaxa << " -->" << endl
    << "    <taxa id=\"taxa\">" << endl;
    
    for (int taxonIter = 0; taxonIter < numTaxa; taxonIter++) {
        BEAST_xml_code
        << "        <taxon id=\"" << taxaAlignment[taxonIter][0] << "\"/>" << endl;
    }
    BEAST_xml_code
    << "    </taxa>" << endl << endl;
}

void BEASTXML::writeTreePrior (ofstream & BEAST_xml_code, string const& treePrior,
    bool const& starterTreePresent, string const& starterTree)
{
    if (treePrior == "bd") {
        BEAST_xml_code
        << "<!-- *** PRIOR ON NODE AGES *** -->" << endl
        << "    <!-- A prior on the distribution node heights defined given a birth-death speciation process -->" << endl
        << "    <birthDeathModel id=\"birthDeath\" units=\"substitutions\">" << endl
        << "        <birthMinusDeathRate>" << endl
        << "            <parameter id=\"birthDeath.BminusDRate\" value=\"1.0\" lower=\"0.0\" upper=\"Infinity\"/>" << endl
        << "        </birthMinusDeathRate>" << endl
        << "        <relativeDeathRate>" << endl
        << "            <parameter id=\"birthDeath.DoverB\" value=\"0.5\" lower=\"0.0\" upper=\"1.0\"/>" << endl
        << "        </relativeDeathRate>" << endl
        << "    </birthDeathModel>" << endl << endl << endl;
        
        if (!starterTreePresent) {
            BEAST_xml_code
            << "    <!-- This is a simple constant population size coalescent model  -->" << endl
            << "    <!-- that is used to generate an initial tree for the chain. -->" << endl
            << "    <constantSize id=\"initialDemo\" units=\"substitutions\">" << endl
            << "        <populationSize>" << endl
            << "            <parameter id=\"initialDemo.popSize\" value=\"100.0\"/>" << endl
            << "        </populationSize>" << endl
            << "    </constantSize>" << endl << endl
            << "    <!-- Generate a random starting tree under the coalescent process -->" << endl
            << "    <coalescentTree id=\"startingTree\" rootHeight=\"0.3\">" << endl
            << "        <taxa idref=\"taxa\"/>" << endl
            << "        <constantSize idref=\"initialDemo\"/>" << endl
            << "    </coalescentTree>" << endl << endl << endl;
        }
    } else if (treePrior == "yule") { // yule
        BEAST_xml_code
        << "<!-- *** PRIOR ON NODE AGES *** -->" << endl
        << "    <!-- A prior on the distribution node heights defined given a Yule speciation process (a pure birth process) -->" << endl
        << "    <yuleModel id=\"yule\" units=\"substitutions\">" << endl
        << "        <birthRate>" << endl
        << "            <parameter id=\"yule.birthRate\" value=\"1.0\" lower=\"0.0\" upper=\"Infinity\"/>" << endl
        << "        </birthRate>" << endl
        << "    </yuleModel>" << endl << endl << endl;
        
        if (!starterTreePresent) {
            BEAST_xml_code
            << "    <!-- This is a simple constant population size coalescent model  -->" << endl
            << "    <!-- that is used to generate an initial tree for the chain. -->" << endl
            << "    <constantSize id=\"initialDemo\" units=\"substitutions\">" << endl
            << "        <populationSize>" << endl
            << "            <parameter id=\"initialDemo.popSize\" value=\"100.0\"/>" << endl
            << "        </populationSize>" << endl
            << "    </constantSize>" << endl << endl
            << "    <!-- Generate a random starting tree under the coalescent process -->" << endl
            << "    <coalescentTree id=\"startingTree\" rootHeight=\"0.3\">" << endl
            << "        <taxa idref=\"taxa\"/>" << endl
            << "        <constantSize idref=\"initialDemo\"/>" << endl
            << "    </coalescentTree>" << endl << endl << endl;
        }
    } else if (treePrior == "concoal") {
        BEAST_xml_code
        << "<!-- *** PRIOR ON NODE AGES *** -->" << endl
        << "    <!-- A prior on the distribution node heights defined given a constant-size coalescent process -->" << endl
        << "    <constantSize id=\"constant\" units=\"substitutions\">" << endl
        << "    <populationSize>" << endl
        << "        <parameter id=\"constant.popSize\" value=\"0.3\" lower=\"0.0\" upper=\"Infinity\"/>" << endl
        << "    </populationSize>" << endl
        << "    </constantSize>" << endl << endl << endl;
        
        if (!starterTreePresent) {
            BEAST_xml_code
            << "    <!-- Generate a random starting tree under the coalescent process            -->" << endl
            << "    <coalescentTree id=\"startingTree\" rootHeight=\"0.3\">" << endl
            << "        <taxa idref=\"taxa\"/>" << endl
            << "        <constantSize idref=\"constant\"/>" << endl
            << "    </coalescentTree>" << endl << endl << endl;
        }
    } else if (treePrior == "expcoal") {
        BEAST_xml_code
        << "<!-- *** PRIOR ON NODE AGES *** -->" << endl
        << "    <!-- A prior on the distribution node heights defined given a exponential-growth coalescent process -->" << endl
        << "    <exponentialGrowth id=\"exponential\" units=\"substitutions\">" << endl
        << "    <populationSize>" << endl
        << "        <parameter id=\"exponential.popSize\" value=\"0.3\" lower=\"0.0\" upper=\"Infinity\"/>" << endl
        << "    </populationSize>" << endl
        << "    <growthRate>" << endl
        << "        <parameter id=\"exponential.growthRate\" value=\"3.0E-4\" lower=\"-Infinity\" upper=\"Infinity\"/>" << endl
        << "    </growthRate>" << endl
        << "    </exponentialGrowth>" << endl << endl << endl;
        
        if (!starterTreePresent) {
            BEAST_xml_code
            << "    <!-- Generate a random starting tree under the coalescent process            -->" << endl
            << "    <coalescentTree id=\"startingTree\" rootHeight=\"0.3\">" << endl
            << "        <taxa idref=\"taxa\"/>" << endl
            << "        <constantSize idref=\"exponential\"/>" << endl
            << "    </coalescentTree>" << endl << endl << endl;
        }
    } else if (treePrior == "logcoal") {
        BEAST_xml_code
        << "<!-- *** PRIOR ON NODE AGES *** -->" << endl
        << "    <!-- A prior on the distribution node heights defined given a logistic-growth coalescent process -->" << endl
        << "    <logisticGrowth id=\"logistic\" units=\"substitutions\">" << endl
        << "    <populationSize>" << endl
        << "        <parameter id=\"logistic.popSize\" value=\"0.3\" lower=\"0.0\" upper=\"Infinity\"/>" << endl
        << "    </populationSize>" << endl
        << "    <growthRate>" << endl
        << "        <parameter id=\"logistic.growthRate\" value=\"3.0E-4\" lower=\"-Infinity\" upper=\"Infinity\"/>" << endl
        << "    </growthRate>" << endl
        << "    <t50>" << endl
        << "        <parameter id=\"logistic.t50\" value=\"1.0\" lower=\"0.0\" upper=\"Infinity\"/>" << endl
        << "    </t50>" << endl
        << "    </logisticGrowth>" << endl << endl << endl;
        
        if (!starterTreePresent) {
            BEAST_xml_code
            << "    <constantSize id=\"initialDemo\" units=\"substitutions\">" << endl
            << "    <populationSize>" << endl
            << "        <parameter idref=\"logistic.popSize\"/>" << endl
            << "    </populationSize>" << endl
            << "    </constantSize>" << endl << endl
            << "    <!-- Generate a random starting tree under the coalescent process -->" << endl
            << "    <coalescentTree id=\"startingTree\" rootHeight=\"0.3\">" << endl
            << "        <taxa idref=\"taxa\"/>" << endl
            << "        <constantSize idref=\"initialDemo\"/>" << endl
            << "    </coalescentTree>" << endl << endl << endl;
        }
    }
        
    if (starterTreePresent) {
        BEAST_xml_code
        << "<!-- *** STARTING TREE - MUST BE COMPATIBLE WITH MONOPHYLY/TEMPORAL CONSTRAINTS OR OR INITIAL STATE OF MODEL WILL HAVE ZERO PROBABILITY *** -->" << endl
        << "<newick id=\"startingTree\" units=\"years\">" << endl
        << starterTree << endl
        << "</newick>" << endl << endl;
    }
}

void BEASTXML::writePartitionInformation (ofstream & BEAST_xml_code) {
    BEAST_xml_code << endl
    << "<!-- *** DEFINE PARTITIONS *** -->" << endl
    << "    <patterns id=\"patterns\" from=\"1\">" << endl
    << "        <alignment idref=\"alignment\"/>" << endl
    << "    </patterns>" << endl << endl;
}

void BEASTXML::writeTreeModel (ofstream & BEAST_xml_code, string const& treePrior) {
    BEAST_xml_code
    << "<!-- *** CONSTRUCT TREE MODEL *** -->" << endl
    << "    <treeModel id=\"treeModel\">" << endl
    << "        <tree idref=\"startingTree\"/>" << endl
    << "        <rootHeight>" << endl
    << "            <parameter id=\"treeModel.rootHeight\"/>" << endl
    << "        </rootHeight>" << endl
    << "        <nodeHeights internalNodes=\"true\">" << endl
    << "            <parameter id=\"treeModel.internalNodeHeights\"/>" << endl
    << "        </nodeHeights>" << endl
    << "        <nodeHeights internalNodes=\"true\" rootNode=\"true\">" << endl
    << "            <parameter id=\"treeModel.allInternalNodeHeights\"/>" << endl
    << "        </nodeHeights>" << endl
    << "    </treeModel>" << endl << endl;
    
    if (treePrior == "bd") {
        BEAST_xml_code
        << "    <speciationLikelihood id=\"speciation\">" << endl
        << "        <model>" << endl
        << "            <birthDeathModel idref=\"birthDeath\"/>" << endl
        << "        </model>" << endl
        << "        <speciesTree>" << endl
        << "            <treeModel idref=\"treeModel\"/>" << endl
        << "        </speciesTree>" << endl
        << "    </speciationLikelihood>" << endl << endl << endl;
    } else if (treePrior == "yule") {
        BEAST_xml_code
        << "    <speciationLikelihood id=\"speciation\">" << endl
        << "        <model>" << endl
        << "            <yuleModel idref=\"yule\"/>" << endl
        << "        </model>" << endl
        << "        <speciesTree>" << endl
        << "            <treeModel idref=\"treeModel\"/>" << endl
        << "        </speciesTree>" << endl
        << "    </speciationLikelihood>" << endl << endl << endl;
    } else if (treePrior == "concoal") {
        BEAST_xml_code
        << "    <!-- Generate a coalescent likelihood -->" << endl
        << "    <coalescentLikelihood id=\"coalescent\">" << endl
        << "        <model>" << endl
        << "            <constantSize idref=\"constant\"/>" << endl
        << "        </model>" << endl
        << "        <populationTree>" << endl
        << "            <treeModel idref=\"treeModel\"/>" << endl
        << "        </populationTree>" << endl
        << "    </coalescentLikelihood>" << endl << endl << endl;
    } else if (treePrior == "expcoal") {
        BEAST_xml_code
        << "    <!-- Generate a coalescent likelihood -->" << endl
        << "    <coalescentLikelihood id=\"coalescent\">" << endl
        << "        <model>" << endl
        << "            <exponentialGrowth idref=\"exponential\"/>" << endl
        << "        </model>" << endl
        << "        <populationTree>" << endl
        << "            <treeModel idref=\"treeModel\"/>" << endl
        << "        </populationTree>" << endl
        << "    </coalescentLikelihood>" << endl << endl << endl;
    } else if (treePrior == "logcoal") {
        BEAST_xml_code
        << "    <!-- Generate a coalescent likelihood -->" << endl
        << "    <coalescentLikelihood id=\"coalescent\">" << endl
        << "        <model>" << endl
        << "            <logisticGrowth idref=\"logistic\"/>" << endl
        << "        </model>" << endl
        << "        <populationTree>" << endl
        << "            <treeModel idref=\"treeModel\"/>" << endl
        << "        </populationTree>" << endl
        << "    </coalescentLikelihood>" << endl << endl << endl;
    }
}


void BEASTXML::writeClockModel (ofstream & BEAST_xml_code, int const& numTaxa, string const& clockFlavour) {
    BEAST_xml_code << "<!-- *** DEFINE CLOCK MODEL *** -->" << endl;

    if (clockFlavour == "ucln") {
        BEAST_xml_code
        << "    <!-- The uncorrelated relaxed clock (Drummond, Ho, Phillips & Rambaut, 2006) -->" << endl
        << "    <discretizedBranchRates id=\"branchRates\">" << endl
        << "        <treeModel idref=\"treeModel\"/>" << endl
        << "        <distribution>" << endl
        << "            <logNormalDistributionModel meanInRealSpace=\"true\">" << endl
        << "                <mean>" << endl
        << "                    <parameter id=\"ucld.mean\" value=\"0.001\" lower=\"0.0\" upper=\"10.0\"/>" << endl
        << "                </mean>" << endl
        << "                <stdev>" << endl
        << "                    <parameter id=\"ucld.stdev\" value=\"0.1\" lower=\"0.0\" upper=\"10.0\"/>" << endl
        << "                </stdev>" << endl
        << "            </logNormalDistributionModel>" << endl
        << "        </distribution>" << endl
        << "        <rateCategories>" << endl
        << "            <parameter id=\"branchRates.categories\" dimension=\"" << (2 * numTaxa) - 2 << "\"/>" << endl
        << "        </rateCategories>" << endl
        << "    </discretizedBranchRates>" << endl
        << endl
        << "    <rateStatistic id=\"meanRate\" name=\"meanRate\" mode=\"mean\" internal=\"true\" external=\"true\">" << endl
        << "        <treeModel idref=\"treeModel\"/>" << endl
        << "        <discretizedBranchRates idref=\"branchRates\"/>" << endl
        << "    </rateStatistic>" << endl
        << endl
        << "    <rateStatistic id=\"coefficientOfVariation\" name=\"coefficientOfVariation\" mode=\"coefficientOfVariation\" internal=\"true\" external=\"true\">" << endl
        << "        <treeModel idref=\"treeModel\"/>" << endl
        << "        <discretizedBranchRates idref=\"branchRates\"/>" << endl
        << "    </rateStatistic>" << endl
        << endl
        << "    <rateCovarianceStatistic id=\"covariance\" name=\"covariance\">" << endl
        << "        <treeModel idref=\"treeModel\"/>" << endl
        << "        <discretizedBranchRates idref=\"branchRates\"/>" << endl
        << "    </rateCovarianceStatistic>" << endl << endl;
    } else if (clockFlavour == "uced") {
        BEAST_xml_code
        << "    <!-- The uncorrelated relaxed clock (Drummond, Ho, Phillips & Rambaut, 2006) -->" << endl
        << "    <discretizedBranchRates id=\"branchRates\">" << endl
        << "        <treeModel idref=\"treeModel\"/>" << endl
        << "        <distribution>" << endl
        << "            <exponentialDistributionModel>" << endl
        << "                <mean>" << endl
        << "                    <parameter id=\"uced.mean\" value=\"1.0\" lower=\"0.0\" upper=\"1000000.0\"/>" << endl
        << "                </mean>" << endl
        << "            </exponentialDistributionModel>" << endl
        << "        </distribution>" << endl
        << "        <rateCategories>" << endl
        << "            <parameter id=\"branchRates.categories\" dimension=\"" << (2 * numTaxa) - 2 << "\"/>" << endl
        << "        </rateCategories>" << endl
        << "    </discretizedBranchRates>" << endl
        << "    <rateStatistic id=\"meanRate\" name=\"meanRate\" mode=\"mean\" internal=\"true\" external=\"true\">" << endl
        << "        <treeModel idref=\"treeModel\"/>" << endl
        << "        <discretizedBranchRates idref=\"branchRates\"/>" << endl
        << "    </rateStatistic>" << endl
        << "    <rateStatistic id=\"coefficientOfVariation\" name=\"coefficientOfVariation\" mode=\"coefficientOfVariation\" internal=\"true\" external=\"true\">" << endl
        << "        <treeModel idref=\"treeModel\"/>" << endl
        << "        <discretizedBranchRates idref=\"branchRates\"/>" << endl
        << "    </rateStatistic>" << endl
        << "    <rateCovarianceStatistic id=\"covariance\" name=\"covariance\">" << endl
        << "        <treeModel idref=\"treeModel\"/>" << endl
        << "        <discretizedBranchRates idref=\"branchRates\"/>" << endl
        << "    </rateCovarianceStatistic>" << endl << endl;
    } else if (clockFlavour == "strict") {
        BEAST_xml_code
        << "    <!-- The strict clock (Uniform rates across branches) -->" << endl
        << "        <strictClockBranchRates id=\"branchRates\">" << endl
        << "            <rate>" << endl
        << "                <parameter id=\"clock.rate\" value=\"1.0\"/>" << endl
        << "            </rate>" << endl
        << "        </strictClockBranchRates>" << endl << endl;
    } else if (clockFlavour == "randlocal") {
        BEAST_xml_code
        << "    <!-- The random local clock model (Drummond & Suchard, 2010) -->" << endl
        << "    <randomLocalClockModel id=\"branchRates\" ratesAreMultipliers=\"false\">" << endl
        << "        <treeModel idref=\"treeModel\"/>" << endl
        << "        <rates>" << endl
        << "            <parameter id=\"localClock.relativeRates\"/>" << endl
        << "        </rates>" << endl
        << "        <rateIndicator>" << endl
        << "            <parameter id=\"localClock.changes\"/>" << endl
        << "        </rateIndicator>" << endl
        << "        <clockRate>" << endl
        << "            <parameter id=\"clock.rate\" value=\"1.0\" lower=\"0.0\"/>" << endl
        << "        </clockRate>" << endl
        << "    </randomLocalClockModel>" << endl
        << "    <sumStatistic id=\"rateChanges\" name=\"rateChangeCount\" elementwise=\"true\">" << endl
        << "        <parameter idref=\"localClock.changes\"/>" << endl
        << "    </sumStatistic>" << endl
        << "    <rateStatistic id=\"meanRate\" name=\"meanRate\" mode=\"mean\" internal=\"true\" external=\"true\">" << endl
        << "        <treeModel idref=\"treeModel\"/>" << endl
        << "        <randomLocalClockModel idref=\"branchRates\"/>" << endl
        << "    </rateStatistic>" << endl
        << "    <rateStatistic id=\"coefficientOfVariation\" name=\"coefficientOfVariation\" mode=\"coefficientOfVariation\" internal=\"true\" external=\"true\">" << endl
        << "        <treeModel idref=\"treeModel\"/>" << endl
        << "        <randomLocalClockModel idref=\"branchRates\"/>" << endl
        << "    </rateStatistic>" << endl
        << "    <rateCovarianceStatistic id=\"covariance\" name=\"covariance\">" << endl
        << "        <treeModel idref=\"treeModel\"/>" << endl
        << "        <randomLocalClockModel idref=\"branchRates\"/>" << endl
        << "    </rateCovarianceStatistic>" << endl << endl;
    }
}


void BEASTXML::writeSubstitutionModel (ofstream & BEAST_xml_code,
    string const& partitionSubstitutionModel)
{
    BEAST_xml_code << endl
    << "<!-- *** DEFINE SUBSTITUTION MODEL -->"<< endl;
    if (partitionSubstitutionModel == "JC") {
        BEAST_xml_code
        << "    <!-- The JC69 substitution model (Jukes & Cantor, 1969) -->" << endl
        << "    <gtrModel id=\"JC69\">" << endl
        << "        <frequencies>" << endl
        << "            <frequencyModel dataType=\"nucleotide\">" << endl
        << "                <frequencies>" << endl
        << "                    <parameter id=\"JC69.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>" << endl
        << "                </frequencies>" << endl
        << "            </frequencyModel>" << endl
        << "        </frequencies>" << endl
        << "        <rateAC>" << endl
        << "            <parameter id=\"JC69.ac\" value=\"1.0\"/>" << endl
        << "        </rateAC>" << endl
        << "        <rateAG>" << endl
        << "            <parameter id=\"JC69.ag\" value=\"1.0\"/>" << endl
        << "        </rateAG>" << endl
        << "        <rateAT>" << endl
        << "            <parameter id=\"JC69.at\" value=\"1.0\"/>" << endl
        << "        </rateAT>" << endl
        << "        <rateCG>" << endl
        << "            <parameter id=\"JC69.cg\" value=\"1.0\"/>" << endl
        << "        </rateCG>" << endl
        << "        <rateGT>" << endl
        << "            <parameter id=\"JC69.gt\" value=\"1.0\"/>" << endl
        << "        </rateGT>" << endl
        << "    </gtrModel>" << endl << endl;
    } else if (partitionSubstitutionModel == "K80") {
        BEAST_xml_code
        << "    <!-- The K80 substitution model (Hasegawa, Kishino & Yano, 1985) -->" << endl
        << "    <hkyModel id=\"K80\">" << endl
        << "        <frequencies>" << endl
        << "            <frequencyModel dataType=\"nucleotide\">" << endl
        << "                <frequencies>" << endl
        << "                    <parameter id=\"K80.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>" << endl
        << "                </frequencies>" << endl
        << "            </frequencyModel>" << endl
        << "        </frequencies>" << endl
        << "        <kappa>" << endl
        << "            <parameter id=\"K80.kappa\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </kappa>" << endl
        << "    </hkyModel>" << endl << endl;
    } else if (partitionSubstitutionModel == "HKY") {
        BEAST_xml_code
        << "    <!-- The HKY substitution model (Hasegawa, Kishino & Yano, 1985) -->" << endl
        << "    <hkyModel id=\"HKY\">" << endl
        << "        <frequencies>" << endl
        << "            <frequencyModel dataType=\"nucleotide\">" << endl
        << "                <frequencies>" << endl
        << "                    <parameter id=\"HKY.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>" << endl
        << "                </frequencies>" << endl
        << "            </frequencyModel>" << endl
        << "        </frequencies>" << endl
        << "        <kappa>" << endl
        << "            <parameter id=\"HKY.kappa\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </kappa>" << endl
        << "    </hkyModel>" << endl << endl;
    } else if (partitionSubstitutionModel == "TrNef") {
        BEAST_xml_code            
        << "    <!-- The Tamura-Nei 1993 (TrNef) substitution model  -->" << endl
        << "    <gtrModel id=\"TrNef\">" << endl
        << "        <frequencies>" << endl
        << "            <frequencyModel dataType=\"nucleotide\">" << endl
        << "                <frequencies>" << endl
        << "                    <parameter id=\"TrNef.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>" << endl
        << "                </frequencies>" << endl
        << "            </frequencyModel>" << endl
        << "        </frequencies>" << endl
        << "        <rateAC>" << endl
        << "            <parameter id=\"TrNef.transversion\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateAC>" << endl
        << "        <rateAG>" << endl
        << "            <parameter id=\"TrNef.ag\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateAG>" << endl
        << "        <rateAT>" << endl
        << "            <parameter idref=\"TrNef.transversion\"/>" << endl
        << "        </rateAT>" << endl
        << "        <rateCG>" << endl
        << "            <parameter idref=\"TrNef.transversion\"/>" << endl
        << "        </rateCG>" << endl
        << "        <rateGT>" << endl
        << "            <parameter idref=\"TrNef.transversion\"/>" << endl
        << "        </rateGT>" << endl
        << "    </gtrModel>" << endl << endl;
    } else if (partitionSubstitutionModel == "TrN") {
        BEAST_xml_code            
        << "    <!-- The Tamura-Nei 1993 (TrN) substitution model  -->" << endl
        << "    <gtrModel id=\"TrN\">" << endl
        << "        <frequencies>" << endl
        << "            <frequencyModel dataType=\"nucleotide\">" << endl
        << "                <frequencies>" << endl
        << "                    <parameter id=\"TrN.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>" << endl
        << "                </frequencies>" << endl
        << "            </frequencyModel>" << endl
        << "        </frequencies>" << endl
        << "        <rateAC>" << endl
        << "            <parameter id=\"TrN.transversion\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateAC>" << endl
        << "        <rateAG>" << endl
        << "            <parameter id=\"TrN.ag\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateAG>" << endl
        << "        <rateAT>" << endl
        << "            <parameter idref=\"TrN.transversion\"/>" << endl
        << "        </rateAT>" << endl
        << "        <rateCG>" << endl
        << "            <parameter idref=\"TrN.transversion\"/>" << endl
        << "        </rateCG>" << endl
        << "        <rateGT>" << endl
        << "            <parameter idref=\"TrN.transversion\"/>" << endl
        << "        </rateGT>" << endl
        << "    </gtrModel>" << endl << endl;
    } else if (partitionSubstitutionModel == "K3P") {
        BEAST_xml_code
        << "    <!-- The Kimura 1981 3-parameter (K3P) substitution model -->" << endl
        << "    <gtrModel id=\"K3P\">" << endl
        << "        <frequencies>" << endl
        << "            <frequencyModel dataType=\"nucleotide\">" << endl
        << "                <frequencies>" << endl
        << "                    <parameter id=\"K3P.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>" << endl
        << "                </frequencies>" << endl
        << "            </frequencyModel>" << endl
        << "        </frequencies>" << endl
        << "        <rateAC>" << endl
        << "            <parameter id=\"K3P.purine2pyrimidine\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateAC>" << endl
        << "        <rateAG>" << endl
        << "            <parameter id=\"K3P.ag\" value=\"1.0\"/>" << endl
        << "        </rateAG>" << endl
        << "        <rateAT>" << endl
        << "            <parameter id=\"K3P.pyrimidine2purine\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateAT>" << endl
        << "        <rateCG>" << endl
        << "            <parameter idref=\"K3P.pyrimidine2purine\"/>" << endl
        << "        </rateCG>" << endl
        << "        <rateGT>" << endl
        << "            <parameter idref=\"K3P.purine2pyrimidine\"/>" << endl
        << "        </rateGT>" << endl
        << "    </gtrModel>" << endl << endl;
    } else if (partitionSubstitutionModel == "K3Puf") {
        BEAST_xml_code
        << "    <!-- The Kimura 1981 3-parameter (K3Puf) substitution model -->" << endl
        << "    <gtrModel id=\"K3Puf\">" << endl
        << "        <frequencies>" << endl
        << "            <frequencyModel dataType=\"nucleotide\">" << endl
        << "                <frequencies>" << endl
        << "                    <parameter id=\"K3Puf.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>" << endl
        << "                </frequencies>" << endl
        << "            </frequencyModel>" << endl
        << "        </frequencies>" << endl
        << "        <rateAC>" << endl
        << "            <parameter id=\"K3Puf.purine2pyrimidine\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateAC>" << endl
        << "        <rateAG>" << endl
        << "            <parameter id=\"K3Puf.ag\" value=\"1.0\"/>" << endl
        << "        </rateAG>" << endl
        << "        <rateAT>" << endl
        << "            <parameter id=\"K3Puf.pyrimidine2purine\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateAT>" << endl
        << "        <rateCG>" << endl
        << "            <parameter idref=\"K3Puf.pyrimidine2purine\"/>" << endl
        << "        </rateCG>" << endl
        << "        <rateGT>" << endl
        << "            <parameter idref=\"K3Puf.purine2pyrimidine\"/>" << endl
        << "        </rateGT>" << endl
        << "    </gtrModel>" << endl << endl;
    } else if (partitionSubstitutionModel == "TIMef") {
        BEAST_xml_code
        << "    <!-- The transitional (TIMef) substitution model -->" << endl
        << "    <gtrModel id=\"TIMef\">" << endl
        << "        <frequencies>" << endl
        << "            <frequencyModel dataType=\"nucleotide\">" << endl
        << "                <frequencies>" << endl
        << "                    <parameter id=\"TIMef.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>" << endl
        << "                </frequencies>" << endl
        << "            </frequencyModel>" << endl
        << "        </frequencies>" << endl
        << "        <rateAC>" << endl
        << "            <parameter id=\"TIMef.purine2pyrimidine\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateAC>" << endl
        << "        <rateAG>" << endl
        << "            <parameter id=\"TIMef.ag\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateAG>" << endl
        << "        <rateAT>" << endl
        << "            <parameter id=\"TIMef.pyrimidine2purine\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateAT>" << endl
        << "        <rateCG>" << endl
        << "            <parameter idref=\"TIMef.pyrimidine2purine\"/>" << endl
        << "        </rateCG>" << endl
        << "        <rateGT>" << endl
        << "            <parameter idref=\"TIMef.purine2pyrimidine\"/>" << endl
        << "        </rateGT>" << endl
        << "    </gtrModel>" << endl << endl;
    } else if (partitionSubstitutionModel == "TIM") {
        BEAST_xml_code
        << "    <!-- The transitional (TIM) substitution model -->" << endl
        << "    <gtrModel id=\"TIM\">" << endl
        << "        <frequencies>" << endl
        << "            <frequencyModel dataType=\"nucleotide\">" << endl
        << "                <frequencies>" << endl
        << "                    <parameter id=\"TIM.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>" << endl
        << "                </frequencies>" << endl
        << "            </frequencyModel>" << endl
        << "        </frequencies>" << endl
        << "        <rateAC>" << endl
        << "            <parameter id=\"TIM.purine2pyrimidine\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateAC>" << endl
        << "        <rateAG>" << endl
        << "            <parameter id=\"TIM.ag\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateAG>" << endl
        << "        <rateAT>" << endl
        << "            <parameter id=\"TIM.pyrimidine2purine\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateAT>" << endl
        << "        <rateCG>" << endl
        << "            <parameter idref=\"TIM.pyrimidine2purine\"/>" << endl
        << "        </rateCG>" << endl
        << "        <rateGT>" << endl
        << "            <parameter idref=\"TIM.purine2pyrimidine\"/>" << endl
        << "        </rateGT>" << endl
        << "    </gtrModel>" << endl << endl;
    } else if (partitionSubstitutionModel == "TVMef") {
        BEAST_xml_code
        << "    <!-- The transversional (TVMef) substitution model -->" << endl
        << "    <gtrModel id=\"TVMef\">" << endl
        << "        <frequencies>" << endl
        << "            <frequencyModel dataType=\"nucleotide\">" << endl
        << "                <frequencies>" << endl
        << "                    <parameter id=\"TVMef.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>" << endl
        << "                </frequencies>" << endl
        << "            </frequencyModel>" << endl
        << "        </frequencies>" << endl
        << "        <rateAC>" << endl
        << "            <parameter id=\"TVMef.ac\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateAC>" << endl
        << "        <rateAG>" << endl
        << "            <parameter id=\"TVMef.ag\" value=\"1.0\"/>" << endl
        << "        </rateAG>" << endl
        << "        <rateAT>" << endl
        << "            <parameter id=\"TVMef.at\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateAT>" << endl
        << "        <rateCG>" << endl
        << "            <parameter id=\"TVMef.cg\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateCG>" << endl
        << "        <rateGT>" << endl
        << "            <parameter id=\"TVMef.gt\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateGT>" << endl
        << "    </gtrModel>" << endl << endl;
    } else if (partitionSubstitutionModel == "TVM") {
        BEAST_xml_code
        << "    <!-- The transversional (TVM) substitution model  -->" << endl
        << "    <gtrModel id=\"TVM\">" << endl
        << "        <frequencies>" << endl
        << "            <frequencyModel dataType=\"nucleotide\">" << endl
        << "                <frequencies>" << endl
        << "                    <parameter id=\"TVM.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>" << endl
        << "                </frequencies>" << endl
        << "            </frequencyModel>" << endl
        << "        </frequencies>" << endl
        << "        <rateAC>" << endl
        << "            <parameter id=\"TVM.ac\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateAC>" << endl
        << "        <rateAG>" << endl
        << "            <parameter id=\"TVM.ag\" value=\"1.0\"/>" << endl
        << "        </rateAG>" << endl
        << "        <rateAT>" << endl
        << "            <parameter id=\"TVM.at\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateAT>" << endl
        << "        <rateCG>" << endl
        << "            <parameter id=\"TVM.cg\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateCG>" << endl
        << "        <rateGT>" << endl
        << "            <parameter id=\"TVM.gt\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateGT>" << endl
        << "    </gtrModel>" << endl << endl;
    } else if (partitionSubstitutionModel == "SYM") {
        BEAST_xml_code
        << "    <!-- The symmetric (SYM) substitution model -->" << endl
        << "    <gtrModel id=\"SYM\">" << endl
        << "        <frequencies>" << endl
        << "            <frequencyModel dataType=\"nucleotide\">" << endl
        << "                <frequencies>" << endl
        << "                    <parameter id=\"SYM.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>" << endl
        << "                </frequencies>" << endl
        << "            </frequencyModel>" << endl
        << "        </frequencies>" << endl
        << "        <rateAC>" << endl
        << "            <parameter id=\"SYM.ac\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateAC>" << endl
        << "        <rateAG>" << endl
        << "            <parameter id=\"SYM.ag\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateAG>" << endl
        << "        <rateAT>" << endl
        << "            <parameter id=\"SYM.at\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateAT>" << endl
        << "        <rateCG>" << endl
        << "            <parameter id=\"SYM.cg\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateCG>" << endl
        << "        <rateGT>" << endl
        << "            <parameter id=\"SYM.gt\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateGT>" << endl
        << "    </gtrModel>" << endl << endl;
    } else if (partitionSubstitutionModel == "GTR") {
        BEAST_xml_code
        << "    <!-- The general time reversible (GTR) substitution model -->" << endl
        << "    <gtrModel id=\"GTR\">" << endl
        << "        <frequencies>" << endl
        << "            <frequencyModel dataType=\"nucleotide\">" << endl
        << "                <frequencies>" << endl
        << "                    <parameter id=\"GTR.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>" << endl
        << "                </frequencies>" << endl
        << "            </frequencyModel>" << endl
        << "        </frequencies>" << endl
        << "        <rateAC>" << endl
        << "            <parameter id=\"GTR.ac\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateAC>" << endl
        << "        <rateAG>" << endl
        << "            <parameter id=\"GTR.ag\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateAG>" << endl
        << "        <rateAT>" << endl
        << "            <parameter id=\"GTR.at\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateAT>" << endl
        << "        <rateCG>" << endl
        << "            <parameter id=\"GTR.cg\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateCG>" << endl
        << "        <rateGT>" << endl
        << "            <parameter id=\"GTR.gt\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
        << "        </rateGT>" << endl
        << "    </gtrModel>" << endl << endl;
    }
    BEAST_xml_code << endl;
}

// not used at the moment
/*
void BEASTXML::writeSubstitutionModels (ofstream & BEAST_xml_code,
    vector <string> const& partitionSubstitutionModels)
{
    BEAST_xml_code << endl
    << "<!-- *** DEFINE SUBSTITUTION MODEL -->"<< endl;
    int counter = 0;
    for (vector <string>::const_iterator iter = partitionSubstitutionModels.begin(); iter != partitionSubstitutionModels.end(); iter++) {
        if (partitionSubstitutionModels[counter] == "JC") {
            BEAST_xml_code
            << "    <!-- The JC69 substitution model (Jukes & Cantor, 1969) -->" << endl
            << "    <gtrModel id=\"JC69\">" << endl
            << "        <frequencies>" << endl
            << "            <frequencyModel dataType=\"nucleotide\">" << endl
            << "                <frequencies>" << endl
            << "                    <parameter id=\"JC69.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>" << endl
            << "                </frequencies>" << endl
            << "            </frequencyModel>" << endl
            << "        </frequencies>" << endl
            << "        <rateAC>" << endl
            << "            <parameter id=\"JC69.ac\" value=\"1.0\"/>" << endl
            << "        </rateAC>" << endl
            << "        <rateAG>" << endl
            << "            <parameter id=\"JC69.ag\" value=\"1.0\"/>" << endl
            << "        </rateAG>" << endl
            << "        <rateAT>" << endl
            << "            <parameter id=\"JC69.at\" value=\"1.0\"/>" << endl
            << "        </rateAT>" << endl
            << "        <rateCG>" << endl
            << "            <parameter id=\"JC69.cg\" value=\"1.0\"/>" << endl
            << "        </rateCG>" << endl
            << "        <rateGT>" << endl
            << "            <parameter id=\"JC69.gt\" value=\"1.0\"/>" << endl
            << "        </rateGT>" << endl
            << "    </gtrModel>" << endl << endl;
        } else if (partitionSubstitutionModels[counter] == "K80") {
            BEAST_xml_code
            << "    <!-- The K80 substitution model (Hasegawa, Kishino & Yano, 1985) -->" << endl
            << "    <hkyModel id=\"K80\">" << endl
            << "        <frequencies>" << endl
            << "            <frequencyModel dataType=\"nucleotide\">" << endl
            << "                <frequencies>" << endl
            << "                    <parameter id=\"K80.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>" << endl
            << "                </frequencies>" << endl
            << "            </frequencyModel>" << endl
            << "        </frequencies>" << endl
            << "        <kappa>" << endl
            << "            <parameter id=\"K80.kappa\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </kappa>" << endl
            << "    </hkyModel>" << endl << endl;
        } else if (partitionSubstitutionModels[counter] == "HKY") {
            BEAST_xml_code
            << "    <!-- The HKY substitution model (Hasegawa, Kishino & Yano, 1985) -->" << endl
            << "    <hkyModel id=\"HKY\">" << endl
            << "        <frequencies>" << endl
            << "            <frequencyModel dataType=\"nucleotide\">" << endl
            << "                <frequencies>" << endl
            << "                    <parameter id=\"HKY.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>" << endl
            << "                </frequencies>" << endl
            << "            </frequencyModel>" << endl
            << "        </frequencies>" << endl
            << "        <kappa>" << endl
            << "            <parameter id=\"HKY.kappa\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </kappa>" << endl
            << "    </hkyModel>" << endl << endl;
        } else if (partitionSubstitutionModels[counter] == "TrNef") {
            BEAST_xml_code            
            << "    <!-- The Tamura-Nei 1993 (TrNef) substitution model  -->" << endl
            << "    <gtrModel id=\"TrNef\">" << endl
            << "        <frequencies>" << endl
            << "            <frequencyModel dataType=\"nucleotide\">" << endl
            << "                <frequencies>" << endl
            << "                    <parameter id=\"TrNef.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>" << endl
            << "                </frequencies>" << endl
            << "            </frequencyModel>" << endl
            << "        </frequencies>" << endl
            << "        <rateAC>" << endl
            << "            <parameter id=\"TrNef.transversion\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateAC>" << endl
            << "        <rateAG>" << endl
            << "            <parameter id=\"TrNef.ag\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateAG>" << endl
            << "        <rateAT>" << endl
            << "            <parameter idref=\"TrNef.transversion\"/>" << endl
            << "        </rateAT>" << endl
            << "        <rateCG>" << endl
            << "            <parameter idref=\"TrNef.transversion\"/>" << endl
            << "        </rateCG>" << endl
            << "        <rateGT>" << endl
            << "            <parameter idref=\"TrNef.transversion\"/>" << endl
            << "        </rateGT>" << endl
            << "    </gtrModel>" << endl << endl;
        } else if (partitionSubstitutionModels[counter] == "TrN") {
            BEAST_xml_code            
            << "    <!-- The Tamura-Nei 1993 (TrN) substitution model  -->" << endl
            << "    <gtrModel id=\"TrN\">" << endl
            << "        <frequencies>" << endl
            << "            <frequencyModel dataType=\"nucleotide\">" << endl
            << "                <frequencies>" << endl
            << "                    <parameter id=\"TrN.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>" << endl
            << "                </frequencies>" << endl
            << "            </frequencyModel>" << endl
            << "        </frequencies>" << endl
            << "        <rateAC>" << endl
            << "            <parameter id=\"TrN.transversion\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateAC>" << endl
            << "        <rateAG>" << endl
            << "            <parameter id=\"TrN.ag\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateAG>" << endl
            << "        <rateAT>" << endl
            << "            <parameter idref=\"TrN.transversion\"/>" << endl
            << "        </rateAT>" << endl
            << "        <rateCG>" << endl
            << "            <parameter idref=\"TrN.transversion\"/>" << endl
            << "        </rateCG>" << endl
            << "        <rateGT>" << endl
            << "            <parameter idref=\"TrN.transversion\"/>" << endl
            << "        </rateGT>" << endl
            << "    </gtrModel>" << endl << endl;
        } else if (partitionSubstitutionModels[counter] == "K3P") {
            BEAST_xml_code
            << "    <!-- The Kimura 1981 3-parameter (K3P) substitution model -->" << endl
            << "    <gtrModel id=\"K3P\">" << endl
            << "        <frequencies>" << endl
            << "            <frequencyModel dataType=\"nucleotide\">" << endl
            << "                <frequencies>" << endl
            << "                    <parameter id=\"K3P.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>" << endl
            << "                </frequencies>" << endl
            << "            </frequencyModel>" << endl
            << "        </frequencies>" << endl
            << "        <rateAC>" << endl
            << "            <parameter id=\"K3P.purine2pyrimidine\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateAC>" << endl
            << "        <rateAG>" << endl
            << "            <parameter id=\"K3P.ag\" value=\"1.0\"/>" << endl
            << "        </rateAG>" << endl
            << "        <rateAT>" << endl
            << "            <parameter id=\"K3P.pyrimidine2purine\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateAT>" << endl
            << "        <rateCG>" << endl
            << "            <parameter idref=\"K3P.pyrimidine2purine\"/>" << endl
            << "        </rateCG>" << endl
            << "        <rateGT>" << endl
            << "            <parameter idref=\"K3P.purine2pyrimidine\"/>" << endl
            << "        </rateGT>" << endl
            << "    </gtrModel>" << endl << endl;
        } else if (partitionSubstitutionModels[counter] == "K3Puf") {
            BEAST_xml_code
            << "    <!-- The Kimura 1981 3-parameter (K3Puf) substitution model -->" << endl
            << "    <gtrModel id=\"K3Puf\">" << endl
            << "        <frequencies>" << endl
            << "            <frequencyModel dataType=\"nucleotide\">" << endl
            << "                <frequencies>" << endl
            << "                    <parameter id=\"K3Puf.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>" << endl
            << "                </frequencies>" << endl
            << "            </frequencyModel>" << endl
            << "        </frequencies>" << endl
            << "        <rateAC>" << endl
            << "            <parameter id=\"K3Puf.purine2pyrimidine\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateAC>" << endl
            << "        <rateAG>" << endl
            << "            <parameter id=\"K3Puf.ag\" value=\"1.0\"/>" << endl
            << "        </rateAG>" << endl
            << "        <rateAT>" << endl
            << "            <parameter id=\"K3Puf.pyrimidine2purine\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateAT>" << endl
            << "        <rateCG>" << endl
            << "            <parameter idref=\"K3Puf.pyrimidine2purine\"/>" << endl
            << "        </rateCG>" << endl
            << "        <rateGT>" << endl
            << "            <parameter idref=\"K3Puf.purine2pyrimidine\"/>" << endl
            << "        </rateGT>" << endl
            << "    </gtrModel>" << endl << endl;
        } else if (partitionSubstitutionModels[counter] == "TIMef") {
            BEAST_xml_code
            << "    <!-- The transitional (TIMef) substitution model -->" << endl
            << "    <gtrModel id=\"TIMef\">" << endl
            << "        <frequencies>" << endl
            << "            <frequencyModel dataType=\"nucleotide\">" << endl
            << "                <frequencies>" << endl
            << "                    <parameter id=\"TIMef.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>" << endl
            << "                </frequencies>" << endl
            << "            </frequencyModel>" << endl
            << "        </frequencies>" << endl
            << "        <rateAC>" << endl
            << "            <parameter id=\"TIMef.purine2pyrimidine\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateAC>" << endl
            << "        <rateAG>" << endl
            << "            <parameter id=\"TIMef.ag\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateAG>" << endl
            << "        <rateAT>" << endl
            << "            <parameter id=\"TIMef.pyrimidine2purine\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateAT>" << endl
            << "        <rateCG>" << endl
            << "            <parameter idref=\"TIMef.pyrimidine2purine\"/>" << endl
            << "        </rateCG>" << endl
            << "        <rateGT>" << endl
            << "            <parameter idref=\"TIMef.purine2pyrimidine\"/>" << endl
            << "        </rateGT>" << endl
            << "    </gtrModel>" << endl << endl;
        } else if (partitionSubstitutionModels[counter] == "TIM") {
            BEAST_xml_code
            << "    <!-- The transitional (TIM) substitution model -->" << endl
            << "    <gtrModel id=\"TIM\">" << endl
            << "        <frequencies>" << endl
            << "            <frequencyModel dataType=\"nucleotide\">" << endl
            << "                <frequencies>" << endl
            << "                    <parameter id=\"TIM.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>" << endl
            << "                </frequencies>" << endl
            << "            </frequencyModel>" << endl
            << "        </frequencies>" << endl
            << "        <rateAC>" << endl
            << "            <parameter id=\"TIM.purine2pyrimidine\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateAC>" << endl
            << "        <rateAG>" << endl
            << "            <parameter id=\"TIM.ag\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateAG>" << endl
            << "        <rateAT>" << endl
            << "            <parameter id=\"TIM.pyrimidine2purine\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateAT>" << endl
            << "        <rateCG>" << endl
            << "            <parameter idref=\"TIM.pyrimidine2purine\"/>" << endl
            << "        </rateCG>" << endl
            << "        <rateGT>" << endl
            << "            <parameter idref=\"TIM.purine2pyrimidine\"/>" << endl
            << "        </rateGT>" << endl
            << "    </gtrModel>" << endl << endl;
        } else if (partitionSubstitutionModels[counter] == "TVMef") {
            BEAST_xml_code
            << "    <!-- The transversional (TVMef) substitution model -->" << endl
            << "    <gtrModel id=\"TVMef\">" << endl
            << "        <frequencies>" << endl
            << "            <frequencyModel dataType=\"nucleotide\">" << endl
            << "                <frequencies>" << endl
            << "                    <parameter id=\"TVMef.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>" << endl
            << "                </frequencies>" << endl
            << "            </frequencyModel>" << endl
            << "        </frequencies>" << endl
            << "        <rateAC>" << endl
            << "            <parameter id=\"TVMef.ac\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateAC>" << endl
            << "        <rateAG>" << endl
            << "            <parameter id=\"TVMef.ag\" value=\"1.0\"/>" << endl
            << "        </rateAG>" << endl
            << "        <rateAT>" << endl
            << "            <parameter id=\"TVMef.at\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateAT>" << endl
            << "        <rateCG>" << endl
            << "            <parameter id=\"TVMef.cg\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateCG>" << endl
            << "        <rateGT>" << endl
            << "            <parameter id=\"TVMef.gt\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateGT>" << endl
            << "    </gtrModel>" << endl << endl;
        } else if (partitionSubstitutionModels[counter] == "TVM") {
            BEAST_xml_code
            << "    <!-- The transversional (TVM) substitution model  -->" << endl
            << "    <gtrModel id=\"TVM\">" << endl
            << "        <frequencies>" << endl
            << "            <frequencyModel dataType=\"nucleotide\">" << endl
            << "                <frequencies>" << endl
            << "                    <parameter id=\"TVM.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>" << endl
            << "                </frequencies>" << endl
            << "            </frequencyModel>" << endl
            << "        </frequencies>" << endl
            << "        <rateAC>" << endl
            << "            <parameter id=\"TVM.ac\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateAC>" << endl
            << "        <rateAG>" << endl
            << "            <parameter id=\"TVM.ag\" value=\"1.0\"/>" << endl
            << "        </rateAG>" << endl
            << "        <rateAT>" << endl
            << "            <parameter id=\"TVM.at\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateAT>" << endl
            << "        <rateCG>" << endl
            << "            <parameter id=\"TVM.cg\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateCG>" << endl
            << "        <rateGT>" << endl
            << "            <parameter id=\"TVM.gt\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateGT>" << endl
            << "    </gtrModel>" << endl << endl;
        } else if (partitionSubstitutionModels[counter] == "SYM") {
            BEAST_xml_code
            << "    <!-- The symmetric (SYM) substitution model -->" << endl
            << "    <gtrModel id=\"SYM\">" << endl
            << "        <frequencies>" << endl
            << "            <frequencyModel dataType=\"nucleotide\">" << endl
            << "                <frequencies>" << endl
            << "                    <parameter id=\"SYM.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>" << endl
            << "                </frequencies>" << endl
            << "            </frequencyModel>" << endl
            << "        </frequencies>" << endl
            << "        <rateAC>" << endl
            << "            <parameter id=\"SYM.ac\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateAC>" << endl
            << "        <rateAG>" << endl
            << "            <parameter id=\"SYM.ag\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateAG>" << endl
            << "        <rateAT>" << endl
            << "            <parameter id=\"SYM.at\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateAT>" << endl
            << "        <rateCG>" << endl
            << "            <parameter id=\"SYM.cg\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateCG>" << endl
            << "        <rateGT>" << endl
            << "            <parameter id=\"SYM.gt\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateGT>" << endl
            << "    </gtrModel>" << endl << endl;
        } else if (partitionSubstitutionModels[counter] == "GTR") {
            BEAST_xml_code
            << "    <!-- The general time reversible (GTR) substitution model -->" << endl
            << "    <gtrModel id=\"GTR\">" << endl
            << "        <frequencies>" << endl
            << "            <frequencyModel dataType=\"nucleotide\">" << endl
            << "                <frequencies>" << endl
            << "                    <parameter id=\"GTR.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>" << endl
            << "                </frequencies>" << endl
            << "            </frequencyModel>" << endl
            << "        </frequencies>" << endl
            << "        <rateAC>" << endl
            << "            <parameter id=\"GTR.ac\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateAC>" << endl
            << "        <rateAG>" << endl
            << "            <parameter id=\"GTR.ag\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateAG>" << endl
            << "        <rateAT>" << endl
            << "            <parameter id=\"GTR.at\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateAT>" << endl
            << "        <rateCG>" << endl
            << "            <parameter id=\"GTR.cg\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateCG>" << endl
            << "        <rateGT>" << endl
            << "            <parameter id=\"GTR.gt\" value=\"1.0\" lower=\"1.0E-8\" upper=\"100.0\"/>" << endl
            << "        </rateGT>" << endl
            << "    </gtrModel>" << endl << endl;
        }
        counter++;
    }
    BEAST_xml_code << endl;
}
*/

void BEASTXML::writeSiteModel (ofstream & BEAST_xml_code, string const& partitionSubstitutionModel,
    string const& partitionSiteModel)
{
    BEAST_xml_code
    << "<!-- *** DEFINE AMONG-SITE HETEROGENEITY (SITE MODEL) *** -->" << endl;
    if (partitionSubstitutionModel == "JC") {
        BEAST_xml_code
        << "    <siteModel id=\"siteModel\">" << endl
        << "        <substitutionModel>" << endl
        << "            <gtrModel idref=\"JC69\"/>" << endl
        << "        </substitutionModel>" << endl;
    } else if (partitionSubstitutionModel == "K80") {
        BEAST_xml_code
        << "    <siteModel id=\"siteModel\">" << endl
        << "        <substitutionModel>" << endl
        << "            <hkyModel idref=\"K80\"/>" << endl
        << "        </substitutionModel>" << endl;
    } else if (partitionSubstitutionModel == "HKY") {
        BEAST_xml_code
        << "    <siteModel id=\"siteModel\">" << endl
        << "        <substitutionModel>" << endl
        << "            <hkyModel idref=\"HKY\"/>" << endl
        << "        </substitutionModel>" << endl;
    } else if (partitionSubstitutionModel == "TrNef") {
        BEAST_xml_code
        << "    <siteModel id=\"siteModel\">" << endl
        << "        <substitutionModel>" << endl
        << "            <gtrModel idref=\"TrNef\"/>" << endl
        << "        </substitutionModel>" << endl;
    } else if (partitionSubstitutionModel == "TrN") {
        BEAST_xml_code
        << "    <siteModel id=\"siteModel\">" << endl
        << "        <substitutionModel>" << endl
        << "            <gtrModel idref=\"TrN\"/>" << endl
        << "        </substitutionModel>" << endl;
    } else if (partitionSubstitutionModel == "K3P") {
        BEAST_xml_code
        << "    <siteModel id=\"siteModel\">" << endl
        << "        <substitutionModel>" << endl
        << "            <gtrModel idref=\"K3P\"/>" << endl
        << "        </substitutionModel>" << endl;
    } else if (partitionSubstitutionModel == "K3Puf") {
        BEAST_xml_code
        << "    <siteModel id=\"siteModel\">" << endl
        << "        <substitutionModel>" << endl
        << "            <gtrModel idref=\"K3Puf\"/>" << endl
        << "        </substitutionModel>" << endl;
    } else if (partitionSubstitutionModel == "TIMef") {
        BEAST_xml_code
        << "    <siteModel id=\"siteModel\">" << endl
        << "        <substitutionModel>" << endl
        << "            <gtrModel idref=\"TIMef\"/>" << endl
        << "        </substitutionModel>" << endl;
    } else if (partitionSubstitutionModel == "TIM") {
        BEAST_xml_code
        << "    <siteModel id=\"siteModel\">" << endl
        << "        <substitutionModel>" << endl
        << "            <gtrModel idref=\"TIM\"/>" << endl
        << "        </substitutionModel>" << endl;
    } else if (partitionSubstitutionModel == "TVMef") {
        BEAST_xml_code
        << "    <siteModel id=\"siteModel\">" << endl
        << "        <substitutionModel>" << endl
        << "            <gtrModel idref=\"TVMef\"/>" << endl
        << "        </substitutionModel>" << endl;
    } else if (partitionSubstitutionModel == "TVM") {
        BEAST_xml_code
        << "    <siteModel id=\"siteModel\">" << endl
        << "        <substitutionModel>" << endl
        << "            <gtrModel idref=\"TVM\"/>" << endl
        << "        </substitutionModel>" << endl;
    } else if (partitionSubstitutionModel == "SYM") {
        BEAST_xml_code
        << "    <siteModel id=\"siteModel\">" << endl
        << "        <substitutionModel>" << endl
        << "            <gtrModel idref=\"SYM\"/>" << endl
        << "        </substitutionModel>" << endl;
    } else if (partitionSubstitutionModel == "GTR") {
        BEAST_xml_code
        << "    <siteModel id=\"siteModel\">" << endl
        << "        <substitutionModel>" << endl
        << "            <gtrModel idref=\"GTR\"/>" << endl
        << "        </substitutionModel>" << endl;
    }
    
    if (partitionSiteModel == "IG" || partitionSiteModel == "G")
    {
        BEAST_xml_code
        << "        <gammaShape gammaCategories=\"4\">" << endl
        << "            <parameter id=\"alpha\" value=\"0.5\" lower=\"0.0\" upper=\"1000.0\"/>" << endl
        << "        </gammaShape>" << endl;
    }
    if (partitionSiteModel == "IG" || partitionSiteModel == "I")
    {
        BEAST_xml_code
        << "        <proportionInvariant>" << endl
        << "            <parameter id=\"pInv\" value=\"0.25\" lower=\"0.0\" upper=\"1.0\"/>" << endl
        << "        </proportionInvariant>" << endl;
    }
    BEAST_xml_code
    << "    </siteModel>" << endl << endl;
}

// not used at the moment
/*
void BEASTXML::writeSiteModels (ofstream & BEAST_xml_code, vector <string> const& partitionSubstitutionModels,
    vector <string> const& partitionSiteModels)
{
    int counter = 0;
    BEAST_xml_code
    << "<!-- *** DEFINE AMONG-SITE HETEROGENEITY (SITE MODEL) *** -->" << endl;
    for (vector<string>::const_iterator iter = partitionSubstitutionModels.begin(); iter != partitionSubstitutionModels.end(); iter++) {
        if (partitionSubstitutionModels[counter] == "JC") {
            BEAST_xml_code
            << "    <siteModel id=\"siteModel\">" << endl
            << "        <substitutionModel>" << endl
            << "            <gtrModel idref=\"JC69\"/>" << endl
            << "        </substitutionModel>" << endl;
        } else if (partitionSubstitutionModels[counter] == "K80") {
            BEAST_xml_code
            << "    <siteModel id=\"siteModel\">" << endl
            << "        <substitutionModel>" << endl
            << "            <hkyModel idref=\"K80\"/>" << endl
            << "        </substitutionModel>" << endl;
        } else if (partitionSubstitutionModels[counter] == "HKY") {
            BEAST_xml_code
            << "    <siteModel id=\"siteModel\">" << endl
            << "        <substitutionModel>" << endl
            << "            <hkyModel idref=\"HKY\"/>" << endl
            << "        </substitutionModel>" << endl;
        } else if (partitionSubstitutionModels[counter] == "TrNef") {
            BEAST_xml_code
            << "    <siteModel id=\"siteModel\">" << endl
            << "        <substitutionModel>" << endl
            << "            <gtrModel idref=\"TrNef\"/>" << endl
            << "        </substitutionModel>" << endl;
        } else if (partitionSubstitutionModels[counter] == "TrN") {
            BEAST_xml_code
            << "    <siteModel id=\"siteModel\">" << endl
            << "        <substitutionModel>" << endl
            << "            <gtrModel idref=\"TrN\"/>" << endl
            << "        </substitutionModel>" << endl;
        } else if (partitionSubstitutionModels[counter] == "K3P") {
            BEAST_xml_code
            << "    <siteModel id=\"siteModel\">" << endl
            << "        <substitutionModel>" << endl
            << "            <gtrModel idref=\"K3P\"/>" << endl
            << "        </substitutionModel>" << endl;
        } else if (partitionSubstitutionModels[counter] == "K3Puf") {
            BEAST_xml_code
            << "    <siteModel id=\"siteModel\">" << endl
            << "        <substitutionModel>" << endl
            << "            <gtrModel idref=\"K3Puf\"/>" << endl
            << "        </substitutionModel>" << endl;
        } else if (partitionSubstitutionModels[counter] == "TIMef") {
            BEAST_xml_code
            << "    <siteModel id=\"siteModel\">" << endl
            << "        <substitutionModel>" << endl
            << "            <gtrModel idref=\"TIMef\"/>" << endl
            << "        </substitutionModel>" << endl;
        } else if (partitionSubstitutionModels[counter] == "TIM") {
            BEAST_xml_code
            << "    <siteModel id=\"siteModel\">" << endl
            << "        <substitutionModel>" << endl
            << "            <gtrModel idref=\"TIM\"/>" << endl
            << "        </substitutionModel>" << endl;
        } else if (partitionSubstitutionModels[counter] == "TVMef") {
            BEAST_xml_code
            << "    <siteModel id=\"siteModel\">" << endl
            << "        <substitutionModel>" << endl
            << "            <gtrModel idref=\"TVMef\"/>" << endl
            << "        </substitutionModel>" << endl;
        } else if (partitionSubstitutionModels[counter] == "TVM") {
            BEAST_xml_code
            << "    <siteModel id=\"siteModel\">" << endl
            << "        <substitutionModel>" << endl
            << "            <gtrModel idref=\"TVM\"/>" << endl
            << "        </substitutionModel>" << endl;
        } else if (partitionSubstitutionModels[counter] == "SYM") {
            BEAST_xml_code
            << "    <siteModel id=\"siteModel\">" << endl
            << "        <substitutionModel>" << endl
            << "            <gtrModel idref=\"SYM\"/>" << endl
            << "        </substitutionModel>" << endl;
        } else if (partitionSubstitutionModels[counter] == "GTR") {
            BEAST_xml_code
            << "    <siteModel id=\"siteModel\">" << endl
            << "        <substitutionModel>" << endl
            << "            <gtrModel idref=\"GTR\"/>" << endl
            << "        </substitutionModel>" << endl;
        }
        
        if (partitionSiteModels[counter] == "IG" || partitionSiteModels[counter] == "G")
        {
            BEAST_xml_code
            << "        <gammaShape gammaCategories=\"4\">" << endl
            << "            <parameter id=\"alpha\" value=\"0.5\" lower=\"0.0\" upper=\"1000.0\"/>" << endl
            << "        </gammaShape>" << endl;
        }
        if (partitionSiteModels[counter] == "IG" || partitionSiteModels[counter] == "I")
        {
            BEAST_xml_code
            << "        <proportionInvariant>" << endl
            << "            <parameter id=\"pInv\" value=\"0.25\" lower=\"0.0\" upper=\"1.0\"/>" << endl
            << "        </proportionInvariant>" << endl;
        }
        BEAST_xml_code
        << "    </siteModel>" << endl << endl;
        counter++;
    }
}
*/

void BEASTXML::writeTreeLikelihoods (ofstream & BEAST_xml_code, string const& clockFlavour) {
    BEAST_xml_code << endl
    << "<!-- *** DEFINE TREE LIKELIHOOD *** -->" << endl
    << "    <treeLikelihood id=\"treeLikelihood\">" << endl
    << "        <patterns idref=\"patterns\"/>" << endl
    << "        <treeModel idref=\"treeModel\"/>" << endl
    << "        <siteModel idref=\"siteModel\"/>" << endl;
    if (clockFlavour == "ucln" || clockFlavour == "uced") {
        BEAST_xml_code
        << "        <discretizedBranchRates idref=\"branchRates\"/>" << endl;
    } else if (clockFlavour == "strict") {
        BEAST_xml_code
        << "        <strictClockBranchRates idref=\"branchRates\"/>" << endl;
    } else if (clockFlavour == "randlocal") {
        BEAST_xml_code
        << "        <randomLocalClockModel idref=\"branchRates\"/>" << endl;
    }
    BEAST_xml_code
    << "    </treeLikelihood>" << endl << endl;
}

void BEASTXML::writeOperators (ofstream & BEAST_xml_code, string const& treePrior, bool const& manipulateTreeTopology,
    string const& partitionSubstitutionModel, string const& partitionSiteModel, string const& clockFlavour, int const& numTaxa)
{
    BEAST_xml_code << endl
    << "<!-- *** DEFINE OPERATORS *** -->" << endl
    <<     "    <operators id=\"operators\">" << endl;
    if (partitionSubstitutionModel == "K80") {
        BEAST_xml_code    
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"K80.kappa\"/>" << endl
        << "        </scaleOperator>" << endl;
    } else if (partitionSubstitutionModel == "HKY") {
        BEAST_xml_code    
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"HKY.kappa\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <deltaExchange delta=\"0.01\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"HKY.frequencies\"/>" << endl
        << "        </deltaExchange>" << endl;
    } else if (partitionSubstitutionModel == "TrNef") {
        BEAST_xml_code
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"TrNef.transversion\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"TrNef.ag\"/>" << endl
        << "        </scaleOperator>" << endl;
    } else if (partitionSubstitutionModel == "TrN") {
        BEAST_xml_code
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"TrN.transversion\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"TrN.ag\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <deltaExchange delta=\"0.01\" weight=\"1\">" << endl
        << "            <parameter idref=\"TrN.frequencies\"/>" << endl
        << "        </deltaExchange>" << endl;
    } else if (partitionSubstitutionModel == "K3P") {
        BEAST_xml_code
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"K3P.purine2pyrimidine\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"K3P.pyrimidine2purine\"/>" << endl
        << "        </scaleOperator>" << endl;
    } else if (partitionSubstitutionModel == "K3Puf") {
        BEAST_xml_code
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"K3Puf.purine2pyrimidine\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"K3Puf.pyrimidine2purine\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <deltaExchange delta=\"0.01\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"K3Puf.frequencies\"/>" << endl
        << "        </deltaExchange>" << endl;
    } else if (partitionSubstitutionModel == "TIMef") {
        BEAST_xml_code
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"TIMef.purine2pyrimidine\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"TIMef.ag\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"TIMef.pyrimidine2purine\"/>" << endl
        << "        </scaleOperator>" << endl;
    } else if (partitionSubstitutionModel == "TIM") {
        BEAST_xml_code
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"TIM.purine2pyrimidine\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"TIM.ag\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"TIM.pyrimidine2purine\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <deltaExchange delta=\"0.01\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"TIM.frequencies\"/>" << endl
        << "        </deltaExchange>" << endl;
    } else if (partitionSubstitutionModel == "TVMef") {
        BEAST_xml_code
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"TVMef.ac\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"TVMef.at\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"TVMef.cg\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"TVMef.gt\"/>" << endl
        << "        </scaleOperator>" << endl;
    } else if (partitionSubstitutionModel == "TVM") {
        BEAST_xml_code
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"TVM.ac\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"TVM.at\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"TVM.cg\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"TVM.gt\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <deltaExchange delta=\"0.01\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"TVM.frequencies\"/>" << endl
        << "        </deltaExchange>" << endl;
    } else if (partitionSubstitutionModel == "SYM") {
        BEAST_xml_code
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"SYM.ac\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"SYM.ag\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"SYM.at\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"SYM.cg\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"SYM.gt\"/>" << endl
        << "        </scaleOperator>" << endl;;
    } else if (partitionSubstitutionModel == "GTR") {
        BEAST_xml_code
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"GTR.ac\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"GTR.ag\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"GTR.at\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"GTR.cg\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"GTR.gt\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <deltaExchange delta=\"0.01\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"GTR.frequencies\"/>" << endl
        << "        </deltaExchange>" << endl;
    }
    if (partitionSiteModel == "IG" || partitionSiteModel == "G")
    {
        BEAST_xml_code
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"alpha\"/>" << endl
        << "        </scaleOperator>" << endl;
    }
    if (partitionSiteModel == "IG" || partitionSiteModel == "I")
    {
        BEAST_xml_code
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
        << "            <parameter idref=\"pInv\"/>" << endl
        << "        </scaleOperator>" << endl;
    }
    BEAST_xml_code << endl;
    
    BEAST_xml_code
    << "<!-- *** CLOCK OPERATORS *** -->" << endl;
    
// NOTE: removed randomWalkIntegerOperator operator on branch rates. Apparently not necessary.
    
    if (clockFlavour == "ucln") {
        BEAST_xml_code
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"5\">" << endl
        << "            <parameter idref=\"ucld.mean\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"5\">" << endl
        << "            <parameter idref=\"ucld.stdev\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <swapOperator size=\"1\" weight=\"" << numTaxa << "\" autoOptimize=\"false\">" << endl
        << "            <parameter idref=\"branchRates.categories\"/>" << endl
        << "        </swapOperator>" << endl
//         << "        <randomWalkIntegerOperator windowSize=\"1\" weight=\"10\">" << endl
//         << "            <parameter idref=\"branchRates.categories\"/>" << endl
//         << "        </randomWalkIntegerOperator>" << endl
        << "        <uniformIntegerOperator weight=\"10\">" << endl
        << "            <parameter idref=\"branchRates.categories\"/>" << endl
        << "        </uniformIntegerOperator>" << endl
        << "        <upDownOperator scaleFactor=\"0.75\" weight=\"5\">" << endl
        << "            <up>" << endl
        << "                <compoundParameter idref=\"ucld.mean\"/>" << endl
         << "            </up>" << endl
        << "            <down>" << endl
        << "                <parameter idref=\"treeModel.allInternalNodeHeights\"/>" << endl
        << "            </down>" << endl
        << "        </upDownOperator>" << endl << endl;
    } else if (clockFlavour == "uced") {
        BEAST_xml_code
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"5\">" << endl
        << "            <parameter idref=\"uced.mean\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <swapOperator size=\"1\" weight=\"" << numTaxa << "\" autoOptimize=\"false\">" << endl
        << "            <parameter idref=\"branchRates.categories\"/>" << endl
        << "        </swapOperator>" << endl
//         << "        <randomWalkIntegerOperator windowSize=\"1\" weight=\"10\">" << endl
//         << "            <parameter idref=\"branchRates.categories\"/>" << endl
//         << "        </randomWalkIntegerOperator>" << endl
        << "        <uniformIntegerOperator weight=\"10\">" << endl
        << "            <parameter idref=\"branchRates.categories\"/>" << endl
        << "        </uniformIntegerOperator>" << endl
        << "        <upDownOperator scaleFactor=\"0.75\" weight=\"5\">" << endl
        << "            <up>" << endl
        << "                <parameter idref=\"uced.mean\"/>" << endl
        << "            </up>" << endl
        << "            <down>" << endl
        << "                <parameter idref=\"treeModel.allInternalNodeHeights\"/>" << endl
        << "            </down>" << endl
        << "        </upDownOperator>" << endl << endl;
    } else if (clockFlavour == "randlocal") {
        BEAST_xml_code
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"5\">" << endl
        << "            <parameter idref=\"clock.rate\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"15\">" << endl
        << "            <parameter idref=\"localClock.relativeRates\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <bitFlipOperator weight=\"15\">" << endl
        << "            <parameter idref=\"localClock.changes\"/>" << endl
        << "        </bitFlipOperator>" << endl
        << "        <upDownOperator scaleFactor=\"0.75\" weight=\"5\">" << endl
        << "            <up>" << endl
        << "                <parameter idref=\"clock.rate\"/>" << endl
        << "            </up>" << endl
        << "            <down>" << endl
        << "                <parameter idref=\"treeModel.allInternalNodeHeights\"/>" << endl
        << "            </down>" << endl
        << "        </upDownOperator>" << endl << endl;
    }
    
    if (treePrior == "bd") {
        BEAST_xml_code
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"3\">" << endl
        << "            <parameter idref=\"birthDeath.BminusDRate\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"3\">" << endl
        << "            <parameter idref=\"birthDeath.DoverB\"/>" << endl
        << "        </scaleOperator>" << endl << endl;
    } else if (treePrior == "yule") {
        BEAST_xml_code
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"3\">" << endl
        << "            <parameter idref=\"yule.birthRate\"/>" << endl
        << "        </scaleOperator>" << endl;
    } else if (treePrior == "concoal") {
        BEAST_xml_code
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"3\">" << endl
        << "            <parameter idref=\"constant.popSize\"/>" << endl
        << "        </scaleOperator>" << endl;
    } else if (treePrior == "expcoal") {
        BEAST_xml_code
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"3\">" << endl
        << "            <parameter idref=\"exponential.popSize\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <randomWalkOperator windowSize=\"1.0\" weight=\"3\">" << endl
        << "            <parameter idref=\"exponential.growthRate\"/>" << endl
        << "        </randomWalkOperator>" << endl
        << "        <randomWalkOperator windowSize=\"1.0\" weight=\"3\" boundaryCondition=\"absorbing\">" << endl
        << "            <parameter idref=\"exponential.growthRate\"/>" << endl
        << "        </randomWalkOperator>" << endl;
    } else if (treePrior == "logcoal") {
        BEAST_xml_code
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"3\">" << endl
        << "            <parameter idref=\"logistic.popSize\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <randomWalkOperator windowSize=\"1.0\" weight=\"3\">" << endl
        << "            <parameter idref=\"logistic.growthRate\"/>" << endl
        << "        </randomWalkOperator>" << endl
        << "        <randomWalkOperator windowSize=\"1.0\" weight=\"3\" boundaryCondition=\"absorbing\">" << endl
        << "            <parameter idref=\"logistic.growthRate\"/>" << endl
        << "        </randomWalkOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"3\">" << endl
        << "            <parameter idref=\"logistic.t50\"/>" << endl
        << "        </scaleOperator>" << endl;
    }
    
    BEAST_xml_code
    << "        <scaleOperator scaleFactor=\"0.75\" weight=\"5\">" << endl
    << "            <parameter idref=\"treeModel.rootHeight\"/>" << endl
    << "        </scaleOperator>" << endl
    
    << "        <uniformOperator weight=\"30\">" << endl
    << "            <parameter idref=\"treeModel.internalNodeHeights\"/>" << endl
    << "        </uniformOperator>" << endl << endl;
    
    if (manipulateTreeTopology) {
        BEAST_xml_code
        << "<!-- *** TOPOLOGY MANIPULATION OPERATORS *** -->" << endl
        << "        <subtreeSlide size=\"0.02\" gaussian=\"true\" weight=\"15\">" << endl
        << "            <treeModel idref=\"treeModel\"/>" << endl
        << "        </subtreeSlide>" << endl
        << "        <narrowExchange weight=\"15\">" << endl
        << "            <treeModel idref=\"treeModel\"/>" << endl
        << "        </narrowExchange>" << endl
        << "        <wideExchange weight=\"3\">" << endl
        << "            <treeModel idref=\"treeModel\"/>" << endl
        << "        </wideExchange>" << endl
        << "        <wilsonBalding weight=\"3\">" << endl
        << "            <treeModel idref=\"treeModel\"/>" << endl
        << "        </wilsonBalding>" << endl;
    } else if (!manipulateTreeTopology) {
        BEAST_xml_code << endl
        << "<!-- *** TOPOLOGY MANIPULATION OPERATORS *** -->" << endl
        << "<!-- *** Commented out for current analysis *** -->" << endl
        << "<!--" << endl
        << "        <subtreeSlide size=\"0.02\" gaussian=\"true\" weight=\"15\">" << endl
        << "            <treeModel idref=\"treeModel\"/>" << endl
        << "        </subtreeSlide>" << endl
        << "        <narrowExchange weight=\"150\">" << endl
        << "            <treeModel idref=\"treeModel\"/>" << endl
        << "        </narrowExchange>" << endl
        << "        <wideExchange weight=\"30\">" << endl
        << "            <treeModel idref=\"treeModel\"/>" << endl
        << "        </wideExchange>" << endl
        << "        <wilsonBalding weight=\"30\">" << endl
        << "            <treeModel idref=\"treeModel\"/>" << endl
        << "        </wilsonBalding>" << endl
        << "-->" << endl;
    }
    
    BEAST_xml_code
    << "    </operators>" << endl << endl;
}

// not used at the moment
/*
void BEASTXML::writeOperators (ofstream & BEAST_xml_code, string const& treePrior, bool const& manipulateTreeTopology,
    vector <string> const& partitionSubstitutionModels, vector <string> const& partitionSiteModels,
    string const& clockFlavour, int const& numTaxa)
{
    int counter = 0;
    BEAST_xml_code << endl
    << "<!-- *** DEFINE OPERATORS *** -->" << endl
    <<     "    <operators id=\"operators\">" << endl;
    for (vector<string>::const_iterator iter = partitionSubstitutionModels.begin(); iter != partitionSubstitutionModels.end(); iter++) {
        if (partitionSubstitutionModels[counter] == "K80") {
            BEAST_xml_code    
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"K80.kappa\"/>" << endl
            << "        </scaleOperator>" << endl;
        } else if (partitionSubstitutionModels[counter] == "HKY") {
            BEAST_xml_code    
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"HKY.kappa\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <deltaExchange delta=\"0.01\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"HKY.frequencies\"/>" << endl
            << "        </deltaExchange>" << endl;
        } else if (partitionSubstitutionModels[counter] == "TrNef") {
            BEAST_xml_code
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"TrNef.transversion\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"TrNef.ag\"/>" << endl
            << "        </scaleOperator>" << endl;
        } else if (partitionSubstitutionModels[counter] == "TrN") {
            BEAST_xml_code
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"TrN.transversion\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"TrN.ag\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <deltaExchange delta=\"0.01\" weight=\"1\">" << endl
            << "            <parameter idref=\"TrN.frequencies\"/>" << endl
            << "        </deltaExchange>" << endl;
        } else if (partitionSubstitutionModels[counter] == "K3P") {
            BEAST_xml_code
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"K3P.purine2pyrimidine\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"K3P.pyrimidine2purine\"/>" << endl
            << "        </scaleOperator>" << endl;
        } else if (partitionSubstitutionModels[counter] == "K3Puf") {
            BEAST_xml_code
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"K3Puf.purine2pyrimidine\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"K3Puf.pyrimidine2purine\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <deltaExchange delta=\"0.01\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"K3Puf.frequencies\"/>" << endl
            << "        </deltaExchange>" << endl;
        } else if (partitionSubstitutionModels[counter] == "TIMef") {
            BEAST_xml_code
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"TIMef.purine2pyrimidine\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"TIMef.ag\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"TIMef.pyrimidine2purine\"/>" << endl
            << "        </scaleOperator>" << endl;
        } else if (partitionSubstitutionModels[counter] == "TIM") {
            BEAST_xml_code
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"TIM.purine2pyrimidine\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"TIM.ag\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"TIM.pyrimidine2purine\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <deltaExchange delta=\"0.01\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"TIM.frequencies\"/>" << endl
            << "        </deltaExchange>" << endl;
        } else if (partitionSubstitutionModels[counter] == "TVMef") {
            BEAST_xml_code
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"TVMef.ac\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"TVMef.at\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"TVMef.cg\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"TVMef.gt\"/>" << endl
            << "        </scaleOperator>" << endl;
        } else if (partitionSubstitutionModels[counter] == "TVM") {
            BEAST_xml_code
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"TVM.ac\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"TVM.at\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"TVM.cg\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"TVM.gt\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <deltaExchange delta=\"0.01\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"TVM.frequencies\"/>" << endl
            << "        </deltaExchange>" << endl;
        } else if (partitionSubstitutionModels[counter] == "SYM") {
            BEAST_xml_code
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"SYM.ac\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"SYM.ag\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"SYM.at\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"SYM.cg\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"SYM.gt\"/>" << endl
            << "        </scaleOperator>" << endl;;
        } else if (partitionSubstitutionModels[counter] == "GTR") {
            BEAST_xml_code
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"GTR.ac\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"GTR.ag\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"GTR.at\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"GTR.cg\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"GTR.gt\"/>" << endl
            << "        </scaleOperator>" << endl
            << "        <deltaExchange delta=\"0.01\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"GTR.frequencies\"/>" << endl
            << "        </deltaExchange>" << endl;
        }
        if (partitionSiteModels[counter] == "IG" || partitionSiteModels[counter] == "G")
        {
            BEAST_xml_code
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"alpha\"/>" << endl
            << "        </scaleOperator>" << endl;
        }
        if (partitionSiteModels[counter] == "IG" || partitionSiteModels[counter] == "I")
        {
            BEAST_xml_code
            << "        <scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">" << endl
            << "            <parameter idref=\"pInv\"/>" << endl
            << "        </scaleOperator>" << endl;
        }
        BEAST_xml_code << endl;
        counter++;
    }
    
    BEAST_xml_code
    << "<!-- *** CLOCK OPERATORS *** -->" << endl;
    
    if (clockFlavour == "ucln") {
        BEAST_xml_code
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"5\">" << endl
        << "            <parameter idref=\"ucld.mean\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"5\">" << endl
        << "            <parameter idref=\"ucld.stdev\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <swapOperator size=\"1\" weight=\"" << numTaxa << "\" autoOptimize=\"false\">" << endl
        << "            <parameter idref=\"branchRates.categories\"/>" << endl
        << "        </swapOperator>" << endl
        << "        <randomWalkIntegerOperator windowSize=\"1\" weight=\"10\">" << endl
        << "            <parameter idref=\"branchRates.categories\"/>" << endl
        << "        </randomWalkIntegerOperator>" << endl
        << "        <uniformIntegerOperator weight=\"10\">" << endl
        << "            <parameter idref=\"branchRates.categories\"/>" << endl
        << "        </uniformIntegerOperator>" << endl
        << "        <upDownOperator scaleFactor=\"0.75\" weight=\"5\" autoOptimize=\"true\">" << endl
        << "            <up>" << endl
        << "                <compoundParameter idref=\"ucld.mean\"/>" << endl
         << "            </up>" << endl
        << "            <down>" << endl
        << "                <parameter idref=\"treeModel.allInternalNodeHeights\"/>" << endl
        << "            </down>" << endl
        << "        </upDownOperator>" << endl << endl;
    } else if (clockFlavour == "uced") {
        BEAST_xml_code
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"3\">" << endl
        << "            <parameter idref=\"uced.mean\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <swapOperator size=\"1\" weight=\"10\" autoOptimize=\"false\">" << endl
        << "            <parameter idref=\"branchRates.categories\"/>" << endl
        << "        </swapOperator>" << endl
        << "        <randomWalkIntegerOperator windowSize=\"1\" weight=\"10\">" << endl
        << "            <parameter idref=\"branchRates.categories\"/>" << endl
        << "        </randomWalkIntegerOperator>" << endl
        << "        <uniformIntegerOperator weight=\"10\">" << endl
        << "            <parameter idref=\"branchRates.categories\"/>" << endl
        << "        </uniformIntegerOperator>" << endl
        << "        <upDownOperator scaleFactor=\"0.75\" weight=\"3\">" << endl
        << "            <up>" << endl
        << "                <parameter idref=\"uced.mean\"/>" << endl
        << "            </up>" << endl
        << "            <down>" << endl
        << "                <parameter idref=\"treeModel.allInternalNodeHeights\"/>" << endl
        << "            </down>" << endl
        << "        </upDownOperator>" << endl << endl;
    }
    
    if (treePrior == "bd") {
        BEAST_xml_code
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"3\">" << endl
        << "            <parameter idref=\"birthDeath.BminusDRate\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"3\">" << endl
        << "            <parameter idref=\"birthDeath.DoverB\"/>" << endl
        << "        </scaleOperator>" << endl << endl;
    } else if (treePrior == "yule") {
        BEAST_xml_code
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"3\">" << endl
        << "            <parameter idref=\"yule.birthRate\"/>" << endl
        << "        </scaleOperator>" << endl;
    } else if (treePrior == "concoal") {
        BEAST_xml_code
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"3\">" << endl
        << "            <parameter idref=\"constant.popSize\"/>" << endl
        << "        </scaleOperator>" << endl;
    } else if (treePrior == "expcoal") {
        BEAST_xml_code
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"3\">" << endl
        << "            <parameter idref=\"exponential.popSize\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <randomWalkOperator windowSize=\"1.0\" weight=\"3\">" << endl
        << "            <parameter idref=\"exponential.growthRate\"/>" << endl
        << "        </randomWalkOperator>" << endl
        << "        <randomWalkOperator windowSize=\"1.0\" weight=\"3\" boundaryCondition=\"absorbing\">" << endl
        << "            <parameter idref=\"exponential.growthRate\"/>" << endl
        << "        </randomWalkOperator>" << endl;
    } else if (treePrior == "logcoal") {
        BEAST_xml_code
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"3\">" << endl
        << "            <parameter idref=\"logistic.popSize\"/>" << endl
        << "        </scaleOperator>" << endl
        << "        <randomWalkOperator windowSize=\"1.0\" weight=\"3\">" << endl
        << "            <parameter idref=\"logistic.growthRate\"/>" << endl
        << "        </randomWalkOperator>" << endl
        << "        <randomWalkOperator windowSize=\"1.0\" weight=\"3\" boundaryCondition=\"absorbing\">" << endl
        << "            <parameter idref=\"logistic.growthRate\"/>" << endl
        << "        </randomWalkOperator>" << endl
        << "        <scaleOperator scaleFactor=\"0.75\" weight=\"3\">" << endl
        << "            <parameter idref=\"logistic.t50\"/>" << endl
        << "        </scaleOperator>" << endl;
    }
    
    BEAST_xml_code
    << "        <scaleOperator scaleFactor=\"0.75\" weight=\"5\">" << endl
    << "            <parameter idref=\"treeModel.rootHeight\"/>" << endl
    << "        </scaleOperator>" << endl
    
    << "        <uniformOperator weight=\"30\">" << endl
    << "            <parameter idref=\"treeModel.internalNodeHeights\"/>" << endl
    << "        </uniformOperator>" << endl << endl;
    
    if (manipulateTreeTopology) {
        BEAST_xml_code
        << "<!-- *** TOPOLOGY MANIPULATION OPERATORS *** -->" << endl
        << "        <subtreeSlide size=\"0.02\" gaussian=\"true\" weight=\"15\">" << endl
        << "            <treeModel idref=\"treeModel\"/>" << endl
        << "        </subtreeSlide>" << endl
        << "        <narrowExchange weight=\"15\">" << endl
        << "            <treeModel idref=\"treeModel\"/>" << endl
        << "        </narrowExchange>" << endl
        << "        <wideExchange weight=\"3\">" << endl
        << "            <treeModel idref=\"treeModel\"/>" << endl
        << "        </wideExchange>" << endl
        << "        <wilsonBalding weight=\"3\">" << endl
        << "            <treeModel idref=\"treeModel\"/>" << endl
        << "        </wilsonBalding>" << endl;
    } else if (!manipulateTreeTopology) {
        BEAST_xml_code << endl
        << "<!-- *** TOPOLOGY MANIPULATION OPERATORS *** -->" << endl
        << "<!-- *** Commented out for current analysis *** -->" << endl
        << "<!--" << endl
        << "        <subtreeSlide size=\"0.02\" gaussian=\"true\" weight=\"15\">" << endl
        << "            <treeModel idref=\"treeModel\"/>" << endl
        << "        </subtreeSlide>" << endl
        << "        <narrowExchange weight=\"150\">" << endl
        << "            <treeModel idref=\"treeModel\"/>" << endl
        << "        </narrowExchange>" << endl
        << "        <wideExchange weight=\"30\">" << endl
        << "            <treeModel idref=\"treeModel\"/>" << endl
        << "        </wideExchange>" << endl
        << "        <wilsonBalding weight=\"30\">" << endl
        << "            <treeModel idref=\"treeModel\"/>" << endl
        << "        </wilsonBalding>" << endl
        << "-->" << endl;
    }
    
    BEAST_xml_code
    << "    </operators>" << endl << endl;
}
*/

void BEASTXML::writeMCMCParameters (ofstream & BEAST_xml_code, int const& mcmcLength, string const& clockFlavour, 
    string const& partitionSubstitutionModel, vector <string> const& rootPrior, string const& treePrior)
{
    BEAST_xml_code << endl
    << "<!-- *** MCMC PARAMETERS *** -->" << endl
    << "    <mcmc id=\"mcmc\" chainLength=\"" << mcmcLength << "\" autoOptimize=\"true\">" << endl
    << "        <posterior id=\"posterior\">" << endl
    << "            <prior id=\"prior\">" << endl;
    
// set tree priors. currently uses improper clock priors
    if (treePrior == "bd") {
        BEAST_xml_code
        << "            <!-- Improper uniform clock prior -->" << endl
        << "                <uniformPrior lower=\"0.0\" upper=\"100000.0\">" << endl
        << "                    <parameter idref=\"birthDeath.BminusDRate\"/>" << endl
        << "                </uniformPrior>" << endl
        << "                <uniformPrior lower=\"0.0\" upper=\"1.0\">" << endl
        << "                    <parameter idref=\"birthDeath.DoverB\"/>" << endl
        << "                </uniformPrior>" << endl
        << "                <speciationLikelihood idref=\"speciation\"/>" << endl << endl;
    } else if (treePrior == "yule") {
        BEAST_xml_code
        << "            <!-- Improper uniform clock prior -->" << endl
        << "                <uniformPrior lower=\"0.0\" upper=\"1.0E100\">" << endl
        << "                    <parameter idref=\"yule.birthRate\"/>" << endl
        << "                </uniformPrior>" << endl
        << "                <speciationLikelihood idref=\"speciation\"/>" << endl << endl;
    } else if (treePrior == "concoal") {
        BEAST_xml_code
        << "                <oneOnXPrior>" << endl
        << "                    <parameter idref=\"constant.popSize\"/>" << endl
        << "                </oneOnXPrior>" << endl
        << "                <coalescentLikelihood idref=\"coalescent\"/>" << endl << endl;
    } else if (treePrior == "expcoal") {
        BEAST_xml_code
        << "                <oneOnXPrior>" << endl
        << "                    <parameter idref=\"exponential.popSize\"/>" << endl
        << "                </oneOnXPrior>" << endl
        << "                <laplacePrior mean=\"0.0\" scale=\"30.701134573253945\">" << endl
        << "                    <parameter idref=\"exponential.growthRate\"/>" << endl
        << "                </laplacePrior>" << endl
        << "                <coalescentLikelihood idref=\"coalescent\"/>" << endl << endl;
    } else if (treePrior == "logcoal") {
        BEAST_xml_code
        << "                <oneOnXPrior>" << endl
        << "                    <parameter idref=\"logistic.popSize\"/>" << endl
        << "                </oneOnXPrior>" << endl
        << "                <laplacePrior mean=\"0.0\" scale=\"23.025850929940457\">" << endl
        << "                    <parameter idref=\"logistic.growthRate\"/>" << endl
        << "                </laplacePrior>" << endl
        << "                <gammaPrior shape=\"0.0010\" scale=\"1000.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"logistic.t50\"/>" << endl
        << "                </gammaPrior>" << endl
        << "                <coalescentLikelihood idref=\"coalescent\"/>" << endl << endl;
    }
    
// set clock priors
// *** allow this to be set by user ***
    if (clockFlavour == "ucln") {
        BEAST_xml_code
        << "                <exponentialPrior mean=\"0.3333333333333333\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"ucld.stdev\"/>" << endl
        << "                </exponentialPrior>" << endl << endl;
    } else if (clockFlavour == "randlocal") {
        BEAST_xml_code
        << "                <poissonPrior mean=\"0.6931471805599453\" offset=\"0.0\">" << endl
        << "                    <statistic idref=\"rateChanges\"/>" << endl
        << "                </poissonPrior>" << endl
        << "                <gammaPrior shape=\"0.5\" scale=\"2.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"localClock.relativeRates\"/>" << endl
        << "                </gammaPrior>" << endl << endl;
    }
    
// set root prior (if present)
    if (!rootPrior.empty()) {
        if (rootPrior[0] == "unif") {
            BEAST_xml_code
            << "                <uniformPrior lower=\"" << rootPrior[1] << "\" upper=\"" << rootPrior[2] << "\">" << endl
            << "                    <parameter idref=\"treeModel.rootHeight\"/>" << endl
            << "                </uniformPrior>" << endl << endl;
        } else if (rootPrior[0] == "norm") {
            BEAST_xml_code
            << "                <normalPrior mean=\"" << rootPrior[1] << "\" stdev=\"" << rootPrior[2] << "\">" << endl
            << "                    <parameter idref=\"treeModel.rootHeight\"/>" << endl
            << "                </normalPrior>" << endl << endl;
        } else {
            cout << "Shit. Fucked up somewhere..." << endl;
            exit(1);
        }
    }
    
// set substitution model priors
    if (partitionSubstitutionModel == "K80") {
        BEAST_xml_code
        << "                <logNormalPrior mean=\"1.0\" stdev=\"1.25\" offset=\"0.0\" meanInRealSpace=\"false\">" << endl
        << "                    <parameter idref=\"k80.kappa\"/>" << endl
        << "                </logNormalPrior>" << endl;
    } else if (partitionSubstitutionModel == "HKY") {
// Updated prior here to logNormal
        BEAST_xml_code
        << "                <logNormalPrior mean=\"1.0\" stdev=\"1.25\" offset=\"0.0\" meanInRealSpace=\"false\">" << endl
        << "                    <parameter idref=\"HKY.kappa\"/>" << endl
        << "                </logNormalPrior>" << endl;
    } else if (partitionSubstitutionModel == "TrNef") {
        BEAST_xml_code
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"TrNef.transversion\"/>" << endl
        << "                </gammaPrior>" << endl
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"TrNef.ag\"/>" << endl
        << "                </gammaPrior>" << endl;
    } else if (partitionSubstitutionModel == "TrN") {
        BEAST_xml_code
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"TrN.transversion\"/>" << endl
        << "                </gammaPrior>" << endl
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"TrN.ag\"/>" << endl
        << "                </gammaPrior>" << endl;
    } else if (partitionSubstitutionModel == "K3P") {
        BEAST_xml_code
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"K3P.purine2pyrimidine\"/>" << endl
        << "                </gammaPrior>" << endl
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"K3P.pyrimidine2purine\"/>" << endl
        << "                </gammaPrior>" << endl;
    } else if (partitionSubstitutionModel == "K3Puf") {
        BEAST_xml_code
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"K3Puf.purine2pyrimidine\"/>" << endl
        << "                </gammaPrior>" << endl
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"K3Puf.pyrimidine2purine\"/>" << endl
        << "                </gammaPrior>" << endl;
    } else if (partitionSubstitutionModel == "TIMef") {
        BEAST_xml_code
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"TIMef.purine2pyrimidine\"/>" << endl
        << "                </gammaPrior>" << endl
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"TIMef.ag\"/>" << endl
        << "                </gammaPrior>" << endl
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"TIMef.pyrimidine2purine\"/>" << endl
        << "                </gammaPrior>" << endl;
    } else if (partitionSubstitutionModel == "TIM") {
        BEAST_xml_code
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"TIM.purine2pyrimidine\"/>" << endl
        << "                </gammaPrior>" << endl
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"TIM.ag\"/>" << endl
        << "                </gammaPrior>" << endl
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"TIM.pyrimidine2purine\"/>" << endl
        << "                </gammaPrior>" << endl;
    } else if (partitionSubstitutionModel == "TVMef") {
        BEAST_xml_code
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"TVMef.ac\"/>" << endl
        << "                </gammaPrior>" << endl
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"TVMef.at\"/>" << endl
        << "                </gammaPrior>" << endl
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"TVMef.cg\"/>" << endl
        << "                </gammaPrior>" << endl
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"TVMef.gt\"/>" << endl
        << "                </gammaPrior>" << endl;
    } else if (partitionSubstitutionModel == "TVM") {
        BEAST_xml_code
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"TVM.ac\"/>" << endl
        << "                </gammaPrior>" << endl
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"TVM.at\"/>" << endl
        << "                </gammaPrior>" << endl
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"TVM.cg\"/>" << endl
        << "                </gammaPrior>" << endl
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"TVM.gt\"/>" << endl
        << "                </gammaPrior>" << endl;
    } else if (partitionSubstitutionModel == "SYM") {
        BEAST_xml_code
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"SYM.ac\"/>" << endl
        << "                </gammaPrior>" << endl
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"SYM.ag\"/>" << endl
        << "                </gammaPrior>" << endl
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"SYM.at\"/>" << endl
        << "                </gammaPrior>" << endl
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"SYM.cg\"/>" << endl
        << "                </gammaPrior>" << endl
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"SYM.gt\"/>" << endl
        << "                </gammaPrior>" << endl;
    } else if (partitionSubstitutionModel == "GTR") {
        BEAST_xml_code
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"GTR.ac\"/>" << endl
        << "                </gammaPrior>" << endl
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"GTR.ag\"/>" << endl
        << "                </gammaPrior>" << endl
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"GTR.at\"/>" << endl
        << "                </gammaPrior>" << endl
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"GTR.cg\"/>" << endl
        << "                </gammaPrior>" << endl
        << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"GTR.gt\"/>" << endl
        << "                </gammaPrior>" << endl;
    }
    
    BEAST_xml_code
    << "            </prior>" << endl;
    
    BEAST_xml_code
    << "            <likelihood id=\"likelihood\">" << endl
    << "                <treeLikelihood idref=\"treeLikelihood\"/>" << endl
    << "            </likelihood>" << endl
    << "        </posterior>" << endl
    << "        <operators idref=\"operators\"/>" << endl;
}

// not used at the moment
/*
void BEASTXML::writeMCMCParameters (ofstream & BEAST_xml_code, int const& mcmcLength, string const& clockFlavour, 
    vector <string> const& partitionSubstitutionModels, vector <string> const& rootPrior, string const& treePrior)
{
    BEAST_xml_code << endl
    << "<!-- *** MCMC PARAMETERS *** -->" << endl
    << "    <mcmc id=\"mcmc\" chainLength=\"" << mcmcLength << "\" autoOptimize=\"true\">" << endl
    << "        <posterior id=\"posterior\">" << endl
    << "            <prior id=\"prior\">" << endl;
    
    if (treePrior == "bd" || treePrior == "yule") {
        BEAST_xml_code
        << "                <speciationLikelihood idref=\"speciation\"/>" << endl << endl;
    } else if (treePrior == "concoal") {
        BEAST_xml_code
        << "                <oneOnXPrior>" << endl
        << "                    <parameter idref=\"constant.popSize\"/>" << endl
        << "                </oneOnXPrior>" << endl
        << "                <coalescentLikelihood idref=\"coalescent\"/>" << endl << endl;
    } else if (treePrior == "expcoal") {
        BEAST_xml_code
        << "                <oneOnXPrior>" << endl
        << "                    <parameter idref=\"exponential.popSize\"/>" << endl
        << "                </oneOnXPrior>" << endl
        << "                <laplacePrior mean=\"0.0\" scale=\"30.701134573253945\">" << endl
        << "                    <parameter idref=\"exponential.growthRate\"/>" << endl
        << "                </laplacePrior>" << endl
        << "                <coalescentLikelihood idref=\"coalescent\"/>" << endl << endl;
    } else if (treePrior == "logcoal") {
        BEAST_xml_code
        << "                <oneOnXPrior>" << endl
        << "                    <parameter idref=\"logistic.popSize\"/>" << endl
        << "                </oneOnXPrior>" << endl
        << "                <laplacePrior mean=\"0.0\" scale=\"23.025850929940457\">" << endl
        << "                    <parameter idref=\"logistic.growthRate\"/>" << endl
        << "                </laplacePrior>" << endl
        << "                <gammaPrior shape=\"0.0010\" scale=\"1000.0\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"logistic.t50\"/>" << endl
        << "                </gammaPrior>" << endl
        << "                <coalescentLikelihood idref=\"coalescent\"/>" << endl << endl;
    }
    
    if (clockFlavour == "ucln") {
        BEAST_xml_code
        << "                <exponentialPrior mean=\"0.3333333333333333\" offset=\"0.0\">" << endl
        << "                    <parameter idref=\"ucld.stdev\"/>" << endl
        << "                </exponentialPrior>" << endl << endl;
    }
    if (!rootPrior.empty()) {
        if (rootPrior[0] == "unif") {
            BEAST_xml_code
            << "                <uniformPrior lower=\"" << rootPrior[1] << "\" upper=\"" << rootPrior[2] << "\">" << endl
            << "                    <parameter idref=\"treeModel.rootHeight\"/>" << endl
            << "                </uniformPrior>" << endl << endl;
        } else if (rootPrior[0] == "norm") {
            BEAST_xml_code
            << "                <normalPrior mean=\"" << rootPrior[1] << "\" stdev=\"" << rootPrior[2] << "\">" << endl
            << "                    <parameter idref=\"treeModel.rootHeight\"/>" << endl
            << "                </normalPrior>" << endl << endl;
        } else {
            cout << "Shit. Fucked up somewhere..." << endl;
            exit(1);
        }
    }
    
    int counter = 0;
    for (vector<string>::const_iterator iter = partitionSubstitutionModels.begin(); iter != partitionSubstitutionModels.end(); iter++) {
        if (partitionSubstitutionModels[counter] == "K80") {
            BEAST_xml_code
            << "                <logNormalPrior mean=\"1.0\" stdev=\"1.25\" offset=\"0.0\" meanInRealSpace=\"false\">" << endl
            << "                    <parameter idref=\"k80.kappa\"/>" << endl
            << "                </logNormalPrior>" << endl;
        } else if (partitionSubstitutionModels[counter] == "HKY") {
// Updated prior here to logNormal
            BEAST_xml_code
            << "                <logNormalPrior mean=\"1.0\" stdev=\"1.25\" offset=\"0.0\" meanInRealSpace=\"false\">" << endl
            << "                    <parameter idref=\"HKY.kappa\"/>" << endl
            << "                </logNormalPrior>" << endl;
        } else if (partitionSubstitutionModels[counter] == "TrNef") {
            BEAST_xml_code
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"TrNef.transversion\"/>" << endl
            << "                </gammaPrior>" << endl
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"TrNef.ag\"/>" << endl
            << "                </gammaPrior>" << endl;
        } else if (partitionSubstitutionModels[counter] == "TrN") {
            BEAST_xml_code
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"TrN.transversion\"/>" << endl
            << "                </gammaPrior>" << endl
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"TrN.ag\"/>" << endl
            << "                </gammaPrior>" << endl;
        } else if (partitionSubstitutionModels[counter] == "K3P") {
            BEAST_xml_code
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"K3P.purine2pyrimidine\"/>" << endl
            << "                </gammaPrior>" << endl
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"K3P.pyrimidine2purine\"/>" << endl
            << "                </gammaPrior>" << endl;
        } else if (partitionSubstitutionModels[counter] == "K3Puf") {
            BEAST_xml_code
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"K3Puf.purine2pyrimidine\"/>" << endl
            << "                </gammaPrior>" << endl
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"K3Puf.pyrimidine2purine\"/>" << endl
            << "                </gammaPrior>" << endl;
        } else if (partitionSubstitutionModels[counter] == "TIMef") {
            BEAST_xml_code
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"TIMef.purine2pyrimidine\"/>" << endl
            << "                </gammaPrior>" << endl
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"TIMef.ag\"/>" << endl
            << "                </gammaPrior>" << endl
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"TIMef.pyrimidine2purine\"/>" << endl
            << "                </gammaPrior>" << endl;
        } else if (partitionSubstitutionModels[counter] == "TIM") {
            BEAST_xml_code
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"TIM.purine2pyrimidine\"/>" << endl
            << "                </gammaPrior>" << endl
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"TIM.ag\"/>" << endl
            << "                </gammaPrior>" << endl
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"TIM.pyrimidine2purine\"/>" << endl
            << "                </gammaPrior>" << endl;
        } else if (partitionSubstitutionModels[counter] == "TVMef") {
            BEAST_xml_code
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"TVMef.ac\"/>" << endl
            << "                </gammaPrior>" << endl
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"TVMef.at\"/>" << endl
            << "                </gammaPrior>" << endl
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"TVMef.cg\"/>" << endl
            << "                </gammaPrior>" << endl
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"TVMef.gt\"/>" << endl
            << "                </gammaPrior>" << endl;
        } else if (partitionSubstitutionModels[counter] == "TVM") {
            BEAST_xml_code
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"TVM.ac\"/>" << endl
            << "                </gammaPrior>" << endl
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"TVM.at\"/>" << endl
            << "                </gammaPrior>" << endl
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"TVM.cg\"/>" << endl
            << "                </gammaPrior>" << endl
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"TVM.gt\"/>" << endl
            << "                </gammaPrior>" << endl;
        } else if (partitionSubstitutionModels[counter] == "SYM") {
            BEAST_xml_code
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"SYM.ac\"/>" << endl
            << "                </gammaPrior>" << endl
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"SYM.ag\"/>" << endl
            << "                </gammaPrior>" << endl
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"SYM.at\"/>" << endl
            << "                </gammaPrior>" << endl
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"SYM.cg\"/>" << endl
            << "                </gammaPrior>" << endl
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"SYM.gt\"/>" << endl
            << "                </gammaPrior>" << endl;
        } else if (partitionSubstitutionModels[counter] == "GTR") {
            BEAST_xml_code
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"GTR.ac\"/>" << endl
            << "                </gammaPrior>" << endl
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"GTR.ag\"/>" << endl
            << "                </gammaPrior>" << endl
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"GTR.at\"/>" << endl
            << "                </gammaPrior>" << endl
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"GTR.cg\"/>" << endl
            << "                </gammaPrior>" << endl
            << "                <gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">" << endl
            << "                    <parameter idref=\"GTR.gt\"/>" << endl
            << "                </gammaPrior>" << endl;
        }
        counter++;
    }
    BEAST_xml_code
    << "            </prior>" << endl;
    
    BEAST_xml_code
    << "            <likelihood id=\"likelihood\">" << endl
    << "                <treeLikelihood idref=\"treeLikelihood\"/>" << endl
    << "            </likelihood>" << endl
    << "        </posterior>" << endl
    << "        <operators idref=\"operators\"/>" << endl;
}
*/

void BEASTXML::writeScreenLog (ofstream & BEAST_xml_code, int const& screenSampling, string const& clockFlavour) {
    BEAST_xml_code << endl
    << "<!-- *** PRINT PARAMETERS TO SCREEN *** -->" << endl
    << "        <log id=\"screenLog\" logEvery=\"" << screenSampling << "\">" << endl
    << "            <column label=\"Posterior\" dp=\"4\" width=\"12\">" << endl
    << "                <posterior idref=\"posterior\"/>" << endl
    << "            </column>" << endl
    << "            <column label=\"Prior\" dp=\"4\" width=\"12\">" << endl
    << "                <prior idref=\"prior\"/>" << endl
    << "            </column>" << endl
    << "            <column label=\"Likelihood\" dp=\"4\" width=\"12\">" << endl
    << "                <likelihood idref=\"likelihood\"/>" << endl
    << "            </column>" << endl
    << "            <column label=\"Root Height\" sf=\"6\" width=\"12\">" << endl
    << "                <parameter idref=\"treeModel.rootHeight\"/>" << endl
    << "            </column>" << endl;
    
    if (clockFlavour == "ucln") {
        BEAST_xml_code
        << "            <column label=\"Rate\" sf=\"6\" width=\"12\">" << endl
        << "                <parameter idref=\"ucld.mean\"/>" << endl
        << "            </column>" << endl;
    } else if (clockFlavour == "uced") {
        BEAST_xml_code
        << "            <column label=\"Rate\" sf=\"6\" width=\"12\">" << endl
        << "                <parameter idref=\"uced.mean\"/>" << endl
        << "            </column>" << endl;
    } else if (clockFlavour == "strict" || clockFlavour == "randlocal") {
        BEAST_xml_code
        << "            <column label=\"Rate\" sf=\"6\" width=\"12\">" << endl
        << "                <parameter idref=\"clock.rate\"/>" << endl
        << "            </column>" << endl;
    }
    BEAST_xml_code
    << "        </log>" << endl << endl;
}

void BEASTXML::writeParameterLog (ofstream & BEAST_xml_code, int const& parameterSampling,
    string const& treePrior, string & clockFlavour, string const& partitionSubstitutionModel,
    string const& partitionSiteModel)
{
    string prunedFileName = getRootName(XMLOutFileName);
    BEAST_xml_code
    << "<!-- *** PRINT PARAMETERS TO FILE *** -->" << endl
    << "        <log id=\"fileLog\" logEvery=\"" << parameterSampling << "\" fileName=\"" << prunedFileName << ".log\">" << endl
    << "            <posterior idref=\"posterior\"/>" << endl
    << "            <prior idref=\"prior\"/>" << endl
    << "            <likelihood idref=\"likelihood\"/>" << endl
    << "            <parameter idref=\"treeModel.rootHeight\"/>" << endl;
    
    if (treePrior == "bd") {
        BEAST_xml_code
        << "            <parameter idref=\"birthDeath.BminusDRate\"/>" << endl
        << "            <parameter idref=\"birthDeath.DoverB\"/>" << endl << endl;
    } else if (treePrior == "yule") {
        BEAST_xml_code
        << "            <parameter idref=\"yule.birthRate\"/>" << endl << endl;
    } else if (treePrior == "concoal") {
        BEAST_xml_code
        << "            <parameter idref=\"constant.popSize\"/>" << endl << endl;
    } else if (treePrior == "expcoal") {
        BEAST_xml_code
        << "            <parameter idref=\"exponential.popSize\"/>" << endl
        << "            <parameter idref=\"exponential.growthRate\"/>" << endl << endl;
    } else if (treePrior == "logcoal") {
        BEAST_xml_code
        << "            <parameter idref=\"logistic.popSize\"/>" << endl
        << "            <parameter idref=\"logistic.growthRate\"/>" << endl
        << "            <parameter idref=\"logistic.t50\"/>" << endl << endl;
    }
    
    if (partitionSubstitutionModel == "K80") {
        BEAST_xml_code
        << "            <parameter idref=\"K80.kappa\"/>" << endl
        << "            <parameter idref=\"K80.frequencies\"/>" << endl;
    } else if (partitionSubstitutionModel == "HKY") {
        BEAST_xml_code
        << "            <parameter idref=\"HKY.kappa\"/>" << endl
        << "            <parameter idref=\"HKY.frequencies\"/>" << endl;
    } else if (partitionSubstitutionModel == "TrNef") {
        BEAST_xml_code
        << "            <parameter idref=\"TrNef.transversion\"/>" << endl
        << "            <parameter idref=\"TrNef.ag\"/>" << endl
        << "            <parameter idref=\"TrNef.frequencies\"/>" << endl;
    } else if (partitionSubstitutionModel == "TrN") {
        BEAST_xml_code
        << "            <parameter idref=\"TrN.transversion\"/>" << endl
        << "            <parameter idref=\"TrN.ag\"/>" << endl
        << "            <parameter idref=\"TrN.frequencies\"/>" << endl;
    } else if (partitionSubstitutionModel == "K3P") {
        BEAST_xml_code
        << "            <parameter idref=\"K3P.purine2pyrimidine\"/>" << endl
        << "            <parameter idref=\"K3P.pyrimidine2purine\"/>" << endl
        << "            <parameter idref=\"K3P.frequencies\"/>" << endl;
    } else if (partitionSubstitutionModel == "K3Puf") {
        BEAST_xml_code
        << "            <parameter idref=\"K3Puf.purine2pyrimidine\"/>" << endl
        << "            <parameter idref=\"K3Puf.pyrimidine2purine\"/>" << endl
        << "            <parameter idref=\"K3Puf.frequencies\"/>" << endl;
    } else if (partitionSubstitutionModel == "TIMef") {
        BEAST_xml_code
        << "            <parameter idref=\"TIMef.purine2pyrimidine\"/>" << endl
        << "            <parameter idref=\"TIMef.ag\"/>" << endl
        << "            <parameter idref=\"TIMef.pyrimidine2purine\"/>" << endl
        << "            <parameter idref=\"TIMef.frequencies\"/>" << endl;
    } else if (partitionSubstitutionModel == "TIM") {
        BEAST_xml_code
        << "            <parameter idref=\"TIM.purine2pyrimidine\"/>" << endl
        << "            <parameter idref=\"TIM.ag\"/>" << endl
        << "            <parameter idref=\"TIM.pyrimidine2purine\"/>" << endl
        << "            <parameter idref=\"TIM.frequencies\"/>" << endl;
    } else if (partitionSubstitutionModel == "TVMef") {
        BEAST_xml_code
        << "            <parameter idref=\"TVMef.ac\"/>" << endl
        << "            <parameter idref=\"TVMef.at\"/>" << endl
        << "            <parameter idref=\"TVMef.cg\"/>" << endl
        << "            <parameter idref=\"TVMef.gt\"/>" << endl
        << "            <parameter idref=\"TVMef.frequencies\"/>" << endl;
    } else if (partitionSubstitutionModel == "TVM") {
        BEAST_xml_code
        << "            <parameter idref=\"TVM.ac\"/>" << endl
        << "            <parameter idref=\"TVM.at\"/>" << endl
        << "            <parameter idref=\"TVM.cg\"/>" << endl
        << "            <parameter idref=\"TVM.gt\"/>" << endl
        << "            <parameter idref=\"TVM.frequencies\"/>" << endl;
    } else if (partitionSubstitutionModel == "SYM") {
        BEAST_xml_code
        << "            <parameter idref=\"SYM.ac\"/>" << endl
        << "            <parameter idref=\"SYM.ag\"/>" << endl
        << "            <parameter idref=\"SYM.at\"/>" << endl
        << "            <parameter idref=\"SYM.cg\"/>" << endl
        << "            <parameter idref=\"SYM.gt\"/>" << endl
        << "            <parameter idref=\"SYM.frequencies\"/>" << endl;
    } else if (partitionSubstitutionModel == "GTR") {
        BEAST_xml_code
        << "            <parameter idref=\"GTR.ac\"/>" << endl
        << "            <parameter idref=\"GTR.ag\"/>" << endl
        << "            <parameter idref=\"GTR.at\"/>" << endl
        << "            <parameter idref=\"GTR.cg\"/>" << endl
        << "            <parameter idref=\"GTR.gt\"/>" << endl
        << "            <parameter idref=\"GTR.frequencies\"/>" << endl;
    }
    if (partitionSiteModel == "IG" || partitionSiteModel == "G") {
        BEAST_xml_code
        << "            <parameter idref=\"alpha\"/>" << endl;
    }
    if (partitionSiteModel == "IG" || partitionSiteModel == "I") {
        BEAST_xml_code
        << "            <parameter idref=\"pInv\"/>" << endl;
    }
    
    if (clockFlavour == "ucln") {
        BEAST_xml_code << endl
        << "            <parameter idref=\"ucld.mean\"/>" << endl
        << "            <parameter idref=\"ucld.stdev\"/>" << endl
        << "            <rateStatistic idref=\"meanRate\"/>" << endl
        << "            <rateStatistic idref=\"coefficientOfVariation\"/>" << endl
        << "            <rateCovarianceStatistic idref=\"covariance\"/>" << endl << endl;
    } else if (clockFlavour == "uced") {
        BEAST_xml_code << endl
        << "            <parameter idref=\"uced.mean\"/>" << endl
        << "            <rateStatistic idref=\"meanRate\"/>" << endl
        << "            <rateStatistic idref=\"coefficientOfVariation\"/>" << endl
        << "            <rateCovarianceStatistic idref=\"covariance\"/>" << endl << endl;
    } else if (clockFlavour == "strict") {
        BEAST_xml_code << endl
        << "            <parameter idref=\"clock.rate\"/>" << endl << endl;
    } else if (clockFlavour == "randlocal") {
        BEAST_xml_code << endl
        << "            <parameter idref=\"clock.rate\"/>" << endl
        << "            <sumStatistic idref=\"rateChanges\"/>" << endl
        << "            <rateStatistic idref=\"meanRate\"/>" << endl
        << "            <rateStatistic idref=\"coefficientOfVariation\"/>" << endl
        << "            <rateCovarianceStatistic idref=\"covariance\"/>" << endl << endl;
    }
    
    BEAST_xml_code
    << "            <treeLikelihood idref=\"treeLikelihood\"/>" << endl;
    
    if (treePrior == "bd" || treePrior == "yule") {
        BEAST_xml_code
        << "            <speciationLikelihood idref=\"speciation\"/>" << endl;
    } else if (treePrior == "concoal") {
        BEAST_xml_code
        << "            <parameter idref=\"constant.popSize\"/>" << endl
        << "            <coalescentLikelihood idref=\"coalescent\"/>" << endl;
    } else if (treePrior == "expcoal") {
        BEAST_xml_code
        << "            <parameter idref=\"exponential.popSize\"/>" << endl
        << "            <parameter idref=\"exponential.growthRate\"/>" << endl
        << "            <coalescentLikelihood idref=\"coalescent\"/>" << endl;
    } else if (treePrior == "logcoal") {
        BEAST_xml_code
        << "            <parameter idref=\"logistic.popSize\"/>" << endl
        << "            <parameter idref=\"logistic.growthRate\"/>" << endl
        << "            <parameter idref=\"logistic.t50\"/>" << endl
        << "            <coalescentLikelihood idref=\"coalescent\"/>" << endl;
    }
    
    BEAST_xml_code
    << "        </log>" << endl << endl;
}

// not used at the moment
/*
void BEASTXML::writeParameterLog (ofstream & BEAST_xml_code, int const& parameterSampling,
    string const& treePrior, string & clockFlavour, vector <string> const& partitionSubstitutionModels,
    vector <string> const& partitionSiteModels)
{
    int counter = 0;
    string prunedFileName = getRootName(XMLOutFileName);
    BEAST_xml_code
    << "<!-- *** PRINT PARAMETERS TO FILE *** -->" << endl
    << "        <log id=\"fileLog\" logEvery=\"" << parameterSampling << "\" fileName=\"" << prunedFileName << ".log\">" << endl
    << "            <posterior idref=\"posterior\"/>" << endl
    << "            <prior idref=\"prior\"/>" << endl
    << "            <likelihood idref=\"likelihood\"/>" << endl
    << "            <parameter idref=\"treeModel.rootHeight\"/>" << endl;
    
    if (treePrior == "bd") {
        BEAST_xml_code
        << "            <parameter idref=\"birthDeath.BminusDRate\"/>" << endl
        << "            <parameter idref=\"birthDeath.DoverB\"/>" << endl << endl;
    } else if (treePrior == "yule") {
        BEAST_xml_code
        << "            <parameter idref=\"yule.birthRate\"/>" << endl << endl;
    } else if (treePrior == "concoal") {
        BEAST_xml_code
        << "            <parameter idref=\"constant.popSize\"/>" << endl << endl;
    } else if (treePrior == "expcoal") {
        BEAST_xml_code
        << "            <parameter idref=\"exponential.popSize\"/>" << endl
        << "            <parameter idref=\"exponential.growthRate\"/>" << endl << endl;
    } else if (treePrior == "logcoal") {
        BEAST_xml_code
        << "            <parameter idref=\"logistic.popSize\"/>" << endl
        << "            <parameter idref=\"logistic.growthRate\"/>" << endl
        << "            <parameter idref=\"logistic.t50\"/>" << endl << endl;
    }
    
    for (vector<string>::const_iterator iter = partitionSubstitutionModels.begin(); iter != partitionSubstitutionModels.end(); iter++) {
        if (partitionSubstitutionModels[counter] == "K80") {
            BEAST_xml_code
            << "            <parameter idref=\"K80.kappa\"/>" << endl
            << "            <parameter idref=\"K80.frequencies\"/>" << endl;
        } else if (partitionSubstitutionModels[counter] == "HKY") {
            BEAST_xml_code
            << "            <parameter idref=\"HKY.kappa\"/>" << endl
            << "            <parameter idref=\"HKY.frequencies\"/>" << endl;
        } else if (partitionSubstitutionModels[counter] == "TrNef") {
            BEAST_xml_code
            << "            <parameter idref=\"TrNef.transversion\"/>" << endl
            << "            <parameter idref=\"TrNef.ag\"/>" << endl
            << "            <parameter idref=\"TrNef.frequencies\"/>" << endl;
        } else if (partitionSubstitutionModels[counter] == "TrN") {
            BEAST_xml_code
            << "            <parameter idref=\"TrN.transversion\"/>" << endl
            << "            <parameter idref=\"TrN.ag\"/>" << endl
            << "            <parameter idref=\"TrN.frequencies\"/>" << endl;
        } else if (partitionSubstitutionModels[counter] == "K3P") {
            BEAST_xml_code
            << "            <parameter idref=\"K3P.purine2pyrimidine\"/>" << endl
            << "            <parameter idref=\"K3P.pyrimidine2purine\"/>" << endl
            << "            <parameter idref=\"K3P.frequencies\"/>" << endl;
        } else if (partitionSubstitutionModels[counter] == "K3Puf") {
            BEAST_xml_code
            << "            <parameter idref=\"K3Puf.purine2pyrimidine\"/>" << endl
            << "            <parameter idref=\"K3Puf.pyrimidine2purine\"/>" << endl
            << "            <parameter idref=\"K3Puf.frequencies\"/>" << endl;
        } else if (partitionSubstitutionModels[counter] == "TIMef") {
            BEAST_xml_code
            << "            <parameter idref=\"TIMef.purine2pyrimidine\"/>" << endl
            << "            <parameter idref=\"TIMef.ag\"/>" << endl
            << "            <parameter idref=\"TIMef.pyrimidine2purine\"/>" << endl
            << "            <parameter idref=\"TIMef.frequencies\"/>" << endl;
        } else if (partitionSubstitutionModels[counter] == "TIM") {
            BEAST_xml_code
            << "            <parameter idref=\"TIM.purine2pyrimidine\"/>" << endl
            << "            <parameter idref=\"TIM.ag\"/>" << endl
            << "            <parameter idref=\"TIM.pyrimidine2purine\"/>" << endl
            << "            <parameter idref=\"TIM.frequencies\"/>" << endl;
        } else if (partitionSubstitutionModels[counter] == "TVMef") {
            BEAST_xml_code
            << "            <parameter idref=\"TVMef.ac\"/>" << endl
            << "            <parameter idref=\"TVMef.at\"/>" << endl
            << "            <parameter idref=\"TVMef.cg\"/>" << endl
            << "            <parameter idref=\"TVMef.gt\"/>" << endl
            << "            <parameter idref=\"TVMef.frequencies\"/>" << endl;
        } else if (partitionSubstitutionModels[counter] == "TVM") {
            BEAST_xml_code
            << "            <parameter idref=\"TVM.ac\"/>" << endl
            << "            <parameter idref=\"TVM.at\"/>" << endl
            << "            <parameter idref=\"TVM.cg\"/>" << endl
            << "            <parameter idref=\"TVM.gt\"/>" << endl
            << "            <parameter idref=\"TVM.frequencies\"/>" << endl;
        } else if (partitionSubstitutionModels[counter] == "SYM") {
            BEAST_xml_code
            << "            <parameter idref=\"SYM.ac\"/>" << endl
            << "            <parameter idref=\"SYM.ag\"/>" << endl
            << "            <parameter idref=\"SYM.at\"/>" << endl
            << "            <parameter idref=\"SYM.cg\"/>" << endl
            << "            <parameter idref=\"SYM.gt\"/>" << endl
            << "            <parameter idref=\"SYM.frequencies\"/>" << endl;
        } else if (partitionSubstitutionModels[counter] == "GTR") {
            BEAST_xml_code
            << "            <parameter idref=\"GTR.ac\"/>" << endl
            << "            <parameter idref=\"GTR.ag\"/>" << endl
            << "            <parameter idref=\"GTR.at\"/>" << endl
            << "            <parameter idref=\"GTR.cg\"/>" << endl
            << "            <parameter idref=\"GTR.gt\"/>" << endl
            << "            <parameter idref=\"GTR.frequencies\"/>" << endl;
        }
        if (partitionSiteModels[counter] == "IG" || partitionSiteModels[counter] == "G") {
            BEAST_xml_code
            << "            <parameter idref=\"alpha\"/>" << endl;
        }
        if (partitionSiteModels[counter] == "IG" || partitionSiteModels[counter] == "I") {
            BEAST_xml_code
            << "            <parameter idref=\"pInv\"/>" << endl;
        }
        counter++;
    }
    
    if (clockFlavour == "ucln") {
        BEAST_xml_code << endl
        << "            <parameter idref=\"ucld.mean\"/>" << endl
        << "            <parameter idref=\"ucld.stdev\"/>" << endl
        << "            <rateStatistic idref=\"meanRate\"/>" << endl
        << "            <rateStatistic idref=\"coefficientOfVariation\"/>" << endl
        << "            <rateCovarianceStatistic idref=\"covariance\"/>" << endl << endl;
    } else if (clockFlavour == "uced") {
        BEAST_xml_code << endl
        << "            <parameter idref=\"uced.mean\"/>" << endl
        << "            <rateStatistic idref=\"meanRate\"/>" << endl
        << "            <rateStatistic idref=\"coefficientOfVariation\"/>" << endl
        << "            <rateCovarianceStatistic idref=\"covariance\"/>" << endl << endl;
    } else if (clockFlavour == "strict") {
        BEAST_xml_code << endl
        << "            <parameter idref=\"clock.rate\"/>" << endl << endl;
    }
    
    BEAST_xml_code
    << "            <treeLikelihood idref=\"treeLikelihood\"/>" << endl;
    
    if (treePrior == "bd" || treePrior == "yule") {
        BEAST_xml_code
        << "            <speciationLikelihood idref=\"speciation\"/>" << endl;
    } else if (treePrior == "concoal") {
        BEAST_xml_code
        << "            <parameter idref=\"constant.popSize\"/>" << endl
        << "            <coalescentLikelihood idref=\"coalescent\"/>" << endl;
    } else if (treePrior == "expcoal") {
        BEAST_xml_code
        << "            <parameter idref=\"exponential.popSize\"/>" << endl
        << "            <parameter idref=\"exponential.growthRate\"/>" << endl
        << "            <coalescentLikelihood idref=\"coalescent\"/>" << endl;
    } else if (treePrior == "logcoal") {
        BEAST_xml_code
        << "            <parameter idref=\"logistic.popSize\"/>" << endl
        << "            <parameter idref=\"logistic.growthRate\"/>" << endl
        << "            <parameter idref=\"logistic.t50\"/>" << endl
        << "            <coalescentLikelihood idref=\"coalescent\"/>" << endl;
    }
    
    BEAST_xml_code
    << "        </log>" << endl << endl;
}
*/

void BEASTXML::writeTreeLogs (ofstream & BEAST_xml_code, string const& clockFlavour,
    bool const& logPhylograms, int const& treeSampling)
{
    string prunedFileName = getRootName(XMLOutFileName);
    if (clockFlavour == "ucln" || clockFlavour == "uced") {
        BEAST_xml_code
        << "<!-- *** TREE LOG FILES *** -->" << endl
        << "        <logTree id=\"treeFileLog\" logEvery=\"" << treeSampling << "\" nexusFormat=\"true\" fileName=\"" << prunedFileName << ".time.trees\" sortTranslationTable=\"true\">" << endl
        << "            <treeModel idref=\"treeModel\"/>" << endl
        << "            <discretizedBranchRates idref=\"branchRates\"/>" << endl
        << "            <posterior idref=\"posterior\"/>" << endl
        << "        </logTree>" << endl;
        if (logPhylograms) {
            BEAST_xml_code
            << "        <logTree id=\"substTreeFileLog\" logEvery=\"" << treeSampling << "\" nexusFormat=\"true\" fileName=\"" << prunedFileName << ".subst.trees\" branchLengths=\"substitutions\">" << endl
            << "            <treeModel idref=\"treeModel\"/>" << endl
            << "            <discretizedBranchRates idref=\"branchRates\"/>" << endl
            << "            <posterior idref=\"posterior\"/>" << endl
            << "        </logTree>" << endl;
        }
    } else if (clockFlavour == "strict") {
        BEAST_xml_code
        << "<!-- *** TREE LOG FILES *** -->" << endl
        << "        <logTree id=\"treeFileLog\" logEvery=\"" << treeSampling << "\" nexusFormat=\"true\" fileName=\"" << prunedFileName << ".time.trees\" sortTranslationTable=\"true\">" << endl
        << "            <treeModel idref=\"treeModel\"/>" << endl
        << "            <strictClockBranchRates idref=\"branchRates\"/>" << endl
        << "            <posterior idref=\"posterior\"/>" << endl
        << "        </logTree>" << endl;
        if (logPhylograms) {
            BEAST_xml_code
            << "        <logTree id=\"substTreeFileLog\" logEvery=\"" << treeSampling << "\" nexusFormat=\"true\" fileName=\"" << prunedFileName << ".subst.trees\" branchLengths=\"substitutions\">" << endl
            << "            <treeModel idref=\"treeModel\"/>" << endl
            << "            <strictClockBranchRates idref=\"branchRates\"/>" << endl
            << "        </logTree>" << endl;
        }
    } else if (clockFlavour == "randlocal") {
        BEAST_xml_code
        << "<!-- *** TREE LOG FILES *** -->" << endl
        << "        <logTree id=\"treeFileLog\" logEvery=\"" << treeSampling << "\" nexusFormat=\"true\" fileName=\"" << prunedFileName << ".time.trees\" sortTranslationTable=\"true\">" << endl
        << "            <treeModel idref=\"treeModel\"/>" << endl
        << "            <randomLocalClockModel idref=\"branchRates\"/>" << endl
        << "            <posterior idref=\"posterior\"/>" << endl
        << "        </logTree>" << endl;
        if (logPhylograms) {
            BEAST_xml_code
            << "        <logTree id=\"substTreeFileLog\" logEvery=\"" << treeSampling << "\" nexusFormat=\"true\" fileName=\"" << prunedFileName << ".subst.trees\" branchLengths=\"substitutions\">" << endl
            << "            <treeModel idref=\"treeModel\"/>" << endl
            << "            <randomLocalClockModel idref=\"branchRates\"/>" << endl
            << "        </logTree>" << endl;
        }
    }
    BEAST_xml_code
    << "    </mcmc>" << endl;
}
