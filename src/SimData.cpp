#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>

using namespace std;

#include "General.h"
#include "SimData.h"

extern bool DEBUG;

SimData::SimData (string const& fileName)
: numTaxa(0), numChar(0), interleavedData(false), starterTreePresent(false)
{
    seqFileName = fileName;
    checkValidInputFile(fileName);
    setNumTaxaChar(seqFileName, numTaxa, numChar, interleavedData);
    taxaAlignment = collectTaxaAlignment(seqFileName, numTaxa, numChar, interleavedData);
    
    root = setRootName(seqFileName);
    simModel = getStringElement(seqFileName, '_', 10);
    treeFileName = getTreeName(root);
    starterTree = collectStartingTreePhylip(treeFileName, starterTreePresent);
    
}

void SimData::setNumTaxaChar (string & seqFileName, int & numTaxa, int & numChar,
    bool & interleavedData)
{
    ifstream inputUserFile;
    bool commentLine = false;
    bool whiteSpaceOnly = false;
    bool numTaxaEncountered = false;
    bool numCharEncountered = false;
    bool equalSignEncountered = false;
    bool semicolonEncountered = false;
    bool matrixEncountered = false; // know to stop looking
    
    inputUserFile.open(seqFileName.c_str());
    string line;
    
// Looking for pattern like 'dimensions ntax=53 nchar=16620;'
//     - can be in either order, but must be stated on same line (for now)
//     - no spaces allowed next to equal sign (for now)

// Looking for pattern like 'Format datatype=dna [gap=-] [missing=?] {[interleave=yes] or [interleave]};'
//     - no spaces allowed next to equal sign (for now)

    while (getline(inputUserFile,line) && !matrixEncountered) {
        int stringPosition = 0;
        commentLine = checkCommentLineNexus(line);
        whiteSpaceOnly = checkWhiteSpaceOnly(line);
        if (line.empty() || commentLine || whiteSpaceOnly) {
            continue;
        } else {
            if (checkStringValue(line, "matrix", stringPosition)) {    // Done - won't find the information further down; really only used for 'interleave'
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
                    if (checkStringValue(tempString, "nchar", 0)) {
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
                    if (checkStringValue(tempString, "interleave", 0)) {                                          // OR a space follows interleave declaration i.e. 'interleave [=];'
                        interleavedData = true;
                        semicolonEncountered = true;
                        cout << "Data are in interleaved format." << endl;
                        continue;
                    }
                }
            }
        }
    }
    //cout << endl;
    inputUserFile.close();
}

vector < vector <string> > SimData::collectTaxaAlignment (string & seqFileName, int const& numTaxa,
    int const& numChar, bool const& interleavedData)
{
// PLEASE NOTE: search strategy below uses very strict format assumptions - that which is exported by PAUP*
//     - do not be surprised if this fucks up - it is probably a simple rearrangement of terms

    vector < vector <string> > taxaAlignment;
    ifstream inputAlignment;
    vector <string> tempStringVector;
    bool commentLine = false;
    bool whiteSpaceOnly = false;
    bool matrixEncountered = false;
    bool allCharacterRead = false;
    int numCharRead = 0;
    
    inputAlignment.open(seqFileName.c_str());
    string line;
    
    if (!interleavedData) {
        while (!matrixEncountered) { // Ignore lines until 'matrix' is encountered
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
                        tempStringVector.push_back(taxonName);        // Taxon name
                        tempStringVector.push_back(taxonSequence);    // sequence
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

string SimData::setRootName (string const& stringToBreak) {
// form: b_1_d_0.5_a_0.9_n_100_sim_JC_rep_4.NEX - 12 root elements
    string desired = stringToBreak;
    bool done = false;
    desired = removeStringSuffix(desired, '.', done); // removes last instance of '.'
    return desired;
}

string SimData::getTreeName(string const& rootName) {
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

string SimData::collectStartingTreePhylip (string & fileName, bool & starterTreePresent) {
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
