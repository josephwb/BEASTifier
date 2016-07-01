#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
#include <cstdlib>

using namespace std;

#include "General.h"

extern bool DEBUG;

// general processing and error-checking functions

bool checkValidFloat (string stringToCheck) {
    float tempFloat = 0.0;
    bool validInput = true;
    
    istringstream tempStream(stringToCheck);
    tempStream >> tempFloat;
    
    if (tempStream.fail()) {
        cout << endl << "Invalid input! Looking for a float argument. Exiting." << endl << endl;
        exit(1);
    }
    return validInput;
}

bool checkValidBoolInput (string queryString) {
    int userInput = 0;
    bool validInput = false;
    cout << queryString;
    validInput = false;
    while (!validInput) {
        cout << endl << "Enter (0) to accept this, or (1) to change: ";
        cin >> userInput;
        if (cin.fail()) {
            cin.clear(); 
            cin.ignore(200, '\n');
            cout << endl << "Invalid input! Must be an integer. Try again." << endl;
            cout << queryString;
        } else if (userInput < 0 || userInput > 1) {
            cin.clear(); 
            cin.ignore(200, '\n');
            cout << endl << "Invalid input! Boolean value - must be 0 or 1. Try again." << endl;
            cout << queryString;
        } else {
            validInput = true;
            cin.clear();             // Get rid of any extra characters accidentally entered
            cin.ignore(200, '\n');
        }
    }
    return userInput;
}

bool checkValidInputFile (string fileName) {
    bool validInput = false;
    ifstream tempStream;
    
    tempStream.open(fileName.c_str());
    if (tempStream.fail()) {
        ofstream errorReport("Error.BEASTifier.txt");
        errorReport << "BEASTifier  failed." << endl << "Error: unable to open file '";
        errorReport << fileName << "'" << endl;
        errorReport.close();
        cerr << endl << "BEASTifier failed." << endl << "Error: unable to open file '";
        cerr << fileName << "'. You fucked up, yo. Exiting." <<  endl << endl;
        exit(1);
    } else {
        if (DEBUG) {cout << "Successfully opened file '" << fileName << "'." <<  endl;}
        tempStream.close();
        tempStream.clear();
    }
    return validInput;
}

bool checkValidOutputFile (string & outputFileName, bool const& overwrite) {
    bool testOutBool = true;
    bool fileNameAcceptable = false;
    bool keepFileName = false;
    
// First, check if file already exists, so overwriting can be prevented
    fstream testIn;
    while (!fileNameAcceptable) {
        testIn.open(outputFileName.c_str());
        if (!testIn) {
            testIn.close();
            fileNameAcceptable = true;
        } else if (overwrite) {
            fileNameAcceptable = true;
        } else {
            if (overwrite) {
                fileNameAcceptable = true;
            } else {
                testIn.close();
                cout << endl << "File exists!  Change name (0) or overwrite (1)? ";
                cin >> keepFileName;
                if (!keepFileName) {
                    cout << "Enter new output file name: ";
                    cin >> outputFileName;
                } else {
                    cout << "Overwriting existing file '" << outputFileName << "'." << endl;
                    fileNameAcceptable = true;
                }
            }
        }
    }
    
    ofstream outFile;
    outFile.open(outputFileName.c_str());
    
    if (outFile.fail()) {
        ofstream errorReport("Error.BEASTifier.txt");
        errorReport << "BEASTifier analysis failed." << endl << "Error: unable to open file '";
        errorReport << outputFileName << "'" << endl;
        errorReport.close();

        cerr << endl << "BEASTifier analysis failed." << endl << "Error: unable to open file '";
        cerr << outputFileName << "'" <<  endl;
        testOutBool = false;
        exit(1);
    } else {
        outFile.close();
        outFile.clear();
    }
    return testOutBool;
}

string parseString (string stringToParse, int stringPosition) {
    vector <string> tempVector;
    istringstream tempStream(stringToParse);
    string tempString;
    while (tempStream >> tempString) {
        tempVector.push_back(tempString);
    }
    return tempVector[stringPosition];
}

bool checkStringValue (string stringToParse, string stringToMatch, int stringPosition) {
// Performs case-insenstive string match test
    string testString = parseString(stringToParse, stringPosition);
    if (testString.size() != stringToMatch.size()) {
        return false;
    }
    if (testString == stringToMatch) {
        return true;
    }
    
    for (size_t i = 0; i < testString.size(); ++i) {
        if (toupper(testString[i]) == toupper(stringToMatch[i])) {
            continue;
        } else {
            return false;
        }
    }
    return true;
}

bool checkCharValue (char const& charInput, char const& charToMatch) {
    if (toupper(charInput) == toupper(charToMatch)) {
        return true;;
    } else {
        return false;
    }
}

bool checkWhiteSpaceOnly (string stringToParse) {
    bool whiteSpaceOnly = true;
    vector<string> tempVector;
    istringstream tempStream(stringToParse);
    string tempString;
    while (tempStream >> tempString) {
        if (tempString != "    " && tempString != " ") {
            whiteSpaceOnly = false;
        }
    }
    return whiteSpaceOnly;
}

bool checkCommentLineNexus (string stringToParse) {
    bool commentLine = false;
    char firstCharacter = stringToParse[0];
    if (firstCharacter == '[') {
        commentLine = true;
    }
    return commentLine;
}

int convertStringtoInt (string stringToConvert) {
    int tempInt = 0;
    istringstream tempStream(stringToConvert);
    tempStream >> tempInt;
    
    return tempInt;
}

string convertIntToString (int intToConvert) {
    string tempString;
    stringstream tempStream;
    tempStream << intToConvert;
    tempString = tempStream.str();
    
    return tempString;
}

string removeStringSuffix (string stringToParse, char suffixToRemove, bool & suffixEncountered) {
    string temp;
    vector<char> tempVector;
    int charCounter = 0;
    int suffixStart = 0;
    suffixEncountered = false;
    
    for (string::iterator charIter = stringToParse.begin(); charIter < stringToParse.end(); charIter++ ) {
        tempVector.push_back(*charIter);
        if (*charIter == suffixToRemove) {
            suffixStart = charCounter;    // will record *last* instance of suffix delimiter
            suffixEncountered = true;
        }
        charCounter++;
    }
    if (suffixEncountered) {
        for (int charIter = 0; charIter < suffixStart; charIter++) {
            temp += tempVector[charIter];
        }
    }
    if (!suffixEncountered) {
        temp = stringToParse;
    }
    return temp;
}

string removeStringPrefix(string stringToParse, char characterToRemove) {
    string tempString;
    vector<char> tempVector;
    int charCounter = 0;
    int characterPosition = 0;
    bool characterEncountered = false;
    
    for (string::iterator charIter = stringToParse.begin(); charIter < stringToParse.end(); charIter++ ) {
        tempVector.push_back(*charIter);
        if (!characterEncountered) {    // make sure first occurrence alone is recorded
            if (*charIter == characterToRemove) {
                characterPosition = charCounter;
                characterEncountered = true;
            }
        }
        charCounter++;
    }
    if (characterEncountered) {
        for (int charIter = characterPosition + 1; charIter < charCounter; charIter++) {
            tempString += tempVector[charIter];
        }
    }
    if (!characterEncountered) {
        tempString = stringToParse;
    }
    return tempString;
}

string getLongestName (vector <string> const& elements) {
    string longestName;
    
// Initialize with starting value
    longestName = elements[0];
    
    for (vector <string>::const_iterator elementIter = elements.begin(); elementIter < elements.end(); elementIter++) {
        string currentNameString = *elementIter;
        if (currentNameString.size() > longestName.size()) {
            longestName = currentNameString;
        }
    }
    return longestName;
}

void printFormattingSpaces (string const& longestString, string const& currentString) {
    if (currentString.size() < longestString.size()) {
        string::size_type tempDiffSize;
        tempDiffSize = longestString.size() - currentString.size();
        for (string::size_type iterSpaces = 0; iterSpaces < tempDiffSize; iterSpaces++) {
            cout << " ";
        }
    }
}

string getStringElement (string const& stringToBreak, char const& delimiter, int const& elementPosition) { // of form: b1d0a0.1n25.GTR.2.NEX; very strict
    string desired = stringToBreak;
    
    for (int i = 0; i < elementPosition-1; i++) {
        desired = removeStringPrefix(desired, delimiter);
    }
    
    bool working = true;
    while (working) {
        desired = removeStringSuffix(desired, delimiter, working);
    }
    
    return desired;
}

string getRootName (string const& stringToBreak) {
// form: b_1_d_0.5_a_0.9_n_100_sim_JC_rep_4.NEX - 12 root elements
    string desired = stringToBreak;
    bool done = false;
    desired = removeStringSuffix(desired, '.', done); // removes last instance of '.'
    return desired;
}

vector <string> tokenizeString (string const& stringToParse) {
    vector <string> result;
    istringstream tempStream(stringToParse);
    string tempString;
    while (tempStream >> tempString) {
        result.push_back(tempString);
    }
    return result;
}

vector <string> readFileList (string const& fileName) {
    vector <string> alignments;
    ifstream inputNames;
    inputNames.open(fileName.c_str());
    
    string line;
    
    while (getline(inputNames,line)) {
        if (!checkWhiteSpaceOnly(line)) {
            alignments.push_back(line);
        }
    }
    return (alignments);
}
