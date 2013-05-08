#ifndef _GENERAL_H_
#define _GENERAL_H_

// General use functions
bool checkValidFloat (string stringToCheck);
bool checkValidBoolInput (string queryString);
bool checkValidInputFile (string fileName);
bool checkValidOutputFile (string & outputFileName, bool const& overwrite);
string parseString (string stringToParse, int stringPosition);
bool checkStringValue (string stringToParse, string stringToMatch, int stringPosition);
bool checkCharValue (char const& charInput, char const& charToMatch);
bool checkWhiteSpaceOnly (string stringToParse);
bool checkCommentLineNexus (string stringToParse);
int convertStringtoInt (string stringToConvert);
string convertIntToString (int intToConvert);
string removeStringSuffix (string stringToParse, char suffixToRemove, bool & suffixEncountered);
string removeStringPrefix (string stringToParse, char characterToRemove);
string getLongestName (vector<string> const& elements);
void printFormattingSpaces (string const& longestString, string const& currentString);
string getStringElement (string const& stringToBreak, char const& delimiter, int const& elementPosition);
string getRootName (string const& stringToBreak);
vector <string> tokenizeString (string const& stringToParse);
vector <string> readFileList (string const& fileName);

#endif /* _GENERAL_H_ */
