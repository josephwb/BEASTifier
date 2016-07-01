#ifndef _USER_INTERFACE_H_
#define _USER_INTERFACE_H_

void printProgramInfo();
void processCommandLineArguments (int argc, char *argv[], vector <string> & listFileNames,
    AnalysisSettings & ASet);
void usage ();
void readConfigFile (string const& fileName, AnalysisSettings & ASet,
    vector <string> & listFileNames);
bool checkComment (string const& val);

#endif /* _USER_INTERFACE_H_ */
