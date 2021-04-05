#include <unordered_set>
#include "storebc.h"
#ifndef PARSER_H
#define PARSER_H

// ----------------------------------|
// Functions available in the file   |
// ----------------------------------|

dvector tovec(std::string input);
void read(std::string path,dataStruct &data);

// Functions to read the different input files

void setSurface(readStruct &read,dataStruct &data);
std::unordered_set<int> locSpecies(readStruct &read,dvector coord);
void readSpecies(readStruct &read,dataStruct &data,std::string path);
void readMeshSize(readStruct &read,dataStruct &data,std::string path);
std::string readInput(readStruct &read,dataStruct &data,std::string path);

#endif