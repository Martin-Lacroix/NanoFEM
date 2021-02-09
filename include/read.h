#include "mesh.h"
#ifndef READ_H
#define READ_H

// --------------------------------------------------------------|
// Structure storing temporary parameters from the input files   |
// --------------------------------------------------------------|

struct readStruct{

    // Young modulus and Poisson ratio of each type of element
    
    dvector emptyEv;
    std::vector<dvector> Ev;

    // Number of elements and size of the domain in each dimension

    dvector eSize;
    dvector dSize;
    ivector dLen;
    dvector zero;

    // BC and keeps track of the row of each node in each dimension in perNode

    std::vector<ivector> neighbour;
    std::vector<sdpair> boundary;
    std::vector<ivector> row;
};

// -------------------------------------------|
// Functions available in the read.cpp file   |
// -------------------------------------------|

dvector tovec(std::string input);
meshStruct read(std::string inputPath,std::string meshPath);

// Functions to read the different input files

void setBC(readStruct &read,meshStruct &mesh);
void readInput(readStruct &read,meshStruct &mesh,std::string path);
std::unordered_set<int> locSpecies(readStruct &read, dvector coord);
void readMeshSize(readStruct &read,meshStruct &mesh,std::string path);
void readSpecies(readStruct &read,meshStruct &mesh,std::string path);

#endif