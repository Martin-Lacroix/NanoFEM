#include "mesh.h"
#ifndef READ_H
#define READ_H

// --------------------------------------------------------------|
// Structure storing temporary parameters from the input files   |
// --------------------------------------------------------------|

struct readStruct{

    // Young modulus and Poisson ratio of each type of element

    dvector E;
    dvector v;

    // Number of elements and size of the domain in each dimension

    dvector dSize;
    ivector dLen;

    // BC and keeps track of the row of each node in each dimension in perNode

    std::vector<sdpair> boundary;
    std::vector<ivector> row;
};

// -------------------------------------------|
// Functions available in the read.cpp file   |
// -------------------------------------------|

dvector tovec(std::string input);
void setBC(readStruct &read,meshStruct &mesh,ivector loop);
meshStruct read(std::string inputPath,std::string meshPath);
void readMesh(readStruct &read,meshStruct &mesh,std::string path);
void readInput(readStruct &read,meshStruct &mesh,std::string path);

#endif