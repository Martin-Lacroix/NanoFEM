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
    double cropZ;

    // Face index, axis of the applied stress and its value

    ivector axis[2];
    double value;

    // Index of the clamped and flat surfaces

    ivector flat;
    ivector clamped;

    // Type and axis of the applied stress
    
    std::string type;
    std::string load;

    // Lis of neightbour elements

    std::vector<ivector> neighbour;

    // list of coupled displacement as (face,axis)

    std::vector<std::pair<int,int>> coupled;
};

// -------------------------------------------|
// Functions available in the read.cpp file   |
// -------------------------------------------|

dvector tovec(std::string input);
meshStruct read(std::string inputPath,std::string meshPath);

// Functions to set the boundary conditions

void neumann(readStruct &read,meshStruct &mesh);
void dirShear(readStruct &read,meshStruct &mesh);
void dirTensile(readStruct &read,meshStruct &mesh);


void dirichlet(readStruct &read,meshStruct &mesh);

// Functions to read the different input files

void readInput(readStruct &read,meshStruct &mesh,std::string path);
std::unordered_set<int> locSpecies(readStruct &read, dvector coord);
void readMeshSize(readStruct &read,meshStruct &mesh,std::string path);
void readSpecies(readStruct &read,meshStruct &mesh,std::string path);

#endif