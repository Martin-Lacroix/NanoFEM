#include <unordered_set>
#include "writer.h"
#ifndef PARSER_H
#define PARSER_H

// --------------------------------------------------------------|
// Structure storing temporary parameters from the input files   |
// --------------------------------------------------------------|

struct readStruct{

    // Young modulus and Poisson ratio of each type of element
    
    array3d emptyLmR;
    std::vector<bool> empty;
    std::vector<array3d> LmR;
    std::vector<array3d> LmS;

    // Number of elements and size of the domain in each dimension

    dvector eSize;
    dvector dSize;
    ivector dLen;
    dvector zero;
    double cropZ;
    double Lc;

    // Neumann boundary conditions

    std::string axis;
    double Fval;

    // Stores the boundary condition parameters

    std::vector<ivector> neighbour;
    std::string deformation;
    std::string transverse;
    ivector opposite;
    ivector free;

    // list of coupled and delta displacement (face,axis)

    std::vector<int> delta;
    std::vector<std::pair<int,int>> coupled;
    std::vector<std::pair<int,int>> lockBot;
    std::vector<std::pair<int,int>> uniform;
};

// ----------------------------------|
// Functions available in the file   |
// ----------------------------------|

dvector tovec(std::string input);
readStruct read(std::string path,dataStruct &data);

// Functions to read the different input files

void setSurface(readStruct &read,dataStruct &data);
std::unordered_set<int> locSpecies(readStruct &read,dvector coord);
void readSpecies(readStruct &read,dataStruct &data,std::string path);
void readMeshSize(readStruct &read,dataStruct &data,std::string path);
std::string readInput(readStruct &read,dataStruct &data,std::string path);

#endif