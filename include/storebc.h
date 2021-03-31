#include "mesh.h"
#ifndef STOREBC_H
#define STOREBC_H

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

    // Face index, axis of the applied stress and its value

    ivector axis[2];
    double Fval;

    // Type, value of the applied stress and ndex of the free surface
    
    std::string type;
    std::string load;
    ivector free;

    // Lis of neightbour elements

    std::vector<ivector> neighbour;

    // list of coupled and delta displacement (face,axis)

    std::vector<std::pair<int,int>> coupled;
    std::vector<std::pair<int,int>> lockBot;
    std::vector<std::pair<int,int>> lockTop;
    std::vector<std::pair<int,int>> uniform;
    std::vector<std::pair<int,double>> axial;
    std::vector<std::pair<int,int>> deltaZero;
};

// ----------------------------------|
// Functions available in the file   |
// ----------------------------------|

void neumann(readStruct &read,dataStruct &data);
void dirichlet(readStruct &read,dataStruct &data);

#endif