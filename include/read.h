#include "mesh.h"
#ifndef READ_H
#define READ_H

// ---------------------------------------------------------|
// Structure storing the time data and initial conditions   |
// ---------------------------------------------------------|

struct timeStruct{

    // Saving frequency and total number of time steps

    double dt;
    int nSave;
    int nSteps;

    // Initial solution in displacement

    darray u0;
};

// --------------------------------------------------------------|
// Structure storing temporary parameters from the input files   |
// --------------------------------------------------------------|

struct readStruct{

    // Young modulus and Poisson ratio of each type of element
    
    dvector emptyEv;
    std::vector<bool> empty;
    std::vector<dvector> Ev;
    std::vector<dvector> EvS;

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

    // Axis of the applied stress
    
    std::vector<std::pair<int,double>> axial;

    // Type and axis of the applied stress
    
    std::string type;
    std::string load;

    // Lis of neightbour elements

    std::vector<ivector> neighbour;

    // list of coupled and delta displacement (face,axis)

    std::vector<std::pair<int,int>> coupled;
    std::vector<std::pair<int,int>> lockBot;
    std::vector<std::pair<int,int>> lockTop;
    std::vector<std::pair<int,int>> uniform;
    std::vector<std::pair<int,int>> deltaZero;
};

// -------------------------------------------|
// Functions available in the read.cpp file   |
// -------------------------------------------|

dvector tovec(std::string input);
void read(std::string path[2],dataStruct &data,timeStruct &time);

// Functions to set the boundary conditions

void neumann(readStruct &read,dataStruct &data);
void dirichlet(readStruct &read,dataStruct &data);

// Functions to read the different input files

void setSurface(readStruct &read,dataStruct &data);
void readInput(readStruct &read,dataStruct &data,std::string path);
std::unordered_set<int> locSpecies(readStruct &read, dvector coord);
void readMeshSize(readStruct &read,dataStruct &data,std::string path);
void readSpecies(readStruct &read,dataStruct &data,std::string path);

#endif