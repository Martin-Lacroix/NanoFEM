#include "mesh.h"
#ifndef READ_H
#define READ_H

struct readStruct{

    // E and v are Young modulus and Poisson ratio
    // dom is the space domain of the simulation
    // boundary is the type and value of the BC in each dimension.

    dvector E;
    dvector v;
    dvector dSize;
    ivector dLen;
    std::vector<spair> boundary;

};

dvector tovec(std::string input);
void readMesh(readStruct &read,meshStruct &mesh,std::string path);
void readInput(readStruct &read,meshStruct &mesh,std::string path);
void setBC(readStruct &read,meshStruct &mesh,std::vector<ivector> &row,ivector loop);
meshStruct read(std::string inputPath,std::string meshPath);

#endif