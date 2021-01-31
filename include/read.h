#include "mesh.h"
#ifndef READ_H
#define READ_H

struct readStruct{

    // E and v are Young modulus and Poisson ratio
    // boundary is the type and value of the BC in each dimension.

    dvector E;
    dvector v;
    std::vector<spair> boundary;

};

dvector tovec(std::string input);
void readMesh(readStruct &read,meshStruct &mesh,std::string path);
void readInput(readStruct &read,meshStruct &mesh,std::string path);
meshStruct read(std::string inputPath,std::string meshPath);

#endif