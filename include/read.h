#include "mesh.h"
#ifndef READ_H
#define READ_H

struct readStruct{

    // Dh = stiffness tensor of holes
    // Db = stiffness tensor of the bulk
    // threshold = maximum fraction of a hole
    // boundary = [dimension] [pair of type and value]

    matrix Dh,Db;
    double threshold;
    std::vector<spair> boundary;
};

dvector tovec(std::string input);
void readMesh(readStruct &read,meshStruct &mesh,std::string path);
void readInput(readStruct &read,meshStruct &mesh,std::string path);
meshStruct read(std::string inputPath,std::string meshPath);

#endif