#include "mesh.h"
#ifndef READ_H
#define READ_H

struct readStruct{

    // Dh = stiffness tensor of holes
    // Db = stiffness tensor of bulk material
    // Threshold = maximum fraction of a hole
    // Boundary = types of boundary conditions

    matrix Dh,Db;
    double threshold;
    std::vector<std::string> boundary;
};

void cleanFace(meshStruct &mesh);
dvector tovec(std::string input);
void readMesh(readStruct &read,meshStruct &mesh,std::string path);
void readInput(readStruct &read,meshStruct &mesh,std::string path);
meshStruct read(std::string inputPath,std::string meshPath);

#endif