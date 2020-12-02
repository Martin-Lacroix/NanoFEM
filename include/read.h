#include "mesh.h"
#ifndef READ_H
#define READ_H

struct readStruct{

    // Dirichlet = displacement perpendicular to the faces
    // Domain = maximum of x, y and z coordinates
    // Zero = minimum of w x, y and z coordinates

    dvector dirichlet;
    dvector domain;
    dvector zero;
};

void cleanFace(meshStruct &mesh);
dvector tovec(std::string input);
void readMesh(readStruct &read,meshStruct &mesh,std::string path);
void readInput(readStruct &read,meshStruct &mesh,std::string path);
meshStruct read(std::string inputPath,std::string meshPath);

#endif