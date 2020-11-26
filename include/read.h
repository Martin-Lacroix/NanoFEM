#include "mesh.h"
#include "linalg.h"
#include "stdafx.h"
#ifndef READ_H
#define READ_H

typedef alglib::real_2d_array matrix;
typedef alglib::real_1d_array darray;
typedef alglib::integer_1d_array iarray;

struct readStruct{

    std::vector<iarray> fixed;
    std::vector<darray> force;
};

matrix stiffness(double E,double v);
readStruct readParam(otherStruct &param);
void readAll(meshStruct &mesh, otherStruct &param);

#endif